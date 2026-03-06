from __future__ import annotations

import itertools
import sys
from collections import namedtuple
from collections.abc import Callable
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import NamedTuple, TextIO, cast

import anndata as ad
import numpy as np
import numpy.typing as npt
import pandas as pd
import scanpy as sc
import sklearn
import sklearn.neighbors
import sklearn.preprocessing
from loguru import logger
from scipy import sparse
from zfpkm import zFPKM, zfpkm_plot

from como.data_types import FilteringTechnique, LogLevel, RNAType
from como.migrations import gene_info_migrations
from como.pipelines.identifier import contains_identical_gene_types, determine_gene_type
from como.project import Config
from como.utils import read_file, set_up_logging


class _FilteringOptions(NamedTuple):
    replicate_ratio: float
    batch_ratio: float
    cut_off: float
    high_replicate_ratio: float
    high_batch_ratio: float


class LayoutMethod(Enum):
    """RNA sequencing layout method."""

    paired_end = "paired-end"
    single_end = "single-end"


@dataclass(slots=True)
class _StudyMetrics:
    study: str
    num_samples: int
    count_matrix: pd.DataFrame | sc.AnnData
    eff_length: pd.DataFrame  # The effective length of the fragment from, e.g., Salmon Quantification
    sample_names: list[str]
    layout: list[LayoutMethod]
    entrez_gene_ids: npt.NDArray[np.integer]
    gene_sizes: npt.NDArray[np.integer]
    __normalization_matrix: pd.DataFrame = field(default_factory=pd.DataFrame)
    __z_score_matrix: pd.DataFrame | sc.AnnData | None = field(default=None)
    __high_confidence_entrez_gene_ids: list[str] = field(default_factory=list)

    def __post_init__(self):
        if not self.layout:
            raise ValueError("Layout list cannot be empty")
        if not self.sample_names:
            raise ValueError("Sample names list cannot be empty")

        for layout in self.layout:
            if layout not in LayoutMethod:
                raise ValueError(f"Layout must be 'paired-end' or 'single-end'; got: {layout}")

    @property
    def normalization_matrix(self) -> pd.DataFrame:
        return self.__normalization_matrix

    @normalization_matrix.setter
    def normalization_matrix(self, value: pd.DataFrame) -> None:
        self.__normalization_matrix = value

    @property
    def z_score_matrix(self) -> pd.DataFrame | sc.AnnData | None:
        return self.__z_score_matrix

    @z_score_matrix.setter
    def z_score_matrix(self, value: pd.DataFrame | sc.AnnData) -> None:
        self.__z_score_matrix = value

    @property
    def high_confidence_entrez_gene_ids(self) -> list[str]:
        return self.__high_confidence_entrez_gene_ids

    @high_confidence_entrez_gene_ids.setter
    def high_confidence_entrez_gene_ids(self, values: list[str]) -> None:
        self.__high_confidence_entrez_gene_ids = values


class _ZFPKMResult(NamedTuple):
    zfpkm: pd.Series
    density: Density
    mu: float
    std_dev: float
    fpkm_at_mu: float


class _ReadMatrixResults(NamedTuple):
    metrics: dict[str, _StudyMetrics]
    entrez_gene_ids: list[str]


Density = namedtuple("Density", ["x", "y"])
NamedMetrics = dict[str, _StudyMetrics]


def k_over_a(k: int, a: float) -> Callable[[npt.NDArray], bool]:
    """Filter rows of an array based on the sum of elements being greater than or equal to A at least k times.

    This code is based on the `kOverA` function found in R's `genefilter` package
    https://www.rdocumentation.org/packages/genefilter/versions/1.54.2/topics/kOverA

    :param k: The minimum number of times the sum of elements must be greater than or equal to A.
    :param a: The value to compare the sum of elements to.
    :return: A function that accepts a NumPy array to perform the actual filtering
    """

    def filter_func(row: npt.NDArray) -> bool:
        return bool(np.sum(row >= a) >= k)

    return filter_func


def genefilter(data: pd.DataFrame | npt.NDArray, filter_func: Callable[[npt.NDArray], bool]) -> npt.NDArray:
    """Apply a filter function to the rows of the data and return the filtered array.

    This code is based on the `genefilter` function found in R's `genefilter` package: https://www.rdocumentation.org/packages/genefilter/versions/1.54.2/topics/genefilter

    Arg:
        data: The data to filter, either a Pandas DataFrame or a NumPy array.
        filter_func: THe function to filter the data by

    Returns:
        A NumPy array of the filtered data.
    """
    if not isinstance(data, (pd.DataFrame, np.ndarray)):
        raise TypeError(f"Unsupported data type. Must be a Pandas DataFrame or a NumPy array, got '{type(data)}'")

    return (
        data.apply(filter_func, axis=1).to_numpy(copy=False)
        if isinstance(data, pd.DataFrame)
        else np.apply_along_axis(filter_func, axis=1, arr=data)
    )


def _build_matrix_results(
    matrix: pd.DataFrame | sc.AnnData,
    *,
    gene_info: pd.DataFrame,
    metadata_df: pd.DataFrame,
    normalize_df: pd.DataFrame,
    taxon: int,
    normalize_df_is_fragment_length: bool = True,
) -> tuple[NamedMetrics, list[int]]:
    """Read the counts matrix and returns the results.

    :param matrix: The gene counts matrix to process
    :param metadata_df: The configuration dataframe related to the current context
    :param taxon: The NCBI Taxon ID
    :returns: A dataclass `ReadMatrixResults`
    """
    gene_info = gene_info_migrations(gene_info)
    gene_info = gene_info[~gene_info["entrez_gene_id"].isna()]
    gene_info = gene_info[gene_info["entrez_gene_id"] != "-"]
    gene_info["entrez_gene_id"] = gene_info["entrez_gene_id"].astype(int)

    if not normalize_df.empty:  # This is likely empty when using single-cell data
        normalize_merge_on = normalize_df.columns.intersection(gene_info.columns).to_list()
        if not normalize_merge_on:
            raise ValueError(
                f"No common columns to merge on between effective length dataframe and gene info dataframe. "
                f"Effective length columns: {normalize_df.columns.tolist()}, gene info columns: {gene_info.columns.tolist()}"
            )
        normalize_df = normalize_df.merge(gene_info, on=normalize_merge_on, how="left")
        normalize_df = normalize_df[~normalize_df["entrez_gene_id"].isna()]

    if isinstance(matrix, pd.DataFrame):
        matrix_merge_on = matrix.columns.intersection(gene_info.columns).to_list()
        matrix = matrix.merge(gene_info, on=matrix_merge_on, how="inner")
    elif isinstance(matrix, sc.AnnData):
        if not isinstance(matrix.var, pd.DataFrame):
            raise TypeError(f"Expected matrix.var object to be 'pd.DataFrame', got '{type(matrix.var)}'")

        gene_info = gene_info.sort_values(["entrez_gene_id", "size"], ascending=[True, True]).drop_duplicates(
            subset=["entrez_gene_id"], keep="first"
        )
        types: dict[str, str] = determine_gene_type(matrix.var_names)
        if contains_identical_gene_types(types):
            # set matrix index name to gene type
            matrix.var.index.name = "gene_symbol"

        matrix.var = matrix.var.reset_index(drop=False).merge(gene_info, on="gene_symbol", how="left")
        matrix = matrix[:, matrix.var["entrez_gene_id"].notna()]
        matrix.var = matrix.var.astype(
            {
                "entrez_gene_id": int,
                "taxon_id": int,
                "size": int,
            }
        )

        # new_matrix = matrix[:, new_var["original_index"]].copy()
        # new_matrix.var = new_var
        # new_matrix.var = new_matrix.var.drop(columns=["original_index"])
        # new_matrix.var.reset_index(drop=True)
        # matrix = new_matrix
        # non_duplicates = ~matrix.var.duplicated(subset=matrix.var.columns, keep="first")
        # matrix = matrix[:, non_duplicates].copy()

    metrics: NamedMetrics = {}
    for study in metadata_df["study"].unique():
        study_sample_names: list[str] = metadata_df[metadata_df["study"] == study]["sample_name"].tolist()
        layouts: list[str] = metadata_df[metadata_df["study"] == study]["layout"].tolist()

        if isinstance(matrix, pd.DataFrame):
            subset = matrix.set_index(keys=["entrez_gene_id"], drop=True)
            subset = subset[subset.columns.intersection(study_sample_names)]
            # subset = subset.groupby(subset.index).mean()
            entrez_gene_ids = subset.index.to_numpy(copy=False)

        elif isinstance(matrix, sc.AnnData):
            # matrix.var = matrix.var.set_index(keys=["entrez_gene_id"], drop=True)
            subset = matrix[matrix.obs_names.intersection(study_sample_names)]
            entrez_gene_ids = subset.var["entrez_gene_id"].to_numpy(dtype=int)
        else:
            raise TypeError(f"Matrix must be a pandas DataFrame or scanpy AnnData object, got: '{type(matrix)}'.")

        normalize_df = normalize_df.drop(columns=["ensembl_gene_id", "gene_symbol"], errors="ignore")
        if normalize_df.index.name != "entrez_gene_id" and "entrez_gene_id" in normalize_df.columns:
            normalize_df = normalize_df.set_index("entrez_gene_id", drop=True)
        normalize_df.index = normalize_df.index.astype(int)

        # gene_sizes = gene_info["size"].to_numpy(dtype=int, copy=False)
        metrics[study] = _StudyMetrics(
            count_matrix=subset,
            eff_length=normalize_df[study_sample_names] if not normalize_df.empty else normalize_df,
            sample_names=study_sample_names,
            layout=[LayoutMethod(layout) for layout in layouts],
            num_samples=len(study_sample_names),
            entrez_gene_ids=entrez_gene_ids,
            gene_sizes=gene_info.loc[gene_info["entrez_gene_id"].isin(entrez_gene_ids), "size"].to_numpy(
                dtype=int, copy=False
            ),
            study=study,
        )

    return metrics, gene_info["entrez_gene_id"].astype(int).tolist()


def _calculate_fpkm(metrics: NamedMetrics) -> NamedMetrics:
    """Calculate the Fragments Per Kilobase of transcript per Million mapped reads (FPKM).

    Args:
        metrics: A dictionary of study metrics to calculate FPKM for.

    Returns:
        A dictionary of study metrics with FPKM calculated.
    """
    for study in metrics:
        count_matrix = metrics[study].count_matrix.copy()
        length = metrics[study].eff_length

        if not isinstance(count_matrix, pd.DataFrame):
            raise TypeError("FPKM can only be performed on pandas.DataFrame objects")
        non_overlapping_cols = count_matrix.columns.difference(length.columns)
        if not non_overlapping_cols.empty:
            raise ValueError(
                f"Count matrix columns and effective length columns must match for FPKM calculation. "
                f"Non-overlapping columns: {non_overlapping_cols}"
            )

        mapped_reads = count_matrix.sum()
        metrics[study].normalization_matrix = (count_matrix * 1e9) / (mapped_reads * length)
    return metrics


def cpm_filter(
    *,
    context_name: str,
    metrics: NamedMetrics,
    filtering_options: _FilteringOptions,
    prep: RNAType,
) -> NamedMetrics:
    """Apply Counts Per Million (CPM) filtering to the count matrix for a given sample.

    Args:
        context_name: The name of the context being processed.
        metrics: A dictionary of study metrics to filter.
        filtering_options: Options for filtering the count matrix.
        prep: The RNA preparation type.

    Returns:
        A dictionary of filtered study metrics.
    """
    config = Config()
    n_exp = filtering_options.replicate_ratio
    n_top = filtering_options.high_replicate_ratio
    cut_off = filtering_options.cut_off

    sample: str
    metric: _StudyMetrics
    for sample, metric in metrics.items():
        counts: pd.DataFrame = metric.count_matrix
        entrez_ids: list[str] = metric.entrez_gene_ids
        library_size: pd.DataFrame = cast(pd.DataFrame, counts.sum(axis=1))

        # For library_sizes equal to 0, add 1 to prevent divide by 0
        # This will not impact the final counts per million calculation because the original counts are still 0
        #   thus, (0 / 1) * 1_000_000 = 0
        library_size[library_size == 0] = 1

        output_filepath = Path(config.result_dir, context_name, prep.value, f"CPM_Matrix_{prep.value}_{sample}.csv")
        output_filepath.parent.mkdir(parents=True, exist_ok=True)
        counts_per_million: pd.DataFrame = (counts / library_size) * 1_000_000
        counts_per_million.insert(0, "entrez_gene_ids", pd.Series(entrez_ids))
        logger.debug(f"Writing CPM matrix to {output_filepath}")
        counts_per_million.dropna(inplace=True)
        counts_per_million.to_csv(output_filepath, index=False)

        # TODO: Counts per million is adding ~61,500 columns (equal to the number of genes) for some reason.
        #  Most likely due to multiplying by 1_000_000, not exactly sure why

        min_samples = round(n_exp * len(counts.columns))  # noqa: F841
        top_samples = round(n_top * len(counts.columns))  # noqa: F841
        test_bools = pd.DataFrame({"entrez_gene_ids": entrez_ids})
        for i in range(len(counts_per_million.columns)):
            median_sum = float(np.median(np.sum(counts[:, i])))
            if np.isnan(cut_off):  # noqa: SIM108
                cutoff = float(10e6 / median_sum)
            else:
                cutoff = float(1e6 * cut_off) / median_sum
            test_bools = test_bools.merge(counts_per_million[counts_per_million.iloc[:, i] > cutoff])

    return metrics


def tpm_filter(
    *,
    metrics: NamedMetrics,
    filtering_options: _FilteringOptions,
    tpm_df: pd.DataFrame,
) -> NamedMetrics:
    """Apply quantile-based filtering to the TPM matrix for a given sample.

    :param metrics: A dictionary of study metrics to filter.
    :param filtering_options: Options for filtering the count matrix.
    :param tpm_df: A dataframe containing TPM values for each gene in each sample, used for TPM quantile filtering
    :returns: A dictionary of filtered study metrics.
    """
    # TODO: Write the TPM matrix to disk

    min_sample_expression = filtering_options.replicate_ratio
    high_confidence_sample_expression = filtering_options.high_replicate_ratio
    cut_off = filtering_options.cut_off

    # Map the TPM results directly onto the sample
    for sample, metric in metrics.items():
        metrics[sample].normalization_matrix = tpm_df[metric.sample_names]
        metrics[sample].normalization_matrix.index = tpm_df["ensembl_gene_id"]

    sample: str
    metric: _StudyMetrics
    for sample, metric in metrics.items():
        entrez_ids = metric.entrez_gene_ids
        gene_size = metric.gene_sizes
        tpm_matrix: pd.DataFrame = metric.normalization_matrix

        tpm_quantile = tpm_matrix[tpm_matrix > 0]

        # TODO: Is using quantile correct?
        # quantile_cutoff = np.nanquantile(a=tpm_quantile.values, q=1 - (cut_off / 100), axis=0)
        tpm_above_cutoff = tpm_matrix[tpm_matrix > filtering_options.cut_off]

        # Only keep `entrez_gene_ids` that pass `min_genes`
        min_samples = round(min_sample_expression * len(tpm_matrix.columns))
        min_func = k_over_a(min_samples, filtering_options.cut_off)
        min_genes: npt.NDArray[np.bool] = genefilter(tpm_above_cutoff, min_func)
        metric.entrez_gene_ids = [gene for gene, keep in zip(tpm_matrix.index, min_genes, strict=True) if keep]
        # metric.gene_sizes = np.asarray(gene for gene, keep in zip(gene_size, min_genes, strict=True) if keep)
        # metric.count_matrix = metric.count_matrix.iloc[min_genes, :]
        metrics[sample].normalization_matrix = metrics[sample].normalization_matrix.iloc[min_genes, :]

        # top_samples = round(high_confidence_sample_expression * len(tpm_matrix.columns))
        # top_func = k_over_a(top_samples, filtering_options.cut_off)
        # top_genes: npt.NDArray[np.bool] = genefilter(boolean_expression, top_func)
        # keep_top_genes = [gene for gene, keep in zip(entrez_ids, top_genes, strict=True) if keep]
        # metric.high_confidence_entrez_gene_ids = [
        #     gene for gene, keep in zip(entrez_ids, keep_top_genes, strict=True) if keep
        # ]

        # Calculate z-scores
        log_matrix = np.log1p(metrics[sample].normalization_matrix.values)
        metrics[sample].z_score_matrix = pd.DataFrame(
            sklearn.preprocessing.scale(log_matrix, axis=0),
            columns=metrics[sample].normalization_matrix.columns,
            index=metrics[sample].normalization_matrix.index,
        )

    return metrics


def zfpkm_filter(
    *,
    metrics: NamedMetrics,
    filtering_options: _FilteringOptions,
    calculate_fpkm: bool,
    force_zfpkm_plot: bool,
    output_fpkm_filepath: Path | None,
    output_png_dirpath: Path | None,
    force_negative_to_zero: bool = False,
) -> NamedMetrics:
    """Apply zFPKM filtering to the FPKM matrix for a given sample.

    Args:
        metrics: A dictionary of study metrics to filter.
        filtering_options: Options for filtering the count matrix.
        calculate_fpkm: Whether to calculate FPKM from counts.
        force_zfpkm_plot: Whether to force plotting of zFPKM results even if there are many samples.
        min_peak_height: Minimum peak height for zFPKM peak identification.
        min_peak_distance: Minimum peak distance for zFPKM peak identification.
        output_png_dirpath: Optional directory path to save zFPKM plots.
        force_negative_to_zero: Should negative values be forcibly set to 0?
            This could happen as a result of normalization producing negative near-zero values (e.g., -0.001)

    Returns:
        A dictionary of filtered study metrics.
    """
    min_sample_expression = filtering_options.replicate_ratio
    high_confidence_sample_expression = filtering_options.high_replicate_ratio
    cut_off = filtering_options.cut_off
    metrics = _calculate_fpkm(metrics) if calculate_fpkm else metrics

    if output_fpkm_filepath:
        combined_fpkm = pd.concat([df.normalization_matrix for df in metrics.values()], axis=1)
        combined_fpkm = combined_fpkm.sort_index().reindex(columns=sorted(combined_fpkm.columns))
        combined_fpkm.to_csv(output_fpkm_filepath, index=True)
        context = next(iter(metrics.values())).sample_names[0].partition("_")[0]
        logger.success(f"Wrote FPKM data for {context} to: {output_fpkm_filepath}")

    for metric in metrics.values():
        # if fpkm was not calculated, the normalization matrix will be empty; collect the count matrix instead
        matrix = metric.count_matrix if metric.normalization_matrix.empty else metric.normalization_matrix
        if not isinstance(matrix, pd.DataFrame):
            raise TypeError(f"Expected a pandas.DataFrame for zFPKM filtering, got: '{type(matrix)}'")

        # TODO: 2025-OCT-31: Re-evaluate whether to remove rows with all 0 counts
        # matrix = matrix[matrix.sum(axis=1) > 0]  # remove rows (genes) that have no counts across all samples

        matrix.replace(to_replace=np.nan, value=0.0, inplace=True)
        if force_negative_to_zero:
            matrix[matrix < 0] = 0.0

        zfpkm_df, zfpkm_results = zFPKM(matrix)

        if len(zfpkm_results) > 25 and not force_zfpkm_plot:
            logger.warning(
                "Not plotting zFPKM results because more than 25 plots would be created. "
                "If you would like to plot them anyway, set 'force_zfpkm_plot' to True"
            )
        elif output_png_dirpath is None:
            logger.critical("Output zFPKM PNG filepath is None, set a path to plot zFPKM graphs")
        else:
            sample_name = zfpkm_results[0].name.split("_")[0]  # go from 'control1hr_S1R1' to 'control1hr'
            zfpkm_plot(zfpkm_results, save_filepath=output_png_dirpath / f"{sample_name}_zfpkm_density.png")

        metric.z_score_matrix = zfpkm_df

        # determine which genes are expressed
        min_samples = round(min_sample_expression * len(zfpkm_df.columns))
        min_func = k_over_a(min_samples, cut_off)
        min_genes: npt.NDArray[np.bool] = genefilter(zfpkm_df, min_func)
        metric.entrez_gene_ids = np.asarray(
            [g_id for g_id, keep in zip(zfpkm_df.index, min_genes, strict=True) if keep], dtype=int
        )

        # determine which genes are confidently expressed
        top_samples = round(high_confidence_sample_expression * len(zfpkm_df.columns))
        top_func = k_over_a(top_samples, cut_off)
        top_genes: npt.NDArray[np.bool] = genefilter(zfpkm_df, top_func)
        metric.high_confidence_entrez_gene_ids = [
            gene for gene, keep in zip(zfpkm_df.index, top_genes, strict=True) if keep
        ]

    return metrics


def umi_filter(
    metrics: NamedMetrics,
    filtering_options: _FilteringOptions,
    target_sum: int = 10_000,
    perform_normalization: bool = False,
) -> NamedMetrics:
    """Perform UMI-based filtering.

    UMI filtering uses ScanPy's built-in `sc.pp.scale` (if `perform_normalization=True`)
    Otherwise, this function assumes that data has been pre-normalized+scaled beforehand and will evaluate expressed & highly expressed genes directly

    For each metric's matrix:
        - The rows are genomic identifiers (gene symbol, entrez gene id, ensembl gene id, etc.)
        - The columns are cell identifiers (e.g., barcodes)

    Calculating counts per cell should, therefore, be a column-wise sum (axis=0)


    :param metrics: The metrics to perform UMI filtering on
    :param filtering_options: Options for filtering the count matrix.
    :param target_sum: The target sum for UMI normalization.
    :param perform_normalization: Whether to perform normalization before filtering.

    :returns: The filtered metrics
    """
    min_sample_expression = filtering_options.replicate_ratio
    high_confidence_sample_expression = filtering_options.high_replicate_ratio
    cut_off = filtering_options.cut_off

    if min_sample_expression > 0.20:
        logger.warning(
            "Setting a minimum sample expression greater than ~20% for UMI-based filtering will likely result in very few/no genes being marked as active. "  # noqa: E501
            "Activity values ranging from 10-20% are recommended based on recent literature. "
            f"Got: {min_sample_expression:%} for option 'replicate_ratio'"
        )
    if high_confidence_sample_expression > 0.40:
        logger.warning(
            f"Setting high-confidence expression greater than ~40% for UMI-based filtering will likely result in very few to no genes being marked as highly active. "  # noqa: E501
            "Activity values ranging from 20-30% are recommended based on recent literature. "
            f"Got: {high_confidence_sample_expression:%} for option 'high_replicate_ratio'."
        )

    for metric in metrics.values():
        metric: _StudyMetrics
        if not isinstance(metric.count_matrix, sc.AnnData):
            raise TypeError(f"Expected a scanpy.AnnData for UMI filtering, got: '{type(metric.count_matrix)}'")
        adata: sc.AnnData = metric.count_matrix.copy()

        if perform_normalization:
            if adata.raw is not None:
                adata.X = adata.raw.X.copy()
            sc.pp.filter_cells(adata, min_genes=10)
            sc.pp.filter_genes(adata, min_cells=1)
            sc.pp.normalize_total(adata, target_sum=target_sum)
            sc.pp.log1p(adata)

        metric.z_score_matrix = adata

        adata_x = adata.X
        n_cells, n_genes = adata.shape

        min_samples: float = round(min_sample_expression * n_cells)
        min_func = k_over_a(min_samples, cut_off)
        min_genes_mask = np.zeros(n_genes, dtype=bool)
        for j in range(n_genes):
            col = adata_x.getcol(j).toarray().ravel() if sparse.issparse(adata_x) else adata_x[:, j]
            min_genes_mask[j] = min_func(col)
        # at this point we do not need/want NA entrez IDs
        metric.entrez_gene_ids = adata.var.loc[min_genes_mask, "entrez_gene_id"].dropna().tolist()

        top_samples = round(high_confidence_sample_expression * n_cells)
        top_func = k_over_a(top_samples, cut_off)
        top_genes_mask = np.zeros(n_genes, dtype=bool)
        for j in range(n_genes):
            col = adata_x.getcol(j).toarray().ravel() if sparse.issparse(adata_x) else adata_x[:, j]
            top_genes_mask[j] = top_func(col)
        metric.high_confidence_entrez_gene_ids = adata.var.loc[top_genes_mask, "entrez_gene_id"].dropna().tolist()

    return metrics


def filter_counts(
    *,
    context_name: str,
    metrics: NamedMetrics,
    technique: FilteringTechnique,
    filtering_options: _FilteringOptions,
    prep: RNAType,
    force_zfpkm_plot: bool,
    umi_target_sum: int = 10_000,
    tpm_df: pd.DataFrame | None,
    calculate_fpkm: bool = True,
    umi_perform_normalization: bool = False,
    output_fpkm_filepath: Path | None,
    output_zfpkm_plot_dirpath: Path | None = None,
    force_negative_to_zero: bool = False,
) -> NamedMetrics:
    """Filter the count matrix based on the specified technique.

    :param context_name: The name of the context being processed.
    :param metrics: A dictionary of study metrics to filter.
    :param technique: The filtering technique to use.
    :param filtering_options: Options for filtering the count matrix.
    :param prep: The RNA preparation type.
    :param force_zfpkm_plot: Whether to force plotting of zFPKM results even if there are many samples.
    :param umi_target_sum: The target sum for UMI normalization.
    :param tpm_df: A dataframe containing TPM values for each gene in each sample, used for TPM quantile filtering
    :param calculate_fpkm: If using paired-end data, should FPKM be calculated?
    :param umi_perform_normalization: Whether to perform normalization before UMI filtering.
    :param output_zfpkm_plot_dirpath: Optional filepath to save the zFPKM plot.
    :param output_fpkm_filepath: Optional filepath to save the FPKM matrix (only applies if `calculate_fpkm=True` and `technique=FilteringTechnique.ZFPKM`)
    :param force_negative_to_zero: Should negative values be forcibly set to 0?
            This could happen as a result of normalization producing negative near-zero values (e.g., -0.001)

    :returns: A dictionary of filtered study metrics.
    """  # noqa: E501
    match technique:
        case FilteringTechnique.CPM:
            return cpm_filter(
                context_name=context_name, metrics=metrics, filtering_options=filtering_options, prep=prep
            )
        case FilteringTechnique.TPM:
            if tpm_df is None:
                raise ValueError("TPM DataFrame must be provided for TPM quantile filtering")
            return tpm_filter(metrics=metrics, filtering_options=filtering_options, tpm_df=tpm_df)
        case FilteringTechnique.ZFPKM:
            return zfpkm_filter(
                metrics=metrics,
                filtering_options=filtering_options,
                calculate_fpkm=True,
                force_zfpkm_plot=force_zfpkm_plot,
                output_png_dirpath=output_zfpkm_plot_dirpath,
                output_fpkm_filepath=output_fpkm_filepath,
                force_negative_to_zero=force_negative_to_zero,
            )
        case FilteringTechnique.UMI:
            return umi_filter(
                metrics=metrics,
                filtering_options=filtering_options,
                target_sum=umi_target_sum,
                perform_normalization=umi_perform_normalization,
            )


def _process(
    context_name: str,
    rnaseq_matrix_filepath: Path,
    metadata_df: pd.DataFrame,
    gene_info_df: pd.DataFrame,
    fragment_df: pd.DataFrame,
    tpm_df: pd.DataFrame,
    prep: RNAType,
    taxon: int,
    replicate_ratio: float,
    batch_ratio: float,
    high_replicate_ratio: float,
    high_batch_ratio: float,
    technique: FilteringTechnique,
    cut_off: int | float,
    force_zfpkm_plot: bool,
    umi_target_sum: int,
    umi_perform_normalization: bool,
    output_boolean_activity_filepath: Path,
    output_zscore_normalization_filepath: Path,
    output_fpkm_filepath: Path | None,
    output_zfpkm_plot_dirpath: Path | None,
    force_negative_to_zero: bool,
):
    """Save the results of the RNA-Seq tests to a CSV file."""
    output_boolean_activity_filepath.parent.mkdir(parents=True, exist_ok=True)

    rnaseq_matrix: pd.DataFrame | sc.AnnData = read_file(rnaseq_matrix_filepath, h5ad_as_df=False)
    filtering_options = _FilteringOptions(
        replicate_ratio=replicate_ratio,
        batch_ratio=batch_ratio,
        cut_off=float(cut_off),
        high_replicate_ratio=high_replicate_ratio,
        high_batch_ratio=high_batch_ratio,
    )

    if fragment_df.empty and tpm_df.empty and prep != RNAType.SCRNA:
        raise ValueError("Effective length dataframe and TPM dataframe cannot both be empty for RNA-Seq processing")
    metrics, entrez_gene_ids = _build_matrix_results(
        rnaseq_matrix,
        gene_info=gene_info_df,
        metadata_df=metadata_df,
        normalize_df=fragment_df if not fragment_df.empty else tpm_df,
        taxon=taxon,
        normalize_df_is_fragment_length=not bool(fragment_df.empty),
    )
    metrics = filter_counts(
        context_name=context_name,
        metrics=metrics,
        technique=technique,
        filtering_options=filtering_options,
        prep=prep,
        tpm_df=tpm_df,
        force_zfpkm_plot=force_zfpkm_plot,
        umi_target_sum=umi_target_sum,
        umi_perform_normalization=umi_perform_normalization,
        output_zfpkm_plot_dirpath=output_zfpkm_plot_dirpath,
        output_fpkm_filepath=output_fpkm_filepath,
        force_negative_to_zero=force_negative_to_zero,
    )

    if isinstance(rnaseq_matrix, pd.DataFrame):
        merged_zscores = pd.concat(
            [m.z_score_matrix[m.z_score_matrix.index != "-"] for m in metrics.values()], axis="columns"
        )

        merged_zscores.index.name = (
            "entrez_gene_id"
            if merged_zscores.index.astype(str).str.isdigit().all()
            else "ensembl_gene_id"
            if merged_zscores.index.astype(str).str.startswith("ENS").all()
            else "gene_symbol"
        )

        merged_zscores = merged_zscores.reindex(columns=sorted(merged_zscores.columns))
        merged_zscores = merged_zscores.groupby(by=merged_zscores.index.name).mean()
        merged_zscores.to_csv(output_zscore_normalization_filepath, index=True)
    elif isinstance(rnaseq_matrix, sc.AnnData):
        merged_zscores = ad.concat([m.z_score_matrix for m in metrics.values()], axis="obs")
        merged_zscores.var.index.name = "entrez_gene_id"
        merged_zscores.obs = merged_zscores.obs.reindex(columns=sorted(merged_zscores.obs.columns))
        merged_zscores.write_h5ad(output_zscore_normalization_filepath.with_suffix(".h5ad"))
    logger.success(f"Wrote z-score normalization matrix to {output_zscore_normalization_filepath}")

    expressed_genes: list[str] = list(itertools.chain.from_iterable(m.entrez_gene_ids for m in metrics.values()))
    top_genes: list[str] = list(
        itertools.chain.from_iterable(m.high_confidence_entrez_gene_ids for m in metrics.values())
    )

    expression_frequency = pd.Series(expressed_genes).value_counts()
    expression_df = pd.DataFrame(
        {"entrez_gene_id": expression_frequency.index, "frequency": expression_frequency.values}
    )
    expression_df["prop"] = expression_df["frequency"] / len(metrics)
    expression_df = expression_df[expression_df["prop"] >= filtering_options.batch_ratio]

    top_frequency = pd.Series(top_genes).value_counts()
    top_df = pd.DataFrame({"entrez_gene_id": top_frequency.index, "frequency": top_frequency.values})
    top_df["prop"] = top_df["frequency"] / len(metrics)
    top_df = top_df[top_df["prop"] >= filtering_options.high_batch_ratio]

    entrez_id_series = pd.Series(entrez_gene_ids)
    boolean_matrix = pd.DataFrame(
        data={
            "entrez_gene_id": entrez_gene_ids,
            "expressed": entrez_id_series.isin(expression_df["entrez_gene_id"]).astype(int),
            "high": entrez_id_series.isin(top_df["entrez_gene_id"]).astype(int),
        }
    )

    expressed_count = len(boolean_matrix[boolean_matrix["expressed"] == 1])
    high_confidence_count = len(boolean_matrix[boolean_matrix["high"] == 1])

    # TODO: 2025-OCT-31: commented out dropping entrez ids, should this be kept?
    # boolean_matrix.dropna(subset="entrez_gene_id", inplace=True)
    boolean_matrix = boolean_matrix.groupby("entrez_gene_id", as_index=False).mean()
    boolean_matrix["expressed"] = boolean_matrix["expressed"].copy().round(0)
    boolean_matrix["high"] = boolean_matrix["high"].copy().astype(int)
    boolean_matrix.to_csv(output_boolean_activity_filepath, index=False)
    logger.info(
        f"{context_name} - Found {expressed_count} expressed genes, "
        f"{high_confidence_count} of which are confidently expressed"
    )
    logger.success(f"Wrote boolean matrix to {output_boolean_activity_filepath}")


def rnaseq_gen(  # noqa: C901
    context_name: str,
    input_rnaseq_filepath: Path,
    input_gene_info_filepath: Path,
    prep: RNAType,
    taxon_id: int,
    output_boolean_activity_filepath: Path,
    output_zscore_normalization_filepath: Path,
    input_metadata_filepath_or_df: Path | pd.DataFrame,
    replicate_ratio: float = 0.5,
    high_replicate_ratio: float = 0.8,
    batch_ratio: float = 0.5,
    high_batch_ratio: float = 0.75,
    technique: FilteringTechnique | str = FilteringTechnique.ZFPKM,
    umi_target_sum: int = 10_000,
    output_fpkm_filepath: Path | None = None,
    input_fragment_lengths: Path | None = None,
    input_tpm_filepath: Path | None = None,
    umi_perform_normalization: bool = False,
    cutoff: int | float | None = None,
    force_zfpkm_plot: bool = False,
    log_level: LogLevel = LogLevel.INFO,
    log_location: str | TextIO = sys.stderr,
    output_zfpkm_plot_dirpath: Path | None = None,
    force_negative_counts_to_zero: bool = False,
) -> None:
    """Generate a list of active and high-confidence genes from a gene count matrix.

    Replicates are compared for consensus within the study/batch number according to replicate ratios,
        then study/batch numbers are checked for consensus according to batch ratios.
    The zFPKM method is outlined here: https://pubmed.ncbi.nlm.nih.gov/24215113/

    :param context_name: The name of the context being processed
    :param input_rnaseq_filepath: The filepath to the gene count matrix
    :param input_gene_info_filepath: The filepath to the gene info file
    :param output_boolean_activity_filepath: The filepath to write the output gene count matrix
    :param output_zscore_normalization_filepath: The filepath to write the output z-score normalization matrix
    :param prep: The preparation method
    :param taxon_id: The NCBI Taxon ID
    :param input_metadata_filepath_or_df: The filepath or dataframe containing metadata information
    :param input_fragment_lengths: The filepath to the fragment lengths file, if applicable.
    :param input_tpm_filepath: The filepath to the TPM values file, if applicable.
    :param replicate_ratio: The percentage of replicates that a gene must
        appear in for a gene to be marked as "active" in a batch/study
    :param batch_ratio: The percentage of batches that a gene must appear in for a gene to be marked as 'active"
    :param high_replicate_ratio: The percentage of replicates that a gene must
        appear in for a gene to be marked "highly confident" in its expression in a batch/study
    :param high_batch_ratio: The percentage of batches that a gene must
        appear in for a gene to be marked "highly confident" in its expression
    :param technique: The filtering technique to use
    :param umi_target_sum: The target sum for UMI normalization
    :param output_fpkm_filepath: Optional filepath to write the FPKM matrix
        (required if using zFPKM filtering with paired-end data)
    :param umi_perform_normalization: Should UMI normalization be performed
    :param cutoff: The cutoff value to use for the provided filtering technique
        TPM: This is a TPM-based value; if `replicate_ratio` genes have a value less than this, the gene will be marked as inactive
        zFPKM: This is a zFPKM-based value; if `replicate_ratio` genes have a value less than this, the gene will be marked as inactive
    :param force_zfpkm_plot: If too many samples exist, should plotting be done anyway?
    :param log_level: The level of logging to output
    :param log_location: The location to write logs to
    :param output_zfpkm_plot_dirpath: Optional filepath to save zFPKM plots
    :param force_negative_counts_to_zero: Should negative values be forcibly set to 0?
        This could happen as a result of normalization producing negative near-zero values (e.g., -0.001)

    :return: None
    """  # noqa: E501
    set_up_logging(level=log_level, location=log_location)

    fragment_df = pd.DataFrame()
    technique = FilteringTechnique(technique) if isinstance(technique, str) else technique
    if technique == FilteringTechnique.TPM:
        if input_tpm_filepath is None:
            raise ValueError("Missing argument: `input_tpm_filepath` is required for TPM filtering")
        cutoff = cutoff or 1
    elif technique == FilteringTechnique.CPM:
        if isinstance(cutoff, int) and cutoff < 0:
            raise ValueError("Cutoff must be greater than or equal to 0")
        elif not cutoff:
            cutoff = np.nan
    elif technique == FilteringTechnique.ZFPKM:
        if not input_fragment_lengths:
            raise ValueError("Fragment lengths file is required for zFPKM filtering")
        if not output_fpkm_filepath:
            raise ValueError("Output FPKM filepath is required for zFPKM filtering")
        cutoff = cutoff or -3
        fragment_df = read_file(input_fragment_lengths)

    elif technique == FilteringTechnique.UMI:
        cutoff = cutoff or 1
    else:
        raise ValueError(f"Technique must be one of {','.join(FilteringTechnique)}. Got: {technique.value}")

    if not input_rnaseq_filepath.exists():
        raise FileNotFoundError(f"Input RNA-seq file not found! Searching for: '{input_rnaseq_filepath}'")

    if prep == RNAType.SCRNA and technique.value.lower() != FilteringTechnique.UMI.value.lower():
        logger.warning(
            "Single cell filtration does not normalize and assumes "
            "genes are counted with Unique Molecular Identifiers (UMIs). "
            f"Switching filtering technique from '{technique.value}' to '{FilteringTechnique.UMI.value}'."
        )
        technique = FilteringTechnique.UMI

    if isinstance(input_metadata_filepath_or_df, pd.DataFrame):
        metadata_df = input_metadata_filepath_or_df
    elif isinstance(input_metadata_filepath_or_df, Path):
        if input_metadata_filepath_or_df.suffix not in {".xlsx", ".xls"}:
            raise ValueError(
                f"Expected an excel file with extension of '.xlsx' or '.xls', "
                f"got '{input_metadata_filepath_or_df.suffix}'."
            )
        if not input_metadata_filepath_or_df.exists():
            raise FileNotFoundError(f"Input metadata file not found! Searching for: '{input_metadata_filepath_or_df}'")

        metadata_df = pd.read_excel(input_metadata_filepath_or_df)
    else:
        raise TypeError(
            f"Expected a pandas DataFrame or Path object as metadata, got '{type(input_metadata_filepath_or_df)}'"
        )

    _process(
        context_name=context_name,
        rnaseq_matrix_filepath=input_rnaseq_filepath,
        metadata_df=metadata_df,
        gene_info_df=read_file(input_gene_info_filepath),
        fragment_df=fragment_df,
        tpm_df=read_file(input_tpm_filepath) if input_tpm_filepath else pd.DataFrame(),
        prep=prep,
        taxon=taxon_id,
        replicate_ratio=replicate_ratio,
        batch_ratio=batch_ratio,
        high_replicate_ratio=high_replicate_ratio,
        high_batch_ratio=high_batch_ratio,
        technique=technique,
        cut_off=cutoff,
        force_zfpkm_plot=force_zfpkm_plot,
        umi_target_sum=umi_target_sum,
        umi_perform_normalization=umi_perform_normalization,
        output_fpkm_filepath=output_fpkm_filepath,
        output_boolean_activity_filepath=output_boolean_activity_filepath,
        output_zscore_normalization_filepath=output_zscore_normalization_filepath,
        output_zfpkm_plot_dirpath=output_zfpkm_plot_dirpath,
        force_negative_to_zero=force_negative_counts_to_zero,
    )


if __name__ == "__main__":
    import matplotlib.pyplot as plt

    data = pd.read_csv("/Users/joshl/Downloads/fpkm_example_data/CD8.genes.results.txt", sep="\t")
    data["gene_id"] = data["gene_id"].str.partition(".")[0]
    counts = (
        data[["gene_id", "expected_count"]]
        .copy()
        .set_index("gene_id")
        .sort_index()
        .rename(columns={"expected_count": "actual"})
    )
    eff_len = (
        data[["gene_id", "effective_length"]]
        .copy()
        .set_index("gene_id")
        .sort_index()
        .rename(columns={"effective_length": "actual"})
    )
    expected_fpkm = (
        data[["gene_id", "FPKM"]].copy().set_index("gene_id").sort_index().rename(columns={"FPKM": "expected"})
    )

    metrics = {
        "S1": _StudyMetrics(
            study="S1",
            num_samples=1,
            count_matrix=counts,
            eff_length=eff_len,
            sample_names=[""],
            layout=[LayoutMethod.paired_end],
            entrez_gene_ids=np.ndarray([0]),
            gene_sizes=np.ndarray([0]),
        )
    }
    calculated_fpkm = _calculate_fpkm(metrics)["S1"].normalization_matrix
    calculated_fpkm = calculated_fpkm.round(2)

    joined = calculated_fpkm.join(expected_fpkm, how="inner")
    joined["actual"] = joined["actual"].replace([np.nan, np.inf], 0)

    zfpkm_df, _ = zFPKM(joined, remove_na=True)
    zfpkm_df = zfpkm_df.replace(-np.inf, np.nan)

    fig, axes = cast(tuple[plt.Figure, list[plt.Axes]], plt.subplots(nrows=2, ncols=1))
    axes[0].hist(zfpkm_df["actual"].to_numpy())
    axes[0].set_title("Expected zFPKM")

    axes[1].hist(zfpkm_df["expected"].to_numpy())
    axes[1].set_title("Actual zFPKM")
    fig.tight_layout()
    fig.show()
