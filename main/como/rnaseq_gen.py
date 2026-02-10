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
import boolean
import numpy as np
import numpy.typing as npt
import pandas as pd
import scanpy as sc
import sklearn
import sklearn.neighbors
from anndata.compat import XDataArray
from anndata.experimental.backed import Dataset2D
from fast_bioservices.pipeline import ensembl_to_gene_id_and_symbol, gene_symbol_to_ensembl_and_gene_id
from loguru import logger
from scipy import sparse
from zfpkm import zFPKM, zfpkm_plot

from como.data_types import FilteringTechnique, LogLevel, RNAType
from como.migrations import gene_info_migrations
from como.project import Config
from como.utils import log_and_raise_error, read_file, set_up_logging


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
    fragment_lengths: npt.NDArray[np.floating] | None
    sample_names: list[str]
    layout: list[LayoutMethod]
    entrez_gene_ids: npt.NDArray[np.integer]
    gene_sizes: npt.NDArray[np.integer]
    __normalization_matrix: pd.DataFrame = field(default_factory=pd.DataFrame)
    __z_score_matrix: pd.DataFrame | sc.AnnData | None = field(default=None)
    __high_confidence_entrez_gene_ids: list[str] = field(default_factory=list)

    def __post_init__(self):
        for layout in self.layout:
            if layout not in LayoutMethod:
                log_and_raise_error(
                    f"Layout must be 'paired-end' or 'single-end'; got: {layout}",
                    error=ValueError,
                    level=LogLevel.ERROR,
                )

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
    if not isinstance(data, pd.DataFrame | np.ndarray):
        log_and_raise_error(
            f"Unsupported data type. Must be a Pandas DataFrame or a NumPy array, got '{type(data)}'",
            error=TypeError,
            level=LogLevel.CRITICAL,
        )

    return (
        data.apply(filter_func, axis=1).to_numpy()
        if isinstance(data, pd.DataFrame)
        else np.apply_along_axis(filter_func, axis=1, arr=data)
    )


async def _build_matrix_results(
    matrix: pd.DataFrame | sc.AnnData,
    *,
    gene_info: pd.DataFrame,
    metadata_df: pd.DataFrame,
    fragment_df: pd.DataFrame | None,
    taxon: int,
) -> tuple[NamedMetrics, list[int]]:
    """Read the counts matrix and returns the results.

    :param matrix: The gene counts matrix to process
    :param metadata_df: The configuration dataframe related to the current context
    :param taxon: The NCBI Taxon ID
    :returns: A dataclass `ReadMatrixResults`
    """
    if isinstance(matrix, sc.AnnData):
        if not isinstance(matrix.var, pd.DataFrame):
            raise TypeError("AnnData.var is expected to be a pandas.DataFrame")

        matrix.var = matrix.var.reset_index(drop=False, names=["gene_symbol"])
        conversion = await gene_symbol_to_ensembl_and_gene_id(symbols=matrix.var["gene_symbol"].tolist(), taxon=taxon)
    else:
        if "ensembl_gene_id" not in matrix.columns:
            log_and_raise_error(
                message="'ensembl_gene_id' column not found in the provided DataFrame.",
                error=ValueError,
                level=LogLevel.CRITICAL,
            )
        conversion: pd.DataFrame = await ensembl_to_gene_id_and_symbol(
            ids=matrix["ensembl_gene_id"].tolist(), taxon=taxon
        )
    # If the entrez gene id column is empty, it is indicative that the incorrect taxon id was provided
    if conversion["entrez_gene_id"].eq("-").all():
        logger.critical(
            f"Conversion of Ensembl Gene IDs to Entrez IDs and Gene Symbols was empty - "
            f"is '{taxon}' the correct taxon ID for this data?"
        )
    conversion["ensembl_gene_id"] = conversion["ensembl_gene_id"].str.split(",")
    conversion = conversion.explode("ensembl_gene_id")
    conversion = conversion[conversion["entrez_gene_id"] != "-"]
    conversion["entrez_gene_id"] = conversion["entrez_gene_id"]
    conversion = conversion.reset_index(drop=False)

    # conversion_merge_on should contain at least one of "ensembl_gene_id", "entrez_gene_id", or "gene_symbol"
    conversion_merge_on: list[str] = list(
        set(matrix.columns if isinstance(matrix, pd.DataFrame) else matrix.var.columns) & set(conversion.columns)
    )
    if not conversion_merge_on:
        log_and_raise_error(
            (
                "No columns to merge on, unable to find at least one of `ensembl_gene_id`, `entrez_gene_id`, or `gene_symbol`. "
                "Please check your input files."
            ),
            error=ValueError,
            level=LogLevel.ERROR,
        )

    if isinstance(matrix, pd.DataFrame):
        if "entrez_gene_id" in matrix.columns:
            matrix["entrez_gene_id"] = matrix["entrez_gene_id"].astype(int)
        matrix = matrix.merge(conversion, on=conversion_merge_on, how="left")
    elif isinstance(matrix, sc.AnnData):
        if "entrez_gene_id" in matrix.var.columns:
            matrix.var["entrez_gene_id"] = matrix.var["entrez_gene_id"].astype(int)
        matrix.var = matrix.var.merge(conversion, on=conversion_merge_on, how="left")

    gene_info = gene_info_migrations(gene_info)
    gene_info = gene_info[gene_info["entrez_gene_id"] != "-"]
    gene_info.loc[:, "entrez_gene_id"] = gene_info.loc[:, "entrez_gene_id"].astype(int)

    gene_info_merge_on: list[str] = list(
        set(matrix.columns if isinstance(matrix, pd.DataFrame) else matrix.var.columns) & set(gene_info.columns)
    )

    if "entrez_gene_id" in gene_info_merge_on:
        gene_info = gene_info[~gene_info["entrez_gene_id"].isna()]
        gene_info["entrez_gene_id"] = gene_info["entrez_gene_id"].astype(int)

        if isinstance(matrix, pd.DataFrame):
            matrix = matrix[~matrix["entrez_gene_id"].isna()]
            matrix["entrez_gene_id"] = matrix["entrez_gene_id"].astype(int)
        elif isinstance(matrix, sc.AnnData):
            if isinstance(matrix.var, XDataArray):
                raise TypeError("Expected matrix.var object to be 'pd.DataFrame', got 'anndata.compat.XDataArray'")
            matrix = matrix[:, ~matrix.var["entrez_gene_id"].isna()]
            matrix.var["entrez_gene_id"] = matrix.var["entrez_gene_id"].astype(int)

    if isinstance(matrix, pd.DataFrame):
        matrix = matrix.merge(gene_info, on=gene_info_merge_on, how="inner")
    elif isinstance(matrix, sc.AnnData):
        if not isinstance(matrix.var, pd.DataFrame):
            raise TypeError(f"Expected matrix.var object to be 'pd.DataFrame', got '{type(matrix.var)}'")
        matrix.var["original_index"] = matrix.var.index
        new_var = matrix.var.merge(gene_info, on=gene_info_merge_on, how="inner")
        new_matrix = matrix[:, new_var["original_index"]].copy()
        new_matrix.var = new_var
        new_matrix.var = new_matrix.var.drop(columns=["original_index"])
        new_matrix.var.reset_index(drop=True)
        matrix = new_matrix

        non_duplicates = ~matrix.var.duplicated(subset=matrix.var.columns, keep="first")
        matrix = matrix[:, non_duplicates].copy()

    metrics: NamedMetrics = {}
    for study in metadata_df["study"].unique():
        study_sample_names: list[str] = metadata_df[metadata_df["study"] == study]["sample_name"].tolist()
        layouts: list[str] = metadata_df[metadata_df["study"] == study]["layout"].tolist()

        if isinstance(matrix, pd.DataFrame):
            subset = matrix.set_index(keys=["entrez_gene_id"], drop=True)
            subset = subset[subset.columns.intersection(study_sample_names)]
            subset.index = subset.index.astype(int)
            entrez_gene_ids = subset.index.to_numpy(copy=False)
            gene_sizes = matrix["size"].to_numpy(dtype=int, copy=False)
        elif isinstance(matrix, sc.AnnData):
            # matrix.var = matrix.var.set_index(keys=["entrez_gene_id"], drop=True)
            subset = matrix[matrix.obs_names.intersection(study_sample_names)]
            entrez_gene_ids = subset.var["entrez_gene_id"].to_numpy(dtype=int)
            gene_sizes = subset.var["size"].to_numpy(dtype=int)
        else:
            log_and_raise_error(
                message=f"Matrix must be a pandas DataFrame or scanpy AnnData object, got: '{type(matrix)}'.",
                error=TypeError,
                level=LogLevel.CRITICAL,
            )

        frag_lengths = None
        if fragment_df is not None:
            frag_lengths = fragment_df["effective_length"].to_numpy(dtype=np.float64)
        metrics[study] = _StudyMetrics(
            count_matrix=subset,
            fragment_lengths=frag_lengths,
            sample_names=study_sample_names,
            layout=[LayoutMethod(layout) for layout in layouts],
            num_samples=len(study_sample_names),
            entrez_gene_ids=entrez_gene_ids,
            gene_sizes=gene_sizes,
            study=study,
        )

    return metrics, gene_info["entrez_gene_id"].astype(int).tolist()


def calculate_tpm(metrics: NamedMetrics) -> NamedMetrics:
    """Calculate the Transcripts Per Million (TPM) for each sample in the metrics dictionary.

    Args:
        metrics: A dictionary of study metrics to calculate TPM for.

    Returns:
        A dictionary of study metrics with TPM calculated.
    """
    for sample, metric in metrics.items():
        if isinstance(metric.count_matrix, sc.AnnData):
            adata = metric.count_matrix
            gene_sizes = pd.Series(metric.gene_sizes, index=adata.var_names)
            counts_df = pd.DataFrame(
                data=np.asarray(adata.X.toarray() if sparse.issparse(adata.X) else adata.X),
                index=adata.var_names,
                columns=adata.obs_names,
            )
        else:
            counts_df = metric.count_matrix
            gene_sizes = pd.Series(metric.gene_sizes)

        adjusted_counts = counts_df.add(1e-6)
        tpm_matrix = adjusted_counts.div(gene_sizes, axis=0)  # (count + 1) / gene_length
        tpm_matrix = tpm_matrix.div(tpm_matrix.sum(axis=0), axis=1)  # normalize by total
        tpm_matrix = tpm_matrix.mul(1e6)  # scale to per-million
        metrics[sample].normalization_matrix = tpm_matrix

    return metrics


def _calculate_fpkm(metrics: NamedMetrics, scale: float = 1e6) -> NamedMetrics:
    """Calculate the Fragments Per Kilobase of transcript per Million mapped reads (FPKM) for each i in the metrics dictionary.

    Args:
        metrics: A dictionary of study metrics to calculate FPKM for.
        scale: The scaling factor for normalization (default is 1e6).

    Returns:
        A dictionary of study metrics with FPKM calculated.
    """
    for study in metrics:
        matrix_values: dict[str, npt.NDArray[np.floating]] = {}
        count_matrix = metrics[study].count_matrix
        if not isinstance(count_matrix, pd.DataFrame):
            log_and_raise_error(
                message="FPKM cannot be performed on scanpy.AnnData objects!",
                error=TypeError,
                level=LogLevel.CRITICAL,
            )

        study_counts = count_matrix.to_numpy(dtype=int, copy=False)
        for i in range(metrics[study].num_samples):
            layout = metrics[study].layout[i]
            sample_name = metrics[study].sample_names[i]
            length = metrics[study].fragment_lengths if layout == LayoutMethod.paired_end else metrics[study].gene_sizes
            counts = study_counts[:, i]
            mapped_reads = counts.sum()
            matrix_values[sample_name] = ((counts * 1e9) / (length * mapped_reads)).astype(int)

        metrics[study].normalization_matrix = pd.DataFrame(matrix_values, index=metrics[study].entrez_gene_ids)
    return metrics


def calculate_z_score(metrics: NamedMetrics) -> NamedMetrics:
    """Calculate the z-score for each sample in the metrics dictionary.

    Args:
        metrics: A dictionary of study metrics to calculate z-scores for.

    Returns:
        A dictionary of study metrics with z-scores calculated.
    """
    for sample in metrics:
        log_matrix = np.log(metrics[sample].normalization_matrix)
        z_matrix = pd.DataFrame(
            data=sklearn.preprocessing.scale(log_matrix, axis=1), columns=metrics[sample].sample_names
        )
        metrics[sample].z_score_matrix = z_matrix
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
            if cut_off == "default":  # noqa: SIM108
                cutoff = float(10e6 / median_sum)
            else:
                cutoff = float(1e6 * cut_off) / median_sum
            test_bools = test_bools.merge(counts_per_million[counts_per_million.iloc[:, i] > cutoff])

    return metrics


def tpm_quantile_filter(*, metrics: NamedMetrics, filtering_options: _FilteringOptions) -> NamedMetrics:
    """Apply quantile-based filtering to the TPM matrix for a given sample.

    Args:
        metrics: A dictionary of study metrics to filter.
        filtering_options: Options for filtering the count matrix.

    Returns:
        A dictionary of filtered study metrics.
    """
    # TODO: Write the TPM matrix to disk

    n_exp = filtering_options.replicate_ratio
    n_top = filtering_options.high_replicate_ratio
    cut_off = filtering_options.cut_off
    metrics = calculate_tpm(metrics)

    sample: str
    metric: _StudyMetrics
    for sample, metric in metrics.items():
        entrez_ids = metric.entrez_gene_ids
        gene_size = metric.gene_sizes
        tpm_matrix: pd.DataFrame = metric.normalization_matrix

        min_samples = round(n_exp * len(tpm_matrix.columns))
        top_samples = round(n_top * len(tpm_matrix.columns))

        tpm_quantile = tpm_matrix[tpm_matrix > 0]
        quantile_cutoff = np.quantile(a=tpm_quantile.values, q=1 - (cut_off / 100), axis=0)
        boolean_expression = pd.DataFrame(
            data=tpm_matrix > quantile_cutoff, index=tpm_matrix.index, columns=tpm_matrix.columns
        ).astype(int)

        min_func = k_over_a(min_samples, 0.9)
        top_func = k_over_a(top_samples, 0.9)

        min_genes: npt.NDArray[np.bool] = genefilter(boolean_expression, min_func)
        top_genes: npt.NDArray[np.bool] = genefilter(boolean_expression, top_func)

        # Only keep `entrez_gene_ids` that pass `min_genes`
        metric.entrez_gene_ids = [gene for gene, keep in zip(entrez_ids, min_genes, strict=True) if keep]
        metric.gene_sizes = np.asarray(gene for gene, keep in zip(gene_size, min_genes, strict=True) if keep)
        metric.count_matrix = cast(pd.DataFrame, metric.count_matrix.iloc[min_genes, :])
        metric.normalization_matrix = cast(pd.DataFrame, metrics[sample].normalization_matrix.iloc[min_genes, :])

        keep_top_genes = [gene for gene, keep in zip(entrez_ids, top_genes, strict=True) if keep]
        metric.high_confidence_entrez_gene_ids = [
            gene for gene, keep in zip(entrez_ids, keep_top_genes, strict=True) if keep
        ]

    metrics = calculate_z_score(metrics)

    return metrics


def zfpkm_filter(
    *,
    metrics: NamedMetrics,
    filtering_options: _FilteringOptions,
    calculate_fpkm: bool,
    force_zfpkm_plot: bool,
    min_peak_height: float,
    min_peak_distance: int,
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

    for metric in metrics.values():
        metric: _StudyMetrics
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
            f"Got: {min_sample_expression} for option 'replicate_ratio'"
        )
    if high_confidence_sample_expression > 0.40:
        logger.warning(
            f"Setting high-confidence expression greater than ~40% for UMI-based filtering will likely result in very few to no genes being marked as highly active. "  # noqa: E501
            "Activity values ranging from 20-30% are recommended based on recent literature. "
            f"Got: {high_confidence_sample_expression} for option 'high_replicate_ratio'."
        )

    for metric in metrics.values():
        metric: _StudyMetrics
        if not isinstance(metric.count_matrix, sc.AnnData):
            raise TypeError(f"Expected a scanpy.AnnData for UMI filtering, got: '{type(metric.count_matrix)}'")
        adata: sc.AnnData = metric.count_matrix

        if perform_normalization:
            if adata.raw is not None:
                adata.X = adata.raw.X.copy()
            sc.pp.filter_cells(adata, min_genes=20)
            sc.pp.filter_genes(adata, min_cells=1)
            sc.pp.normalize_total(adata, target_sum=target_sum)
            sc.pp.log1p(adata)
            # sc.pp.scale(adata, max_value=15)  # abs(values)>10 standard deviations away will be set to +/-10

        metric.z_score_matrix = adata

        adata_x = adata.X
        n_cells, n_genes = adata.shape

        min_samples: float = round(min_sample_expression * n_cells)
        min_func = k_over_a(min_samples, cut_off)
        min_genes_mask = np.zeros(n_genes, dtype=bool)
        for j in range(n_genes):
            col = adata_x.getcol(j).toarray().ravel() if sparse.issparse(adata_x) else adata_x[:, j]
            min_genes_mask[j] = min_func(col)
        metric.entrez_gene_ids = (
            adata.var.loc[min_genes_mask, "entrez_gene_id"].dropna().tolist()
        )  # at this point we do not need/want NA entrez IDs

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
    zfpkm_min_peak_height: float,
    zfpkm_min_peak_distance: int,
    umi_target_sum: int = 10_000,
    umi_perform_normalization: bool = False,
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
    :param zfpkm_min_peak_height: Minimum peak height for zFPKM peak identification.
    :param zfpkm_min_peak_distance: Minimum peak distance for zFPKM peak identification.
    :param umi_target_sum: The target sum for UMI normalization.
    :param umi_perform_normalization: Whether to perform normalization before UMI filtering.
    :param output_zfpkm_plot_dirpath: Optional filepath to save the zFPKM plot.
    :param force_negative_to_zero: Should negative values be forcibly set to 0?
            This could happen as a result of normalization producing negative near-zero values (e.g., -0.001)

    :returns: A dictionary of filtered study metrics.
    """
    match technique:
        case FilteringTechnique.CPM:
            return cpm_filter(
                context_name=context_name, metrics=metrics, filtering_options=filtering_options, prep=prep
            )
        case FilteringTechnique.TPM:
            return tpm_quantile_filter(metrics=metrics, filtering_options=filtering_options)
        case FilteringTechnique.ZFPKM:
            return zfpkm_filter(
                metrics=metrics,
                filtering_options=filtering_options,
                calculate_fpkm=True,
                force_zfpkm_plot=force_zfpkm_plot,
                min_peak_height=zfpkm_min_peak_height,
                min_peak_distance=zfpkm_min_peak_distance,
                output_png_dirpath=output_zfpkm_plot_dirpath,
                force_negative_to_zero=force_negative_to_zero,
            )
        case FilteringTechnique.UMI:
            return umi_filter(
                metrics=metrics,
                filtering_options=filtering_options,
                target_sum=umi_target_sum,
                perform_normalization=umi_perform_normalization,
            )
        case _:
            log_and_raise_error(
                f"Technique must be one of {FilteringTechnique}, got '{technique.value}'",
                error=ValueError,
                level=LogLevel.ERROR,
            )


async def _process(
    context_name: str,
    rnaseq_matrix_filepath: Path,
    metadata_df: pd.DataFrame,
    gene_info_df: pd.DataFrame,
    fragment_df: pd.DataFrame | None,
    prep: RNAType,
    taxon: int,
    replicate_ratio: float,
    batch_ratio: float,
    high_replicate_ratio: float,
    high_batch_ratio: float,
    technique: FilteringTechnique,
    cut_off: int | float,
    force_zfpkm_plot: bool,
    zfpkm_min_peak_height: float,
    zfpkm_min_peak_distance: int,
    umi_target_sum: int,
    umi_perform_normalization: bool,
    output_boolean_activity_filepath: Path,
    output_zscore_normalization_filepath: Path,
    output_zfpkm_plot_dirpath: Path | None,
    force_negative_to_zero: bool,
):
    """Save the results of the RNA-Seq tests to a CSV file."""
    output_boolean_activity_filepath.parent.mkdir(parents=True, exist_ok=True)

    rnaseq_matrix: pd.DataFrame | sc.AnnData = _read_file(rnaseq_matrix_filepath, h5ad_as_df=False)
    filtering_options = _FilteringOptions(
        replicate_ratio=replicate_ratio,
        batch_ratio=batch_ratio,
        cut_off=float(cut_off),
        high_replicate_ratio=high_replicate_ratio,
        high_batch_ratio=high_batch_ratio,
    )

    metrics, entrez_gene_ids = await _build_matrix_results(
        rnaseq_matrix,
        gene_info=gene_info_df,
        metadata_df=metadata_df,
        fragment_df=fragment_df,
        taxon=taxon,
    )
    metrics = filter_counts(
        context_name=context_name,
        metrics=metrics,
        technique=technique,
        filtering_options=filtering_options,
        prep=prep,
        force_zfpkm_plot=force_zfpkm_plot,
        zfpkm_min_peak_height=zfpkm_min_peak_height,
        zfpkm_min_peak_distance=zfpkm_min_peak_distance,
        umi_target_sum=umi_target_sum,
        umi_perform_normalization=umi_perform_normalization,
        output_zfpkm_plot_dirpath=output_zfpkm_plot_dirpath,
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
        merged_zscores = merged_zscores.groupby("entrez_gene_id").mean()
        merged_zscores.to_csv(output_zscore_normalization_filepath, index=True)
    elif isinstance(rnaseq_matrix, sc.AnnData):
        merged_zscores = ad.concat([m.z_score_matrix for m in metrics.values()], axis="obs")
        merged_zscores.var.index.name = "entrez_gene_id"
        merged_zscores.obs = merged_zscores.obs.reindex(columns=sorted(merged_zscores.obs.columns))
        merged_zscores.write_h5ad(output_zscore_normalization_filepath.with_suffix(".h5ad"))
    expressed_genes: list[str] = list(itertools.chain.from_iterable(m.entrez_gene_ids for m in metrics.values()))
    top_genes: list[str] = list(
        itertools.chain.from_iterable(m.high_confidence_entrez_gene_ids for m in metrics.values())
    )

    logger.success(f"Wrote z-score normalization matrix to {output_zscore_normalization_filepath}")

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
    boolean_matrix["expressed"] = boolean_matrix["expressed"].copy().astype(int)
    boolean_matrix["high"] = boolean_matrix["high"].copy().astype(int)
    boolean_matrix.to_csv(output_boolean_activity_filepath, index=False)
    logger.info(
        f"{context_name} - Found {expressed_count} expressed genes, {high_confidence_count} of which are confidently expressed"
    )
    logger.success(f"Wrote boolean matrix to {output_boolean_activity_filepath}")


async def rnaseq_gen(  # noqa: C901
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
    zfpkm_min_peak_height: float = 0.02,
    zfpkm_min_peak_distance: int = 1,
    umi_target_sum: int = 10_000,
    input_fragment_lengths: Path | None = None,
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
    :param replicate_ratio: The percentage of replicates that a gene must
        appear in for a gene to be marked as "active" in a batch/study
    :param batch_ratio: The percentage of batches that a gene must appear in for a gene to be marked as 'active"
    :param high_replicate_ratio: The percentage of replicates that a gene must
        appear in for a gene to be marked "highly confident" in its expression in a batch/study
    :param high_batch_ratio: The percentage of batches that a gene must
        appear in for a gene to be marked "highly confident" in its expression
    :param technique: The filtering technique to use
    :param zfpkm_min_peak_height: The height of the zFPKM peak
    :param zfpkm_min_peak_distance: The distance of the zFPKM peak
    :param umi_target_sum: The target sum for UMI normalization
    :param umi_perform_normalization: Should UMI normalization be performed?
    :param cutoff: The cutoff value to use for the provided filtering technique
    :param force_zfpkm_plot: If too many samples exist, should plotting be done anyway?
    :param log_level: The level of logging to output
    :param log_location: The location to write logs to
    :param output_zfpkm_plot_dirpath: Optional filepath to save zFPKM plots
    :param force_negative_counts_to_zero: Should negative values be forcibly set to 0?
        This could happen as a result of normalization producing negative near-zero values (e.g., -0.001)

    :return: None
    """
    set_up_logging(level=log_level, location=log_location)

    technique = FilteringTechnique(technique) if isinstance(technique, str) else technique
    match technique:
        case FilteringTechnique.TPM:
            cutoff: int | float = cutoff or 25
            if cutoff < 1 or cutoff > 100:
                log_and_raise_error(
                    "Quantile must be between 1 - 100",
                    error=ValueError,
                    level=LogLevel.ERROR,
                )

        case FilteringTechnique.CPM:
            if cutoff and cutoff < 0:
                log_and_raise_error(
                    "Cutoff must be greater than or equal to 0",
                    error=ValueError,
                    level=LogLevel.ERROR,
                )
            elif cutoff:
                cutoff = "default"

        case FilteringTechnique.ZFPKM:
            cutoff: int | float = cutoff or -3
        case FilteringTechnique.UMI:
            cutoff: int = cutoff or 1
        case _:
            log_and_raise_error(
                f"Technique must be one of {','.join(FilteringTechnique)}. Got: {technique.value}",
                error=ValueError,
                level=LogLevel.ERROR,
            )

    if not input_rnaseq_filepath.exists():
        log_and_raise_error(
            f"Input RNA-seq file not found! Searching for: '{input_rnaseq_filepath}'",
            error=FileNotFoundError,
            level=LogLevel.ERROR,
        )

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
            log_and_raise_error(
                f"Expected an excel file with extension of '.xlsx' or '.xls', got '{input_metadata_filepath_or_df.suffix}'.",
                error=ValueError,
                level=LogLevel.ERROR,
            )
        if not input_metadata_filepath_or_df.exists():
            log_and_raise_error(
                f"Input metadata file not found! Searching for: '{input_metadata_filepath_or_df}'",
                error=FileNotFoundError,
                level=LogLevel.ERROR,
            )

        metadata_df = pd.read_excel(input_metadata_filepath_or_df)
    else:
        log_and_raise_error(
            f"Expected a pandas DataFrame or Path object as metadata, got '{type(input_metadata_filepath_or_df)}'",
            error=TypeError,
            level=LogLevel.ERROR,
        )

    logger.debug(f"Starting '{context_name}'")
    await _process(
        context_name=context_name,
        rnaseq_matrix_filepath=input_rnaseq_filepath,
        metadata_df=metadata_df,
        gene_info_df=_read_file(input_gene_info_filepath),
        fragment_df=_read_file(input_fragment_lengths),
        prep=prep,
        taxon=taxon_id,
        replicate_ratio=replicate_ratio,
        batch_ratio=batch_ratio,
        high_replicate_ratio=high_replicate_ratio,
        high_batch_ratio=high_batch_ratio,
        technique=technique,
        cut_off=cutoff,
        force_zfpkm_plot=force_zfpkm_plot,
        zfpkm_min_peak_height=zfpkm_min_peak_height,
        zfpkm_min_peak_distance=zfpkm_min_peak_distance,
        umi_target_sum=umi_target_sum,
        umi_perform_normalization=umi_perform_normalization,
        output_boolean_activity_filepath=output_boolean_activity_filepath,
        output_zscore_normalization_filepath=output_zscore_normalization_filepath,
        output_zfpkm_plot_dirpath=output_zfpkm_plot_dirpath,
        force_negative_to_zero=force_negative_counts_to_zero,
    )
