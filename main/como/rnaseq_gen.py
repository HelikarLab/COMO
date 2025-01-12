from __future__ import annotations

import multiprocessing
import sys
import time
from collections import namedtuple
from concurrent.futures import Future, ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from enum import Enum
from io import TextIOWrapper
from pathlib import Path
from typing import Callable, NamedTuple

import numpy as np
import numpy.typing as npt
import pandas as pd
import plotly.graph_objs as go
import scanpy as sc
import sklearn
import sklearn.neighbors
from fast_bioservices.pipeline import ensembl_to_gene_id_and_symbol, gene_symbol_to_ensembl_and_gene_id
from loguru import logger
from pandas import DataFrame
from plotly.subplots import make_subplots
from scipy.signal import find_peaks
from sklearn.neighbors import KernelDensity

from como.data_types import FilteringTechnique, LogLevel, PeakIdentificationParameters, RNAType
from como.migrations import gene_info_migrations
from como.project import Config
from como.utils import _log_and_raise_error, _num_columns, _read_file, _set_up_logging


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


@dataclass
class _StudyMetrics:
    study: str
    num_samples: int
    count_matrix: pd.DataFrame
    fragment_lengths: npt.NDArray[np.float32]
    sample_names: list[str]
    layout: list[LayoutMethod]
    entrez_gene_ids: list[str]
    gene_sizes: npt.NDArray[np.float32]
    __normalization_matrix: pd.DataFrame = field(default_factory=pd.DataFrame)
    __z_score_matrix: pd.DataFrame = field(default_factory=pd.DataFrame)
    __high_confidence_entrez_gene_ids: list[str] = field(default=list)

    def __post_init__(self):
        for layout in self.layout:
            if layout not in LayoutMethod:
                _log_and_raise_error(
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
    def z_score_matrix(self) -> pd.DataFrame:
        return self.__z_score_matrix

    @z_score_matrix.setter
    def z_score_matrix(self, value: pd.DataFrame) -> None:
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
    max_fpkm: float


class _ReadMatrixResults(NamedTuple):
    metrics: dict[str, _StudyMetrics]
    entrez_gene_ids: list[str]


Density = namedtuple("Density", ["x", "y"])
NamedMetrics = dict[str, _StudyMetrics]


def k_over_a(k: int, a: float) -> Callable[[npt.NDArray], bool]:
    """Return a function that filters rows of an array based on the sum of elements being greater than or equal to A at least k times.

    This code is based on the `kOverA` function found in R's `genefilter` package: https://www.rdocumentation.org/packages/genefilter/versions/1.54.2/topics/kOverA

    :param k: The minimum number of times the sum of elements must be greater than or equal to A.
    :param a: The value to compare the sum of elements to.
    :return: A function that accepts a NumPy array to perform the actual filtering
    """  # noqa: E501

    def filter_func(row: npt.NDArray) -> bool:
        return np.sum(row >= a) >= k

    return filter_func


def genefilter(data: pd.DataFrame | npt.NDArray, filter_func: Callable[[npt.NDArray], bool]) -> npt.NDArray:
    """Apply a filter function to the rows of the data and return the filtered array.

    This code is based on the `genefilter` function found in R's `genefilter` package: https://www.rdocumentation.org/packages/genefilter/versions/1.54.2/topics/genefilter

    :param data: The data to filter
    :param filter_func: THe function to filter the data by
    :return: A NumPy array of the filtered data.
    """
    if not isinstance(data, (pd.DataFrame, npt.NDArray)):
        _log_and_raise_error(
            f"Unsupported data type. Must be a Pandas DataFrame or a NumPy array, got '{type(data)}'",
            error=TypeError,
            level=LogLevel.CRITICAL,
        )

    return (
        data.apply(filter_func, axis=1).values
        if isinstance(data, pd.DataFrame)
        else np.apply_along_axis(filter_func, axis=1, arr=data)
    )


async def _build_matrix_results(
    *,
    matrix: pd.DataFrame,
    gene_info: pd.DataFrame,
    metadata_df: pd.DataFrame,
    taxon: int,
) -> _ReadMatrixResults:
    """Read the counts matrix and returns the results.

    :param matrix: The gene counts matrix to process
    :param metadata_df: The configuration dataframe related to the current context
    :param taxon: The NCBI Taxon ID
    :return: A dataclass `ReadMatrixResults`
    """
    matrix.dropna(inplace=True)
    conversion = await ensembl_to_gene_id_and_symbol(ids=matrix["ensembl_gene_id"].tolist(), taxon=taxon)
    conversion["ensembl_gene_id"] = conversion["ensembl_gene_id"].str.split(",")
    conversion = conversion.explode("ensembl_gene_id")
    conversion.reset_index(inplace=True, drop=True)
    matrix = matrix.merge(conversion, on=["ensembl_gene_id", "entrez_gene_id", "gene_symbol"], how="left")

    # Only include Entrez and Ensembl Gene IDs that are present in `gene_info`
    matrix["entrez_gene_id"] = matrix["entrez_gene_id"].str.split("//")
    matrix = matrix.explode("entrez_gene_id")
    matrix = matrix.replace(to_replace="-", value=pd.NA).dropna()
    matrix["entrez_gene_id"] = matrix["entrez_gene_id"].astype(int)

    gene_info = gene_info_migrations(gene_info)
    gene_info = gene_info.replace(to_replace="-", value=pd.NA).dropna()
    gene_info["entrez_gene_id"] = gene_info["entrez_gene_id"].astype(int)

    counts_matrix = matrix.merge(
        gene_info[["entrez_gene_id", "ensembl_gene_id"]],
        on=["entrez_gene_id", "ensembl_gene_id"],
        how="inner",
    )

    gene_info = gene_info.merge(
        counts_matrix[["entrez_gene_id", "ensembl_gene_id"]],
        on=["entrez_gene_id", "ensembl_gene_id"],
        how="inner",
    )

    entrez_gene_ids: list[str] = gene_info["entrez_gene_id"].tolist()
    metrics: NamedMetrics = {}
    for study in metadata_df["study"].unique().tolist():
        study_sample_names = metadata_df[metadata_df["study"] == study]["sample_name"].tolist()
        layouts = metadata_df[metadata_df["study"] == study]["layout"].tolist()
        metrics[study] = _StudyMetrics(
            count_matrix=counts_matrix[counts_matrix.columns.intersection(study_sample_names)],
            fragment_lengths=metadata_df[metadata_df["study"] == study]["fragment_length"].values,
            sample_names=study_sample_names,
            layout=[LayoutMethod(layout) for layout in layouts],
            num_samples=len(study_sample_names),
            entrez_gene_ids=entrez_gene_ids,
            gene_sizes=np.array(gene_info["size"].values).astype(np.float32),
            study=study,
        )
        metrics[study].fragment_lengths[np.isnan(metrics[study].fragment_lengths)] = 0
        metrics[study].count_matrix.index = pd.Index(entrez_gene_ids, name="entrez_gene_id")

    return _ReadMatrixResults(metrics=metrics, entrez_gene_ids=gene_info["entrez_gene_id"].tolist())


def calculate_tpm(metrics: NamedMetrics) -> NamedMetrics:
    """Calculate the Transcripts Per Million (TPM) for each sample in the metrics dictionary."""
    for sample in metrics:
        count_matrix = metrics[sample].count_matrix

        gene_sizes = metrics[sample].gene_sizes

        tpm_matrix = pd.DataFrame(data=None, index=count_matrix.index, columns=count_matrix.columns)
        for i in range(len(count_matrix.columns)):
            values: pd.Series = count_matrix.iloc[:, i] + 1  # Add 1 to prevent division by 0
            rate = np.log(values.tolist()) - np.log(gene_sizes)
            denominator = np.log(np.sum(np.exp(rate)))
            tpm_value = np.exp(rate - denominator + np.log(1e6))
            tpm_matrix.iloc[:, i] = tpm_value
        metrics[sample].normalization_matrix = tpm_matrix

    return metrics


def _calculate_fpkm(metrics: NamedMetrics, scale: int = 1e6) -> NamedMetrics:
    """Calculate the Fragments Per Kilobase of transcript per Million mapped reads (FPKM) for each sample in the metrics dictionary."""  # noqa: E501
    for study in metrics:
        matrix_values = []
        for sample in range(metrics[study].num_samples):
            layout = metrics[study].layout[sample]
            count_matrix: npt.NDArray = metrics[study].count_matrix.iloc[:, sample].values.astype(np.float32)
            gene_lengths = (
                metrics[study].fragment_lengths[sample].astype(np.float32)
                if layout == LayoutMethod.paired_end
                else metrics[study].gene_sizes.astype(np.float32)
            )
            gene_lengths_kb = gene_lengths / 1000.0

            match layout:
                case LayoutMethod.paired_end:  # FPKM
                    total_fragments = count_matrix.sum(axis=0)
                    counts_per_million = total_fragments / scale
                    fragments_per_kilobase = count_matrix / gene_lengths_kb[:, np.newaxis]
                    fragments_per_kilobase_million = fragments_per_kilobase / counts_per_million
                    matrix_values.append(fragments_per_kilobase_million)
                case LayoutMethod.single_end:  # RPKM
                    reads_per_kilobase = count_matrix / gene_lengths_kb[: np.newaxis]
                    total_reads = count_matrix.sum(axis=0)
                    counts_per_million = total_reads / scale
                    reads_per_kilobase_million = reads_per_kilobase / counts_per_million
                    matrix_values.append(reads_per_kilobase_million)
                case _:
                    _log_and_raise_error(
                        (
                            f"Invalid normalization method specified ''. "
                            f"Must be one of '{LayoutMethod.paired_end.value}' or '{LayoutMethod.single_end.value}'."
                        ),
                        error=ValueError,
                        level=LogLevel.ERROR,
                    )

        fpkm_matrix = pd.DataFrame(matrix_values).T  # Transpose is needed because values were appended as rows
        fpkm_matrix = fpkm_matrix[~pd.isna(fpkm_matrix)]
        metrics[study].normalization_matrix = fpkm_matrix
        metrics[study].normalization_matrix.columns = metrics[study].count_matrix.columns

    return metrics


# def _zfpkm_calculation(row: pd.Series, kernel: KernelDensity, peak_parameters: tuple[float, float]) -> _ZFPKMResult:
def _zfpkm_calculation(
    column: pd.Series,
    peak_parameters: PeakIdentificationParameters,
    bandwidth: int = 0.5,
) -> _ZFPKMResult:
    """Log2 Transformations.

    Stabilize the variance in the data to make the distribution more symmetric; this is helpful for Gaussian fitting

    Kernel Density Estimation (kde)
        - Non-parametric method to estimate the probability density function (PDF) of a random variable
        - Estimates the distribution of log2-transformed FPKM values
        - Bandwidth parameter controls the smoothness of the density estimate
        - KDE Explanation
            - A way to smooth a histogram to get a better idea of the underlying distribution of the data
            - Given a set of data points, we want to understand how they are distributed.
                Histograms can be useful, but are sensitive to bin size and number
            - The KDE places a "kernel" - a small symmetric function (i.e., Gaussian curve) - at each data point
            - The "kernel" acts as a weight, giving more weight to points closer to the center of the kernel,
                and less weight to points farther away
            - Kernel functions are summed along each point on the x-axis
            - A smooth curve is created that represents the estimated density of the data

    Peak Finding
        - Identifies that are above a certain height and separated by a minimum distance
        - Represent potential local maxima in the distribution

    Peak Selection
        - The peak with the highest x-value (from log2-FPKM) is chosen as the mean (mu)
            of the "inactive" gene distribution
        - The peak representing unexpressed or inactive genes should be at a lower FPKM
            value compared to the peak representing expressed genes

    Standard Deviation Estimation
         - The mean of log2-FPKM values are greater than the calculated mu
         - Standard deviation is estimated based on the assumption that the right tail of the distribution
            This represents expressed genes) can be approximated by a half-normal distribution

     zFPKM Transformation
        - Centers disbribution around 0 and scales it by the standard deviation.
            This makes it easier to compare gene expression across different samples
        - Represents the number of standard deviations away from the mean of the "inactive" gene distribution
        - Higher zFPKM values indicate higher expression levels relative to the "inactive" peak
        - A zFPKM value of 0 represents the mean of the "inactive" distribution
        - Research shows that a zFPKM value of -3 or greater can be used as
            a threshold for calling a gene as "active" and/or "expressed"
            : https://doi.org/10.1186/1471-2164-14-778
    """
    values = column.values
    refit: KernelDensity = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(values.reshape(-1, 1))  # type: ignore

    x_range = np.linspace(values.min(), values.max(), 2000)
    density = np.exp(refit.score_samples(x_range.reshape(-1, 1)))
    peaks, _ = find_peaks(density, height=peak_parameters.height, distance=peak_parameters.distance)
    peak_positions = x_range[peaks]

    mu = 0
    max_fpkm = 0
    stddev = 1

    if len(peaks) != 0:
        mu = peak_positions.max()
        max_fpkm = density[peaks[np.argmax(peak_positions)]]
        u = values[values > mu].mean()
        stddev = (u - mu) * np.sqrt(np.pi / 2)
    zfpkm = pd.Series((values - mu) / stddev, dtype=np.float32, name=column.name)

    return _ZFPKMResult(zfpkm=zfpkm, density=Density(x_range, density), mu=mu, std_dev=stddev, max_fpkm=max_fpkm)


def zfpkm_transform(
    fpkm_df: pd.DataFrame,
    bandwidth: int = 0.5,
    peak_parameters: tuple[float, float] = (0.02, 1.0),
    update_every_percent: float = 0.1,
) -> tuple[dict[str, _ZFPKMResult], DataFrame]:
    """Perform zFPKM calculation/transformation."""
    if update_every_percent > 1:
        logger.warning(
            f"update_every_percent should be a decimal value between 0 and 1; got: {update_every_percent} - "
            f"will convert to percentage"
        )
        update_every_percent /= 100

    total_samples = _num_columns(fpkm_df)
    update_per_step: int = int(np.ceil(total_samples * update_every_percent))

    # Get at least 1 core and at most cpu_count() - 2
    cores = max(min(multiprocessing.cpu_count() - 2, total_samples), 1)
    logger.debug(
        f"zFPKM transforming {len(fpkm_df.columns)} sample(s) "
        f"containing {len(fpkm_df):,} genes(s) using {cores} core(s)"
    )
    logger.debug(f"Will update every {update_per_step:,} steps (~{update_every_percent:.1%} of {total_samples:,})")

    chunk_time = time.time()
    start_time = time.time()
    log_padding = len(str(f"{total_samples:,}"))
    zfpkm_series: list[pd.Series | None] = [None] * total_samples
    results: dict[str, _ZFPKMResult] = {}

    with ProcessPoolExecutor(max_workers=cores) as pool:
        futures: list[Future[_ZFPKMResult]] = [
            pool.submit(
                _zfpkm_calculation,
                column=fpkm_df[column],
                peak_parameters=peak_parameters,
                bandwidth=bandwidth,
            )
            for column in fpkm_df
        ]
        for i, future in enumerate(as_completed(futures)):
            result = future.result()
            key = str(result.zfpkm.name)
            results[key] = result
            zfpkm_series[i] = result.zfpkm

            if i != 0 and ((i + 1) % update_per_step == 0 or (i + 1) == total_samples):
                current_time = time.time()
                chunk = current_time - chunk_time
                total_time = current_time - start_time
                chunk_num = f"{i + 1:,}"
                logger.debug(
                    f"Processed {chunk_num:>{log_padding}} of {total_samples:,} - "
                    f"chunk took {chunk:.1f} seconds - "
                    f"running for {total_time:.1f} seconds"
                )
                chunk_time = current_time

    zfpkm_df = pd.DataFrame({series.name: series for series in zfpkm_series}, index=fpkm_df.index)
    return results, zfpkm_df


def zfpkm_plot(results, *, plot_xfloor: int = -4, subplot_titles: bool = True):
    """Plot the log2(FPKM) density and fitted Gaussian for each sample.

    :param results: A dictionary of intermediate results from zfpkm_transform.
    :param: subplot_titles: Whether to display facet titles (sample names).
    :param plot_xfloor: Lower limit for the x-axis.
    :param subplot_titles: Whether to display facet titles (sample names).
    """
    to_concat: list[pd.DataFrame] = [None] * len(results)  # type: ignore  # ignoring because None is not of type pd.DataFrame
    for name, result in results.items():
        stddev = result.std_dev
        x = np.array(result.density.x)
        y = np.array(result.density.y)

        fitted = np.exp(-0.5 * ((x - result.mu) / stddev) ** 2) / (stddev * np.sqrt(2 * np.pi))
        max_fpkm = y.max()
        max_fitted = fitted.max()
        scale_fitted = fitted * (max_fpkm / max_fitted)

        to_concat.append(
            pd.DataFrame(
                {
                    "sample_name": [name] * len(x),
                    "log2fpkm": x,
                    "fpkm_density": y,
                    "fitted_density_scaled": scale_fitted,
                }
            )
        )
    mega_df = pd.concat(to_concat, ignore_index=True)
    mega_df.columns = pd.Series(data=["sample_name", "log2fpkm", "fpkm_density", "fitted_density_scaled"])

    mega_df = mega_df.melt(id_vars=["log2fpkm", "sample_name"], var_name="source", value_name="density")
    subplot_titles = list(results.keys()) if subplot_titles else None
    fig = make_subplots(
        rows=len(results),
        cols=1,
        subplot_titles=subplot_titles,
        vertical_spacing=min(0.05, (1 / (len(results) - 1))),
    )

    for i, (name, group) in enumerate(mega_df.groupby("sample_name"), start=1):
        fig.add_trace(
            trace=go.Scatter(x=group["log2fpkm"], y=group["density"], mode="lines", name=name, legendgroup=name),
            row=i,
            col=1,
        )
        fig.update_xaxes(title_text="log2(FPKM)", range=[plot_xfloor, max(group["log2fpkm"].tolist())], row=i, col=1)
        fig.update_yaxes(title_text="density [scaled]", row=i, col=1)
        fig.update_layout(legend_tracegroupgap=0)

    fig.update_layout(height=600 * len(results), width=1000, title_text="zFPKM Plots", showlegend=True)
    fig.write_image("zfpkm_plot.png")


def calculate_z_score(metrics: NamedMetrics) -> NamedMetrics:
    """Calculate the z-score for each sample in the metrics dictionary."""
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
    prep: RNAPrepMethod,
) -> NamedMetrics:
    """Apply Counts Per Million (CPM) filtering to the count matrix for a given sample."""
    config = Config()
    n_exp = filtering_options.replicate_ratio
    n_top = filtering_options.high_replicate_ratio
    cut_off = filtering_options.cut_off

    sample: str
    metric: _StudyMetrics
    for sample, metric in metrics.items():
        counts: pd.DataFrame = metric.count_matrix
        entrez_ids: list[str] = metric.entrez_gene_ids
        library_size: pd.DataFrame = counts.sum(axis=1)

        # For library_sizes equal to 0, add 1 to prevent divide by 0
        # This will not impact the final counts per million calculation because the original counts are still 0
        #   thus, (0 / 1) * 1_000_000 = 0
        library_size[library_size == 0] = 1

        output_filepath = config.result_dir / context_name / prep.value / f"CPM_Matrix_{prep.value}_{sample}.csv"
        output_filepath.parent.mkdir(parents=True, exist_ok=True)
        counts_per_million: pd.DataFrame = (counts / library_size) * 1_000_000
        counts_per_million.insert(0, "entrez_gene_ids", pd.Series(entrez_ids))
        logger.debug(f"Writing CPM matrix to {output_filepath}")
        counts_per_million.to_csv(output_filepath, index=False)

        # TODO: Counts per million is adding ~61,500 columns (equal to the number of genes) for some reason.
        #  Most likely due to multiplying by 1_000_000, not exactly sure why

        min_samples = round(n_exp * len(counts.columns))  # noqa: F841
        top_samples = round(n_top * len(counts.columns))  # noqa: F841
        test_bools = pd.DataFrame({"entrez_gene_ids": entrez_ids})
        for i in range(len(counts_per_million.columns)):
            cutoff = (
                10e6 / (np.median(np.sum(counts[:, i])))
                if cut_off == "default"
                else (1e6 * cut_off) / np.median(np.sum(counts[:, i]))
            )
            test_bools = test_bools.merge(counts_per_million[counts_per_million.iloc[:, i] > cutoff])

    return metrics


def tpm_quantile_filter(*, metrics: NamedMetrics, filtering_options: _FilteringOptions) -> NamedMetrics:
    """Apply quantile-based filtering to the TPM matrix for a given sample."""
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
        quantile_cutoff = np.quantile(
            a=tpm_quantile.values, q=1 - (cut_off / 100), axis=0
        )  # Compute quantile across columns
        boolean_expression = pd.DataFrame(
            data=tpm_matrix > quantile_cutoff, index=tpm_matrix.index, columns=tpm_matrix.columns
        ).astype(int)

        min_func = k_over_a(min_samples, 0.9)
        top_func = k_over_a(top_samples, 0.9)

        min_genes: npt.NDArray[bool] = genefilter(boolean_expression, min_func)
        top_genes: npt.NDArray[bool] = genefilter(boolean_expression, top_func)

        # Only keep `entrez_gene_ids` that pass `min_genes`
        metric.entrez_gene_ids = [gene for gene, keep in zip(entrez_ids, min_genes) if keep]
        metric.gene_sizes = np.array(gene for gene, keep in zip(gene_size, min_genes) if keep)
        metric.count_matrix = metric.count_matrix.iloc[min_genes, :]
        metric.normalization_matrix = metrics[sample].normalization_matrix.iloc[min_genes, :]

        keep_top_genes = [gene for gene, keep in zip(entrez_ids, top_genes) if keep]
        metric.high_confidence_entrez_gene_ids = [gene for gene, keep in zip(entrez_ids, keep_top_genes) if keep]

    metrics = calculate_z_score(metrics)

    return metrics


def zfpkm_filter(*, metrics: NamedMetrics, filtering_options: _FilteringOptions, calculate_fpkm: bool) -> NamedMetrics:
    """Apply zFPKM filtering to the FPKM matrix for a given sample."""
    min_sample_expression = filtering_options.replicate_ratio
    high_confidence_sample_expression = filtering_options.high_replicate_ratio
    cut_off = filtering_options.cut_off
    metrics = _calculate_fpkm(metrics) if calculate_fpkm else metrics

    metric: _StudyMetrics
    for metric in metrics.values():
        # if fpkm was not calculated, the normalization matrix will be empty; collect the count matrix instead
        matrix = metric.count_matrix if metric.normalization_matrix.empty else metric.normalization_matrix
        matrix = matrix[matrix.sum(axis=1) > 0]  # remove rows (genes) that have no counts

        minimums = matrix == 0
        results, zfpkm_df = zfpkm_transform(matrix)
        zfpkm_df[minimums] = -4
        zfpkm_plot(results)

        # determine which genes are expressed
        min_samples = round(min_sample_expression * len(zfpkm_df.columns))
        min_func = k_over_a(min_samples, cut_off)
        min_genes: npt.NDArray[bool] = genefilter(zfpkm_df, min_func)
        metric.entrez_gene_ids = [gene for gene, keep in zip(metric.entrez_gene_ids, min_genes) if keep]

        # determine which genes are confidently expressed
        top_samples = round(high_confidence_sample_expression * len(zfpkm_df.columns))
        top_func = k_over_a(top_samples, cut_off)
        top_genes: npt.NDArray[bool] = genefilter(zfpkm_df, top_func)
        metric.high_confidence_entrez_gene_ids = [gene for gene, keep in zip(metric.entrez_gene_ids, top_genes) if keep]

    return metrics


def filter_counts(
    *,
    context_name: str,
    metrics: NamedMetrics,
    technique: FilteringTechnique,
    filtering_options: _FilteringOptions,
    prep: RNAPrepMethod,
) -> NamedMetrics:
    """Filter the count matrix based on the specified technique."""
    match technique:
        case FilteringTechnique.cpm:
            return cpm_filter(
                context_name=context_name, metrics=metrics, filtering_options=filtering_options, prep=prep
            )
        case FilteringTechnique.tpm:
            return tpm_quantile_filter(metrics=metrics, filtering_options=filtering_options)
        case FilteringTechnique.zfpkm:
            return zfpkm_filter(metrics=metrics, filtering_options=filtering_options, calculate_fpkm=True)
        case FilteringTechnique.umi:
            return zfpkm_filter(metrics=metrics, filtering_options=filtering_options, calculate_fpkm=False)
        case _:
            raise ValueError(f"Technique must be one of {FilteringTechnique}")


async def _save_rnaseq_tests(
    context_name: str,
    rnaseq_matrix: pd.DataFrame,
    metadata_df: pd.DataFrame,
    gene_info_df: pd.DataFrame,
    prep: RNAPrepMethod,
    taxon: int,
    replicate_ratio: float,
    batch_ratio: float,
    high_replicate_ratio: float,
    high_batch_ratio: float,
    technique: FilteringTechnique,
    cut_off: int | float,
    output_boolean_activity_filepath: Path,
    output_zscore_normalization_filepath: Path,
):
    """Save the results of the RNA-Seq tests to a CSV file."""
    filtering_options = _FilteringOptions(
        replicate_ratio=replicate_ratio,
        batch_ratio=batch_ratio,
        cut_off=float(cut_off),
        high_replicate_ratio=high_replicate_ratio,
        high_batch_ratio=high_batch_ratio,
    )

    read_counts_results: _ReadMatrixResults = await _build_matrix_results(
        matrix=rnaseq_matrix,
        gene_info=gene_info_df,
        metadata_df=metadata_df,
        taxon=taxon,
    )
    metrics = read_counts_results.metrics
    entrez_gene_ids = read_counts_results.entrez_gene_ids

    metrics: NamedMetrics = filter_counts(
        context_name=context_name,
        metrics=metrics,
        technique=technique,
        filtering_options=filtering_options,
        prep=prep,
    )

    merged_zscore_df = pd.DataFrame()
    expressed_genes: list[str] = []
    top_genes: list[str] = []
    for metric in metrics.values():
        expressed_genes.extend(metric.entrez_gene_ids)
        top_genes.extend(metric.high_confidence_entrez_gene_ids)
        if metric.normalization_matrix is not None:
            merged_zscore_df = (
                metric.normalization_matrix
                if merged_zscore_df.empty
                else pd.concat(
                    [merged_zscore_df, metric.normalization_matrix],
                    axis=1,
                )
            )
    merged_zscore_df.index = pd.Series(entrez_gene_ids, name="entrez_gene_id")
    merged_zscore_df.to_csv(output_zscore_normalization_filepath, index=True)
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

    boolean_matrix = pd.DataFrame(data={"entrez_gene_id": entrez_gene_ids, "expressed": 0, "high": 0})
    for gene in entrez_gene_ids:
        if gene in expression_df["entrez_gene_id"]:
            boolean_matrix.loc[gene, "expressed"] = 1
        if gene in top_df["entrez_gene_id"]:
            boolean_matrix.loc[gene, "high"] = 1

    expressed_count = len(boolean_matrix[boolean_matrix["expressed"] == 1])
    high_confidence_count = len(boolean_matrix[boolean_matrix["high"] == 1])

    boolean_matrix.to_csv(output_boolean_activity_filepath, index=False)
    logger.info(
        f"{context_name} - Found {expressed_count} expressed and {high_confidence_count} confidently expressed genes"
    )
    logger.success(f"Wrote boolean matrix to {output_boolean_activity_filepath}")


async def _create_metadata_df(path: Path) -> pd.DataFrame:
    if path.suffix not in {".xls", ".xlsx"}:
        raise ValueError(
            f"Expected an excel file with extension of '.xlsx' or '.xls', got '{path.suffix}'. "
            f"Attempted to process: {path}"
        )
    return pd.read_excel(path)


async def rnaseq_gen(
    context_name: str,
    input_rnaseq_filepath: Path,
    input_gene_info_filepath: Path,
    prep: RNAPrepMethod,
    taxon_id: int,
    output_boolean_activity_filepath: Path,
    output_zscore_normalization_filepath: Path,
    input_metadata_filepath_or_df: Path | pd.DataFrame,
    replicate_ratio: float = 0.5,
    high_replicate_ratio: float = 1.0,
    batch_ratio: float = 0.5,
    high_batch_ratio: float = 1.0,
    technique: FilteringTechnique | str = FilteringTechnique.tpm,
    cutoff: int | float | None = None,
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
    :param replicate_ratio: The percentage of replicates that a gene must
        appear in for a gene to be marked as "active" in a batch/study
    :param batch_ratio: The percentage of batches that a gene must appear in for a gene to be marked as 'active"
    :param high_replicate_ratio: The percentage of replicates that a gene must
        appear in for a gene to be marked "highly confident" in its expression in a batch/study
    :param high_batch_ratio: The percentage of batches that a gene must
        appear in for a gene to be marked "highly confident" in its expression
    :param technique: The filtering technique to use
    :param cutoff: The cutoff value to use for the provided filtering technique
    :return: None
    """
    technique = (
        FilteringTechnique.from_string(str(technique.lower())) if isinstance(technique, (str, int)) else technique
    )

    match technique:
        case FilteringTechnique.tpm:
            cutoff = cutoff or 25
            if cutoff < 1 or cutoff > 100:
                raise ValueError("Quantile must be between 1 - 100")

        case FilteringTechnique.cpm:
            if cutoff and cutoff < 0:
                raise ValueError("Cutoff must be greater than 0")
            elif cutoff:
                cutoff = "default"

        case FilteringTechnique.zfpkm:
            cutoff = cutoff or -3
        case FilteringTechnique.umi:
            pass
        case _:
            raise ValueError(f"Technique must be one of {FilteringTechnique}")

    if not input_rnaseq_filepath.exists():
        raise FileNotFoundError(f"Input RNA-seq file not found! Searching for: '{input_rnaseq_filepath}'")

    if prep == RNAPrepMethod.SCRNA:
        technique = FilteringTechnique.umi
        logger.warning(
            "Single cell filtration does not normalize and assumes "
            "gene counts are counted with Unique Molecular Identifiers (UMIs). "
            "Setting filtering technique to UMI now."
        )

    logger.debug(f"Starting '{context_name}'")
    output_boolean_activity_filepath.parent.mkdir(parents=True, exist_ok=True)
    metadata_df = (
        input_metadata_filepath_or_df
        if isinstance(input_metadata_filepath_or_df, pd.DataFrame)
        else await _create_metadata_df(input_metadata_filepath_or_df)
    )
    await _save_rnaseq_tests(
        context_name=context_name,
        rnaseq_matrix=await _read_counts(input_rnaseq_filepath),
        metadata_df=metadata_df,
        gene_info_df=pd.read_csv(input_gene_info_filepath),
        prep=prep,
        taxon=taxon_id,
        replicate_ratio=replicate_ratio,
        batch_ratio=batch_ratio,
        high_replicate_ratio=high_replicate_ratio,
        high_batch_ratio=high_batch_ratio,
        technique=technique,
        cut_off=cutoff,
        output_boolean_activity_filepath=output_boolean_activity_filepath,
        output_zscore_normalization_filepath=output_zscore_normalization_filepath,
    )
