from __future__ import annotations

import multiprocessing
import sys
import time
from collections import namedtuple
from collections.abc import Callable
from concurrent.futures import Future, ProcessPoolExecutor, as_completed
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import NamedTuple, TextIO

import matplotlib.pyplot as plt
import numpy as np
import numpy.typing as npt
import pandas as pd
import seaborn as sns
import sklearn
import sklearn.neighbors
import sklearn.preprocessing
from fast_bioservices.pipeline import ensembl_to_gene_id_and_symbol, gene_symbol_to_ensembl_and_gene_id
from loguru import logger
from pandas import DataFrame
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


@dataclass(slots=True)
class _StudyMetrics:
    study: str
    num_samples: int
    count_matrix: pd.DataFrame
    fragment_lengths: npt.NDArray[float]
    sample_names: list[str]
    layout: list[LayoutMethod]
    entrez_gene_ids: npt.NDArray[int]
    gene_sizes: npt.NDArray[int]
    __normalization_matrix: pd.DataFrame = field(default_factory=pd.DataFrame)
    __z_score_matrix: pd.DataFrame = field(default_factory=pd.DataFrame)
    __high_confidence_entrez_gene_ids: list[str] = field(default_factory=list)

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
    """

    def filter_func(row: npt.NDArray) -> bool:
        return np.sum(row >= a) >= k

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
    if not isinstance(data, pd.DataFrame | npt.NDArray):
        _log_and_raise_error(
            f"Unsupported data type. Must be a Pandas DataFrame or a NumPy array, got '{type(data)}'",
            error=TypeError,
            level=LogLevel.CRITICAL,
        )

    return data.apply(filter_func, axis=1).values if isinstance(data, pd.DataFrame) else np.apply_along_axis(filter_func, axis=1, arr=data)


async def _build_matrix_results(
    *,
    matrix: pd.DataFrame,
    gene_info: pd.DataFrame,
    metadata_df: pd.DataFrame,
    taxon: int,
) -> _ReadMatrixResults:
    """Read the counts matrix and returns the results.

    Arg:
        matrix: The gene counts matrix to process
        metadata_df: The configuration dataframe related to the current context
        taxon: The NCBI Taxon ID

    Returns:
        A dataclass `ReadMatrixResults`
    """
    matrix.dropna(subset="ensembl_gene_id", inplace=True)
    conversion = await ensembl_to_gene_id_and_symbol(ids=matrix["ensembl_gene_id"].tolist(), taxon=taxon)
    conversion["ensembl_gene_id"] = conversion["ensembl_gene_id"].str.split(",")
    conversion = conversion.explode("ensembl_gene_id")
    conversion.reset_index(inplace=True, drop=True)

    # merge_on should contain at least one of "ensembl_gene_id", "entrez_gene_id", or "gene_symbol"
    merge_on: list[str] = list(set(matrix.columns).intersection(conversion.columns))
    if not merge_on:
        _log_and_raise_error(
            (
                "No columns to merge on, unable to find at least one of `ensembl_gene_id`, `entrez_gene_id`, or `gene_symbol`. "
                "Please check your input files."
            ),
            error=ValueError,
            level=LogLevel.ERROR,
        )
    matrix = matrix.merge(conversion, on=merge_on, how="left")
    if "entrez_gene_id" in matrix.columns:
        matrix["entrez_gene_id"] = matrix["entrez_gene_id"].astype(int)
    matrix = matrix.merge(conversion, on=merge_on, how="left")
    # drop rows that have `0` in `entrez_gene_id` column
    matrix = matrix[matrix["entrez_gene_id"] != 0].reset_index(drop=True, inplace=False)

    gene_info = gene_info[gene_info["entrez_gene_id"] != 0].reset_index(drop=True, inplace=False)
    gene_info = gene_info_migrations(gene_info)
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

    entrez_gene_ids: npt.NDArray[int] = gene_info["entrez_gene_id"].to_numpy()
    metrics: NamedMetrics = {}
    for study in metadata_df["study"].unique():
        study_sample_names = metadata_df[metadata_df["study"] == study]["sample_name"].tolist()
        layouts = metadata_df[metadata_df["study"] == study]["layout"].tolist()
        metrics[study] = _StudyMetrics(
            count_matrix=counts_matrix[counts_matrix.columns.intersection(study_sample_names)],
            fragment_lengths=metadata_df[metadata_df["study"] == study]["fragment_length"].values.astype(float),
            sample_names=study_sample_names,
            layout=[LayoutMethod(layout) for layout in layouts],
            num_samples=len(study_sample_names),
            entrez_gene_ids=entrez_gene_ids,
            gene_sizes=gene_info["size"].values.astype(int),
            study=study,
        )
        metrics[study].fragment_lengths[np.isnan(metrics[study].fragment_lengths)] = 0
        metrics[study].count_matrix.index = pd.Index(entrez_gene_ids, name="entrez_gene_id")

    return _ReadMatrixResults(metrics=metrics, entrez_gene_ids=gene_info["entrez_gene_id"].tolist())


def calculate_tpm(metrics: NamedMetrics) -> NamedMetrics:
    """Calculate the Transcripts Per Million (TPM) for each sample in the metrics dictionary.

    Args:
        metrics: A dictionary of study metrics to calculate TPM for.

    Returns:
        A dictionary of study metrics with TPM calculated.
    """
    for sample in metrics:
        count_matrix = metrics[sample].count_matrix

        gene_sizes = pd.Series(metrics[sample].gene_sizes, index=count_matrix.index)
        adjusted_counts = count_matrix.add(1e-6)

        tpm_matrix = adjusted_counts.divide(gene_sizes, axis=0)  # (count + 1) / gene_length
        tpm_matrix = tpm_matrix.div(tpm_matrix.sum(axis=0), axis=1)  # normalize by total
        tpm_matrix = tpm_matrix.mul(1e6)  # scale to per-million
        metrics[sample].normalization_matrix = tpm_matrix

    return metrics


def _calculate_fpkm(metrics: NamedMetrics, scale: int = 1e6) -> NamedMetrics:
    """Calculate the Fragments Per Kilobase of transcript per Million mapped reads (FPKM) for each sample in the metrics dictionary.

    Args:
        metrics: A dictionary of study metrics to calculate FPKM for.
        scale: The scaling factor for normalization (default is 1e6).

    Returns:
        A dictionary of study metrics with FPKM calculated.
    """
    for study in metrics:
        matrix_values: list[npt.NDArray[float]] = []
        for sample in range(metrics[study].num_samples):
            layout = metrics[study].layout[sample]
            count_matrix: npt.NDArray[int] = metrics[study].count_matrix.iloc[:, sample].values

            gene_lengths = metrics[study].fragment_lengths[sample] if layout == LayoutMethod.paired_end else metrics[study].gene_sizes.astype(int)
            gene_lengths_kb = gene_lengths / 1000.0

            match layout:
                case LayoutMethod.paired_end:  # FPKM
                    total_fragments: npt.NDArray[int] = count_matrix.sum(axis=0)
                    if total_fragments == 0:
                        fragments_per_kilobase_million: float = np.nan
                    else:
                        counts_per_million: npt.NDArray[float] = total_fragments / scale
                        fragments_per_kilobase: npt.NDArray[float] = count_matrix / counts_per_million
                        fragments_per_kilobase_million: npt.NDArray[float] = fragments_per_kilobase / gene_lengths_kb
                    matrix_values.append(fragments_per_kilobase_million)
                case LayoutMethod.single_end:  # RPKM
                    reads_per_kilobase = count_matrix / gene_lengths_kb
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

        # Transpose is needed because values were appended as rows
        fpkm_matrix = pd.DataFrame(matrix_values).T
        fpkm_matrix.index = metrics[study].count_matrix.index
        fpkm_matrix.columns = metrics[study].sample_names

        fpkm_matrix = fpkm_matrix[~pd.isna(fpkm_matrix)]
        metrics[study].normalization_matrix = fpkm_matrix
        metrics[study].normalization_matrix.columns = metrics[study].count_matrix.columns

    return metrics


def _zfpkm_calculation(
    column: pd.Series,
    peak_parameters: PeakIdentificationParameters,
    bandwidth: float = 1.0,
    epsilon: float = 1e-10,
) -> _ZFPKMResult:
    """Log2 Transformations.

    Stabilize the variance in the data to make the distribution more symmetric; this is helpful for Gaussian fitting

    Kernel Density Estimation (kde)
        - SciKit Learn: https://scikit-learn.org/stable/modules/density.html
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

    Args:
        column: A Pandas Series containing FPKM values for a single sample.
        peak_parameters: Parameters for peak identification in zFPKM calculation.
        bandwidth: The bandwidth for kernel density estimation in zFPKM calculation.
        epsilon: A small value to add to FPKM values to prevent log2-divide-by-0 errors

    Returns:
        A named tuple containing:
            - zfpkm: A Pandas Series of zFPKM values for the input sample.
            - density: A named tuple containing the x and y values of the KDE.
            - mu: The mean of the "inactive" gene distribution.
            - std_dev: The estimated standard deviation of the "inactive" gene distribution.
            - max_fpkm: The maximum FPKM value at the identified peak.
    """
    log2values: npt.NDArray[float] = np.log2(column.values + epsilon)

    # 1D KDE *always* has one feature and many samples. scikit-learn expects data in the format of (n_samples, n_features)
    # Thus, we use `.reshape(-1, 1)` because we know there is a single feature
    # Even though the data is FPKM of many genes for a single sample, it is still one feature over many samples
    # `-1` indicates the unknown dimension (the number of samples)
    # `1` indicates the known dimension (the number of genes [also known as features])
    # https://scikit-learn.org/stable/auto_examples/neighbors/plot_kde_1d.html#sphx-glr-auto-examples-neighbors-plot-kde-1d-py
    kde: KernelDensity = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(log2values.reshape(-1, 1))

    x_range: npt.NDArray[float] = np.linspace(log2values.min(), log2values.max(), 2000).reshape(-1, 1)
    density: npt.NDArray[float] = np.exp(kde.score_samples(x_range))
    peaks, _ = find_peaks(density, height=peak_parameters.height, distance=peak_parameters.distance)
    peak_positions = x_range[peaks]

    mu = 0
    max_fpkm = 0
    stddev = 1
    if len(peaks) != 0:
        mu = peak_positions.max()
        max_fpkm = density[peaks[np.argmax(peak_positions)]]
        u = log2values[log2values > mu].mean()
        stddev = (u - mu) * np.sqrt(np.pi / 2)
    zfpkm = pd.Series((log2values - mu) / stddev, dtype=float, name=column.name)

    return _ZFPKMResult(zfpkm=zfpkm, density=Density(x_range, density), mu=mu, std_dev=stddev, max_fpkm=max_fpkm)


def zfpkm_transform(
    fpkm_df: pd.DataFrame,
    peak_parameters: PeakIdentificationParameters,
    bandwidth: float,
    update_every_percent: float = 0.1,
) -> tuple[dict[str, _ZFPKMResult], DataFrame]:
    """Perform zFPKM calculation/transformation.

    Args:
        fpkm_df: A DataFrame containing FPKM values with genes as rows and samples as columns.
        peak_parameters: Parameters for peak identification in zFPKM calculation.
        bandwidth: The bandwidth for kernel density estimation in zFPKM calculation.
        update_every_percent: Frequency of progress updates as a decimal between 0 and 1 (e.g., 0.1 for every 10%).

    Returns:
        A tuple containing:
            - A dictionary of intermediate results for each sample.
            - A DataFrame of zFPKM values with the same shape as the input fpkm_df.
    """
    if update_every_percent > 1:
        logger.warning(f"update_every_percent should be a decimal value between 0 and 1; got: {update_every_percent} - will convert to percentage")
        update_every_percent /= 100

    total_samples = _num_columns(fpkm_df)
    update_per_step: int = int(np.ceil(total_samples * update_every_percent))

    # Get at least 1 core and at most cpu_count() - 2
    cores = max(min(multiprocessing.cpu_count() - 2, total_samples), 1)
    logger.debug(f"zFPKM transforming {len(fpkm_df.columns)} sample(s) containing {len(fpkm_df):,} genes(s) using {cores} core(s)")
    logger.debug(f"Will update every {update_per_step:,} steps (~{update_every_percent:.1%} of {total_samples:,})")

    chunk_time = time.time()
    start_time = time.time()
    log_padding = len(str(f"{total_samples:,}"))
    zfpkm_series: list[pd.Series] = []
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
            zfpkm_series.append(result.zfpkm)

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


def zfpkm_plot(results: dict[str, _ZFPKMResult], *, output_png_dirpath: Path, plot_xfloor: int = -4, subplot_titles: bool = True) -> None:
    """Plot the log2(FPKM) density and fitted Gaussian for each sample.

    Args:
        results: A dictionary of intermediate results from zfpkm_transform.
        output_png_dirpath: Output directory location
        subplot_titles: Whether to display facet titles (sample names).
        plot_xfloor: Lower limit for the x-axis.
        subplot_titles: Whether to display facet titles (sample names).

    """
    to_concat: list[pd.DataFrame] = [None] * len(results)  # type: ignore  # ignoring because None is not of type pd.DataFrame
    for i, name, result in enumerate(results.items()):
        stddev: float = float(result.std_dev)
        x: npt.NDArray[float] = result.density.x.flatten()
        y: npt.NDArray[float] = result.density.y.flatten()

        fitted: npt.NDArray[float] = np.exp(-0.5 * ((x - result.mu) / stddev) ** 2) / (stddev * np.sqrt(2 * np.pi))
        max_fpkm: float = float(y.max())
        max_fitted: float = float(fitted.max())
        scale_fitted: float = fitted * max_fpkm / max_fitted
        to_concat[i] = pd.DataFrame({"sample_name": name, "log2fpkm": x, "fpkm_density": y, "fitted_density_scaled": scale_fitted})

    mega_df = pd.concat(to_concat, ignore_index=True)
    mega_df.columns = pd.Series(data=["sample_name", "log2fpkm", "fpkm_density", "fitted_density_scaled"])
    mega_df = mega_df.melt(id_vars=["log2fpkm", "sample_name"], var_name="source", value_name="density")

    axes: list[plt.Axes]
    _, axes = plt.subplots(nrows=len(results), ncols=1, figsize=(8, 4 * len(results)))
    if len(results) == 1:
        axes = [axes]

    for i, sample_name in enumerate(results):
        sample_data = mega_df[mega_df["sample_name"] == sample_name]
        axis = axes[i]

        for source_type in sample_data["source"].unique():
            group = sample_data[sample_data["source"] == source_type]
            sns.lineplot(data=group, x="log2fpkm", y="density", label=source_type, ax=axis)

        if subplot_titles:
            axis.set_title(f"Sample: {sample_name}")
        axis.set_xlim(plot_xfloor, sample_data["log2fpkm"].max())
        axis.set_xlabel("log2(FPKM)")
        axis.set_ylabel("density [scaled]")
        axis.legend(title="Source")

    output_png_dirpath.mkdir(parents=True, exist_ok=True)
    sample_name: str = next(iter(results.keys()))[:-2]  # Go from 'control1hr_S1R1' to 'control1hr_S1'
    plt.tight_layout()
    plt.savefig(Path(output_png_dirpath, f"{sample_name}_zfpkm_density.png"))


def calculate_z_score(metrics: NamedMetrics) -> NamedMetrics:
    """Calculate the z-score for each sample in the metrics dictionary.

    Args:
        metrics: A dictionary of study metrics to calculate z-scores for.

    Returns:
        A dictionary of study metrics with z-scores calculated.
    """
    for sample in metrics:
        log_matrix = np.log(metrics[sample].normalization_matrix)
        z_matrix = pd.DataFrame(data=sklearn.preprocessing.scale(log_matrix, axis=1), columns=metrics[sample].sample_names)
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
        library_size: pd.DataFrame = counts.sum(axis=1)

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
            median_sum = np.float64(np.median(np.sum(counts[:, i])))
            if cut_off == "default":  # noqa: SIM108
                cutoff = np.float64(10e6) / median_sum
            else:
                cutoff = np.float64(1e6 * cut_off) / median_sum
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
        quantile_cutoff = np.quantile(a=tpm_quantile.values, q=1 - (cut_off / 100), axis=0)  # Compute quantile across columns
        boolean_expression = pd.DataFrame(data=tpm_matrix > quantile_cutoff, index=tpm_matrix.index, columns=tpm_matrix.columns).astype(int)

        min_func = k_over_a(min_samples, 0.9)
        top_func = k_over_a(top_samples, 0.9)

        min_genes: npt.NDArray[bool] = genefilter(boolean_expression, min_func)
        top_genes: npt.NDArray[bool] = genefilter(boolean_expression, top_func)

        # Only keep `entrez_gene_ids` that pass `min_genes`
        metric.entrez_gene_ids = [gene for gene, keep in zip(entrez_ids, min_genes, strict=True) if keep]
        metric.gene_sizes = np.array(gene for gene, keep in zip(gene_size, min_genes, strict=True) if keep)
        metric.count_matrix = metric.count_matrix.iloc[min_genes, :]
        metric.normalization_matrix = metrics[sample].normalization_matrix.iloc[min_genes, :]

        keep_top_genes = [gene for gene, keep in zip(entrez_ids, top_genes, strict=True) if keep]
        metric.high_confidence_entrez_gene_ids = [gene for gene, keep in zip(entrez_ids, keep_top_genes, strict=True) if keep]

    metrics = calculate_z_score(metrics)

    return metrics


def zfpkm_filter(
    *,
    metrics: NamedMetrics,
    filtering_options: _FilteringOptions,
    calculate_fpkm: bool,
    force_zfpkm_plot: bool,
    peak_parameters: PeakIdentificationParameters,
    bandwidth: float,
    output_png_dirpath: Path | None,
) -> NamedMetrics:
    """Apply zFPKM filtering to the FPKM matrix for a given sample.

    Args:
        metrics: A dictionary of study metrics to filter.
        filtering_options: Options for filtering the count matrix.
        calculate_fpkm: Whether to calculate FPKM from counts.
        force_zfpkm_plot: Whether to force plotting of zFPKM results even if there are many samples.
        peak_parameters: Parameters for peak identification in zFPKM calculation.
        bandwidth: The bandwidth for kernel density estimation in zFPKM calculation.
        output_png_dirpath: Optional directory path to save the zFPKM plot.

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
        matrix = matrix[matrix.sum(axis=1) > 0]  # remove rows (genes) that have no counts across all samples

        results, zfpkm_df = zfpkm_transform(matrix, peak_parameters=peak_parameters, bandwidth=bandwidth)
        zfpkm_df[(matrix == 0) | (zfpkm_df.isna())] = -4

        if len(results) > 10 and not force_zfpkm_plot:
            logger.warning(
                "Not plotting zFPKM results because more than 10 plots would be created. "
                "If you would like to plot them anyway, set 'force_zfpkm_plot' to True"
            )
        elif output_png_dirpath is None:
            logger.critical("Output zFPKM PNG filepath is None, set a path to plot zFPKM graphs")
        else:
            zfpkm_plot(results, output_png_dirpath=output_png_dirpath)

        metric.z_score_matrix = zfpkm_df

        # determine which genes are expressed
        min_samples = round(min_sample_expression * len(zfpkm_df.columns))
        min_func = k_over_a(min_samples, cut_off)
        min_genes: npt.NDArray[bool] = genefilter(zfpkm_df, min_func)
        metric.entrez_gene_ids = [gene for gene, keep in zip(zfpkm_df.index, min_genes, strict=True) if keep]

        # determine which genes are confidently expressed
        top_samples = round(high_confidence_sample_expression * len(zfpkm_df.columns))
        top_func = k_over_a(top_samples, cut_off)
        top_genes: npt.NDArray[bool] = genefilter(zfpkm_df, top_func)
        metric.high_confidence_entrez_gene_ids = [gene for gene, keep in zip(zfpkm_df.index, top_genes, strict=True) if keep]

    return metrics


def filter_counts(
    *,
    context_name: str,
    metrics: NamedMetrics,
    technique: FilteringTechnique,
    filtering_options: _FilteringOptions,
    prep: RNAType,
    force_zfpkm_plot: bool,
    peak_parameters: PeakIdentificationParameters,
    bandwidth: float,
    output_zfpkm_plot_dirpath: Path | None = None,
) -> NamedMetrics:
    """Filter the count matrix based on the specified technique.

    Args:
        context_name: The name of the context being processed.
        metrics: A dictionary of study metrics to filter.
        technique: The filtering technique to use.
        filtering_options: Options for filtering the count matrix.
        prep: The RNA preparation type.
        force_zfpkm_plot: Whether to force plotting of zFPKM results even if there are many samples.
        peak_parameters: Parameters for peak identification in zFPKM calculation.
        bandwidth: The bandwidth for kernel density estimation in zFPKM calculation.
        output_zfpkm_plot_dirpath: Optional directory path to save the zFPKM plot.

    Returns:
        A dictionary of filtered study metrics.
    """
    match technique:
        case FilteringTechnique.CPM:
            return cpm_filter(context_name=context_name, metrics=metrics, filtering_options=filtering_options, prep=prep)
        case FilteringTechnique.TPM:
            return tpm_quantile_filter(metrics=metrics, filtering_options=filtering_options)
        case FilteringTechnique.ZFPKM:
            return zfpkm_filter(
                metrics=metrics,
                filtering_options=filtering_options,
                calculate_fpkm=True,
                force_zfpkm_plot=force_zfpkm_plot,
                peak_parameters=peak_parameters,
                bandwidth=bandwidth,
                output_png_dirpath=output_zfpkm_plot_dirpath,
            )
        case FilteringTechnique.UMI:
            # UMI filtering is the same as zFPKM filtering without calculating FPKM
            return zfpkm_filter(
                metrics=metrics,
                filtering_options=filtering_options,
                calculate_fpkm=False,
                force_zfpkm_plot=force_zfpkm_plot,
                peak_parameters=peak_parameters,
                bandwidth=bandwidth,
                output_png_dirpath=output_zfpkm_plot_dirpath,
            )
        case _:
            _log_and_raise_error(
                f"Technique must be one of {FilteringTechnique}, got '{technique.value}'",
                error=ValueError,
                level=LogLevel.ERROR,
            )


async def _process(
    context_name: str,
    rnaseq_matrix_filepath: Path,
    metadata_df: pd.DataFrame,
    gene_info_df: pd.DataFrame,
    prep: RNAType,
    taxon: int,
    replicate_ratio: float,
    batch_ratio: float,
    high_replicate_ratio: float,
    high_batch_ratio: float,
    technique: FilteringTechnique,
    cut_off: int | float,
    force_zfpkm_plot: bool,
    peak_parameters: PeakIdentificationParameters,
    bandwidth: float,
    output_boolean_activity_filepath: Path,
    output_zscore_normalization_filepath: Path,
    output_zfpkm_plot_dirpath: Path | None,
):
    """Save the results of the RNA-Seq tests to a CSV file."""
    output_boolean_activity_filepath.parent.mkdir(parents=True, exist_ok=True)

    rnaseq_matrix: pd.DataFrame = await _read_file(rnaseq_matrix_filepath, h5ad_as_df=True)

    if rnaseq_matrix_filepath.suffix == ".h5ad":
        conversion = await gene_symbol_to_ensembl_and_gene_id(symbols=rnaseq_matrix["gene_symbol"].tolist(), taxon=taxon)
        conversion.reset_index(inplace=True)
        rnaseq_matrix = rnaseq_matrix.merge(conversion, how="left", on="gene_symbol")
        rnaseq_matrix.replace(to_replace=pd.NA, value="-")

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

    metrics: NamedMetrics = filter_counts(
        context_name=context_name,
        metrics=metrics,
        technique=technique,
        filtering_options=filtering_options,
        prep=prep,
        force_zfpkm_plot=force_zfpkm_plot,
        peak_parameters=peak_parameters,
        bandwidth=bandwidth,
        output_zfpkm_plot_dirpath=output_zfpkm_plot_dirpath,
    )

    merged_zscore_df = pd.DataFrame()
    expressed_genes: list[str] = []
    top_genes: list[str] = []
    for metric in metrics.values():
        expressed_genes.extend(metric.entrez_gene_ids)
        top_genes.extend(metric.high_confidence_entrez_gene_ids)

        merged_zscore_df = (
            metric.z_score_matrix
            if merged_zscore_df.empty
            else merged_zscore_df.merge(
                metric.z_score_matrix,
                how="outer",
                left_index=True,
                right_index=True,
            )
        )
    merged_zscore_df[merged_zscore_df.isna()] = -4

    # If any of the normalization metrics are not empty, write the normalized metrics to disk
    if not all(metric.normalization_matrix.empty for metric in metrics.values()):
        merged_zscore_df.to_csv(output_zscore_normalization_filepath, index=True)
        logger.success(f"Wrote z-score normalization matrix to {output_zscore_normalization_filepath}")
    else:
        logger.warning(
            "Not writing z-score normalization matrix because no normalization matrices exist. This is expected if you are using UMI filtering."
        )

    expression_frequency = pd.Series(expressed_genes).value_counts()
    expression_df = pd.DataFrame({"entrez_gene_id": expression_frequency.index, "frequency": expression_frequency.values})
    expression_df["prop"] = expression_df["frequency"] / len(metrics)
    expression_df = expression_df[expression_df["prop"] >= filtering_options.batch_ratio]

    top_frequency = pd.Series(top_genes).value_counts()
    top_df = pd.DataFrame({"entrez_gene_id": top_frequency.index, "frequency": top_frequency.values})
    top_df["prop"] = top_df["frequency"] / len(metrics)
    top_df = top_df[top_df["prop"] >= filtering_options.high_batch_ratio]

    entrez_id_series = pd.Series(read_counts_results.entrez_gene_ids)
    boolean_matrix = pd.DataFrame(
        data={
            "entrez_gene_id": read_counts_results.entrez_gene_ids,
            "expressed": entrez_id_series.isin(expression_df["entrez_gene_id"]).astype(int),
            "high": entrez_id_series.isin(top_df["entrez_gene_id"]).astype(int),
        }
    )

    expressed_count = len(boolean_matrix[boolean_matrix["expressed"] == 1])
    high_confidence_count = len(boolean_matrix[boolean_matrix["high"] == 1])

    boolean_matrix.dropna(subset="entrez_gene_id", inplace=True)
    boolean_matrix.to_csv(output_boolean_activity_filepath, index=False)
    logger.info(f"{context_name} - Found {expressed_count} expressed genes, {high_confidence_count} of which are confidently expressed")
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
    zfpkm_peak_height: float = 0.02,
    zfpkm_peak_distance: float = 1.0,
    zfpkm_bandwidth: float = 1.0,
    cutoff: int | float | None = None,
    force_zfpkm_plot: bool = False,
    log_level: LogLevel = LogLevel.INFO,
    log_location: str | TextIO = sys.stderr,
    output_zfpkm_plot_dirpath: Path | None = None,
) -> None:
    """Generate a list of active and high-confidence genes from a gene count matrix.

    Replicates are compared for consensus within the study/batch number according to replicate ratios,
        then study/batch numbers are checked for consensus according to batch ratios.
    The zFPKM method is outlined here: https://pubmed.ncbi.nlm.nih.gov/24215113/

    Arg:
        context_name: The name of the context being processed
        input_rnaseq_filepath: The filepath to the gene count matrix
        input_gene_info_filepath: The filepath to the gene info file
        output_boolean_activity_filepath: The filepath to write the output gene count matrix
        output_zscore_normalization_filepath: The filepath to write the output z-score normalization matrix
        prep: The preparation method
        taxon_id: The NCBI Taxon ID
        input_metadata_filepath_or_df: The filepath or dataframe containing metadata information
        replicate_ratio: The percentage of replicates that a gene must appear in for a gene to be marked as "active" in a batch/study
        batch_ratio: The percentage of batches that a gene must appear in for a gene to be marked as 'active"
        high_replicate_ratio: The percentage of replicates that a gene must appear in for a gene to be
            marked "highly confident" in its expression in a batch/study
        high_batch_ratio: The percentage of batches that a gene must appear in for a gene to be marked "highly confident" in its expression
        technique: The filtering technique to use
        zfpkm_peak_height: The height of the zFPKM peak
        zfpkm_peak_distance: The distance of the zFPKM peak
        zfpkm_bandwidth: The bandwidth of the zFPKM
        cutoff: The cutoff value to use for the provided filtering technique
        force_zfpkm_plot: If too many samples exist, should plotting be done anyway?
        log_level: The level of logging to output
        log_location: The location to write logs to
        output_zfpkm_plot_dirpath: Optional filepath to save zFPKM plots

    """
    _set_up_logging(level=log_level, location=log_location)

    technique = FilteringTechnique(technique) if isinstance(technique, str) else technique
    match technique:
        case FilteringTechnique.TPM:
            cutoff: int | float = cutoff or 25
            if cutoff < 1 or cutoff > 100:
                _log_and_raise_error(
                    "Quantile must be between 1 - 100",
                    error=ValueError,
                    level=LogLevel.ERROR,
                )

        case FilteringTechnique.CPM:
            if cutoff and cutoff < 0:
                _log_and_raise_error(
                    "Cutoff must be greater than or equal to 0",
                    error=ValueError,
                    level=LogLevel.ERROR,
                )
            elif cutoff:
                cutoff = "default"

        case FilteringTechnique.ZFPKM | FilteringTechnique.UMI:
            cutoff: int | float = cutoff or -3
        case _:
            _log_and_raise_error(
                f"Technique must be one of {','.join(FilteringTechnique)}. Got: {technique.value}",
                error=ValueError,
                level=LogLevel.ERROR,
            )

    if not input_rnaseq_filepath.exists():
        _log_and_raise_error(
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
            _log_and_raise_error(
                f"Expected an excel file with extension of '.xlsx' or '.xls', got '{input_metadata_filepath_or_df.suffix}'.",
                error=ValueError,
                level=LogLevel.ERROR,
            )
        if not input_metadata_filepath_or_df.exists():
            _log_and_raise_error(
                f"Input metadata file not found! Searching for: '{input_metadata_filepath_or_df}'",
                error=FileNotFoundError,
                level=LogLevel.ERROR,
            )

        metadata_df = pd.read_excel(input_metadata_filepath_or_df)
    else:
        _log_and_raise_error(
            f"Expected a pandas DataFrame or Path object as metadata, got '{type(input_metadata_filepath_or_df)}'",
            error=TypeError,
            level=LogLevel.ERROR,
        )

    # metadata_df["fragment_length"] = metadata_df["fragment_length"].astype(np.float32)
    # metadata_df = metadata_df.groupby("sample_name", as_index=False).agg(
    #     {
    #         "sample_name": "first",
    #         "fragment_length": "mean",
    #         "layout": "first",
    #         "strand": "first",
    #         "study": "first",
    #         "library_prep": "first",
    #     }
    # )
    logger.debug(f"Starting '{context_name}'")
    await _process(
        context_name=context_name,
        rnaseq_matrix_filepath=input_rnaseq_filepath,
        metadata_df=metadata_df,
        gene_info_df=await _read_file(input_gene_info_filepath),
        prep=prep,
        taxon=taxon_id,
        replicate_ratio=replicate_ratio,
        batch_ratio=batch_ratio,
        high_replicate_ratio=high_replicate_ratio,
        high_batch_ratio=high_batch_ratio,
        technique=technique,
        cut_off=cutoff,
        force_zfpkm_plot=force_zfpkm_plot,
        peak_parameters=PeakIdentificationParameters(height=zfpkm_peak_height, distance=zfpkm_peak_distance),
        bandwidth=zfpkm_bandwidth,
        output_boolean_activity_filepath=output_boolean_activity_filepath,
        output_zscore_normalization_filepath=output_zscore_normalization_filepath,
        output_zfpkm_plot_dirpath=output_zfpkm_plot_dirpath,
    )
