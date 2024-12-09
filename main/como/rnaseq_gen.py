from __future__ import annotations

import math
import multiprocessing
import time
from collections import namedtuple
from dataclasses import dataclass, field
from enum import Enum
from functools import partial
from multiprocessing.pool import Pool
from pathlib import Path
from typing import Callable, NamedTuple

import numpy as np
import numpy.typing as npt
import pandas as pd
import plotly.graph_objs as go
import scanpy as sc
import sklearn
import sklearn.neighbors
from fast_bioservices.pipeline import ensembl_to_gene_id_and_symbol
from loguru import logger
from pandas import DataFrame
from plotly.subplots import make_subplots
from scipy.signal import find_peaks
from sklearn.neighbors import KernelDensity

from como.migrations import gene_info_migrations
from como.project import Config
from como.types import RNAPrepMethod


class _FilteringOptions(NamedTuple):
    replicate_ratio: float
    batch_ratio: float
    cut_off: float
    high_replicate_ratio: float
    high_batch_ratio: float


class FilteringTechnique(Enum):
    """RNA sequencing filtering capabilities."""

    cpm = "cpm"
    zfpkm = "zfpkm"
    tpm = "quantile"
    umi = "umi"

    @staticmethod
    def from_string(value: str) -> FilteringTechnique:
        """Create a filtering technique object from a string."""
        match value.lower():
            case "cpm":
                return FilteringTechnique.cpm
            case "zfpkm":
                return FilteringTechnique.zfpkm
            case "quantile":
                return FilteringTechnique.tpm
            case "umi":
                return FilteringTechnique.umi
            case _:
                possible_values = [t.value for t in FilteringTechnique]
                raise ValueError(f"Got a filtering technique of '{value}'; should be one of: {possible_values}")


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
                raise ValueError(f"Layout must be 'paired-end' or 'single-end'; got: {layout}")

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
        raise TypeError("Unsupported data type. Must be a Pandas DataFrame or a NumPy array.")

    return (
        data.apply(filter_func, axis=1).values
        if isinstance(data, pd.DataFrame)
        else np.apply_along_axis(filter_func, axis=1, arr=data)
    )


async def _read_counts(path: Path) -> pd.DataFrame:
    if path.suffix not in {".csv", ".h5ad"}:
        raise ValueError(f"Unknown file extension '{path.suffix}'. Valid options are '.csv' or '.h5ad'.")

    matrix: pd.DataFrame
    if path.suffix == ".csv":
        logger.debug(f"Reading CSV file at '{path}'")
        matrix = pd.read_csv(path, header=0)
    elif path.suffix == ".h5ad":
        logger.debug(f"Reading h5ad file at '{path}'")
        # Make sample names the columns and gene data the index
        matrix = sc.read_h5ad(path).to_df().T

    return matrix


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
    gene_info = gene_info_migrations(gene_info)
    conversion = await ensembl_to_gene_id_and_symbol(ids=matrix["ensembl_gene_id"].tolist(), taxon=taxon)
    matrix = matrix.merge(conversion, on="ensembl_gene_id", how="left")

    # Only include Entrez and Ensembl Gene IDs that are present in `gene_info`
    matrix["entrez_gene_id"] = matrix["entrez_gene_id"].str.split("//")
    matrix = matrix.explode("entrez_gene_id")
    matrix = matrix.replace(to_replace="-", value=pd.NA).dropna()
    matrix["entrez_gene_id"] = matrix["entrez_gene_id"].astype(int)

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


def calculate_fpkm(metrics: NamedMetrics) -> NamedMetrics:
    """Calculate the Fragments Per Kilobase of transcript per Million mapped reads (FPKM) for each sample in the metrics dictionary."""  # noqa: E501
    matrix_values = []
    for study in metrics:
        for sample in range(metrics[study].num_samples):
            layout = metrics[study].layout[sample]
            count_matrix: npt.NDArray = metrics[study].count_matrix.iloc[:, sample].values
            gene_size = metrics[study].gene_sizes

            count_matrix = count_matrix.astype(np.float32)
            gene_size = gene_size.astype(np.float32)

            match layout:
                case LayoutMethod.paired_end:  # FPKM
                    mean_fragment_lengths = metrics[study].fragment_lengths[sample]
                    # Ensure non-negative value
                    effective_length = [max(0, size - (mean_fragment_lengths + 1)) for size in gene_size]
                    n = count_matrix.sum()
                    fpkm = ((count_matrix + 1) * 1e9) / (np.array(effective_length) * n)
                    matrix_values.append(fpkm)
                case LayoutMethod.single_end:  # RPKM
                    # Add a pseudocount before log to ensure log(0) does not happen
                    rate = np.log(count_matrix + 1) - np.log(gene_size)
                    exp_rate = np.exp(rate - np.log(np.sum(count_matrix)) + np.log(1e9))
                    matrix_values.append(exp_rate)
                case _:
                    raise ValueError("Invalid normalization method specified")

        fpkm_matrix = pd.DataFrame(matrix_values).T  # Transpose is needed because values were appended as rows
        fpkm_matrix = fpkm_matrix[~pd.isna(fpkm_matrix)]
        metrics[study].normalization_matrix = fpkm_matrix

        metrics[study].normalization_matrix.columns = metrics[study].count_matrix.columns

    return metrics


def _zfpkm_calculation(col: pd.Series, kernel: KernelDensity, peak_parameters: tuple[float, float]) -> _ZFPKMResult:
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
            a threshold for calling a gene as "expressed"
            : https://doi.org/10.1186/1471-2164-14-778
    """
    col_log2: npt.NDArray = np.log2(col + 1)
    col_log2 = np.nan_to_num(col_log2, nan=0)
    refit: KernelDensity = kernel.fit(col_log2.reshape(-1, 1))  # type: ignore

    # kde: KernelDensity = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(col_log2.reshape(-1, 1))
    x_range = np.linspace(col_log2.min(), col_log2.max(), 1000)
    density = np.exp(refit.score_samples(x_range.reshape(-1, 1)))
    peaks, _ = find_peaks(density, height=peak_parameters[0], distance=peak_parameters[1])
    peak_positions = x_range[peaks]

    mu = 0
    max_fpkm = 0
    stddev = 1

    if len(peaks) != 0:
        mu = peak_positions.max()
        max_fpkm = density[peaks[np.argmax(peak_positions)]]
        u = col_log2[col_log2 > mu].mean()
        stddev = (u - mu) * np.sqrt(np.pi / 2)
    zfpkm = pd.Series((col_log2 - mu) / stddev, dtype=np.float32, name=col.name)

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

    total = len(fpkm_df.columns)
    update_per_step: int = int(np.ceil(total * update_every_percent))
    cores = multiprocessing.cpu_count() - 2
    logger.debug(f"Processing {total:,} samples through zFPKM transform using {cores} cores")
    logger.debug(
        f"Will update every {update_per_step:,} steps as this is approximately "
        f"{update_every_percent:.1%} of {total:,}"
    )

    with Pool(processes=cores) as pool:
        kernel = KernelDensity(kernel="gaussian", bandwidth=bandwidth)
        chunksize = int(math.ceil(len(fpkm_df.columns) / (4 * cores)))
        partial_func = partial(_zfpkm_calculation, kernel=kernel, peak_parameters=peak_parameters)
        chunk_time = time.time()
        start_time = time.time()

        log_padding = len(str(f"{total:,}"))
        zfpkm_df = pd.DataFrame(data=0, index=fpkm_df.index, columns=fpkm_df.columns)
        results: dict[str, _ZFPKMResult] = {}
        result: _ZFPKMResult
        for i, result in enumerate(
            pool.imap(
                partial_func,
                (fpkm_df[col] for col in fpkm_df.columns),
                chunksize=chunksize,
            )
        ):
            key = str(result.zfpkm.name)
            results[key] = result
            zfpkm_df[key] = result.zfpkm

            # show updates every X% and at the end, but skip on first iteration
            if i != 0 and (i % update_per_step == 0 or i == total):
                current_time = time.time()
                chunk = current_time - chunk_time
                total_time = current_time - start_time
                formatted = f"{i:,}"
                logger.debug(
                    f"Processed {formatted:>{log_padding}} of {total:,} - "
                    f"chunk took {chunk:.1f} seconds - "
                    f"running for {total_time:.1f} seconds"
                )
                chunk_time = current_time
    return results, zfpkm_df


def zfpkm_plot(results, *, plot_xfloor: int = -4, subplot_titles: bool = True):
    """Plot the log2(FPKM) density and fitted Gaussian for each sample.

    :param results: A dictionary of intermediate results from zfpkm_transform.
    :param: subplot_titles: Whether to display facet titles (sample names).
    :param plot_xfloor: Lower limit for the x-axis.
    :param subplot_titles: Whether to display facet titles (sample names).
    """
    mega_df = pd.DataFrame(columns=["sample_name", "log2fpkm", "fpkm_density", "fitted_density_scaled"])
    for name, result in results.items():
        stddev = result.std_dev
        x = np.array(result.density.x)
        y = np.array(result.density.y)

        fitted = np.exp(-0.5 * ((x - result.mu) / stddev) ** 2) / (stddev * np.sqrt(2 * np.pi))
        max_fpkm = y.max()
        max_fitted = fitted.max()
        scale_fitted = fitted * (max_fpkm / max_fitted)

        df = pd.DataFrame(
            {
                "sample_name": [name] * len(x),
                "log2fpkm": x,
                "fpkm_density": y,
                "fitted_density_scaled": scale_fitted,
            }
        )
        mega_df = pd.concat([mega_df, df], ignore_index=True)

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



def zfpkm_filter(*, metrics: NamedMetrics, filtering_options: _FilteringOptions, calcualte_fpkm: bool) -> NamedMetrics:
    """Apply zFPKM filtering to the FPKM matrix for a given sample."""
    min_sample_expression = filtering_options.replicate_ratio
    high_confidence_sample_expression = filtering_options.high_replicate_ratio
    cut_off = filtering_options.cut_off

    if calcualte_fpkm:
        metrics = calculate_fpkm(metrics)

    metric: _StudyMetrics
    for metric in metrics.values():
        # if fpkm was not calculated, the normalization matrix will be empty; collect the count matrix instead
        matrix = metric.count_matrix if metric.normalization_matrix.empty else metric.normalization_matrix
        matrix = matrix[matrix.sum(axis=1) > 0]

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


async def rnaseq_gen(
    # config_filepath: Path,
    config_filename: str,
    prep: RNAPrepMethod,
    taxon_id: int | str | Taxon,
    replicate_ratio: float = 0.5,
    high_replicate_ratio: float = 1.0,
    batch_ratio: float = 0.5,
    high_batch_ratio: float = 1.0,
    technique: FilteringTechnique | str = FilteringTechnique.tpm,
    cut_off: int | float | None = None,
) -> None:
    """Generate a list of active and high-confidence genes from a gene count matrix.

    Replicates are compared for consensus within the study/batch number according to replicate ratios,
        then study/batch numbers are checked for consensus according to batch ratios.
    The zFPKM method is outlined here: https://pubmed.ncbi.nlm.nih.gov/24215113/

    :param config_filename: The configuration filename to read
    :param prep: The preparation method
    :param taxon_id: The NCBI Taxon ID
    :param replicate_ratio: The percentage of replicates that a gene must
        appear in for a gene to be marked as "active" in a batch/study
    :param batch_ratio: The percentage of batches that a gene must appear in for a gene to be marked as 'active"
    :param high_replicate_ratio: The percentage of replicates that a gene must
        appear in for a gene to be marked "highly confident" in its expression in a batch/study
    :param high_batch_ratio: The percentage of batches that a gene must
        appear in for a gene to be marked "highly confident" in its expression
    :param technique: The filtering technique to use
    :param cut_off: The cutoff value to use for the provided filtering technique
    :return: None
    """
    if isinstance(technique, str):
        technique = FilteringTechnique(technique.lower())
    if isinstance(taxon_id, (str, int)):
        taxon_id = Taxon.from_string(str(taxon_id))

    match technique:
        case FilteringTechnique.tpm:
            cut_off = 25 if cut_off is None else cut_off
            if cut_off < 1 or cut_off > 100:
                raise ValueError("Quantile must be between 1 - 100")

        case FilteringTechnique.cpm:
            if cut_off is not None and cut_off < 0:
                raise ValueError("Cutoff must be greater than 0")
            elif cut_off is None:
                cut_off = "default"

        case FilteringTechnique.zfpkm:
            cut_off = "default" if cut_off is None else cut_off
        case FilteringTechnique.umi:
            pass
        case _:
            raise ValueError(f"Technique must be one of {FilteringTechnique}")

    await _handle_context_batch(
        config_filename=config_filename,
        replicate_ratio=replicate_ratio,
        replicate_ratio_high=high_replicate_ratio,
        batch_ratio=batch_ratio,
        batch_ratio_high=high_batch_ratio,
        technique=technique,
        cut_off=cut_off,
        prep=prep,
        taxon=taxon_id,
    )


    )
