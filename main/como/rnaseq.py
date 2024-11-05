from collections import namedtuple
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Callable, NamedTuple

import numpy as np
import numpy.typing as npt
import pandas as pd
import plotly.graph_objs as go
import scanpy as sc
import sklearn
import sklearn.neighbors
from fast_bioservices import BioDBNet, Output, Taxon
from loguru import logger
from plotly.subplots import make_subplots
from scipy.signal import find_peaks
from scipy.stats import norm
from sklearn.neighbors import KernelDensity

from como.custom_types import RNASeqPreparationMethod
from como.project import Config


class _FilteringOptions(NamedTuple):
    replicate_ratio: float
    batch_ratio: float
    cut_off: float
    high_replicate_ratio: float
    high_batch_ratio: float


class FilteringTechnique(Enum):
    cpm = "cpm"
    zfpkm = "zfpkm"
    tpm = "quantile"
    umi = "umi"

    @staticmethod
    def from_string(value: str) -> "FilteringTechnique":
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
                raise ValueError(f"Filtering technique must be one of {possible_values}; got: {value}")


@dataclass
class _StudyMetrics:
    study: str
    num_samples: int
    count_matrix: pd.DataFrame
    fragment_lengths: npt.NDArray[np.float32]
    sample_names: list[str]
    layout: list[str]
    entrez_gene_ids: list[str]
    gene_sizes: npt.NDArray[np.float32]
    __normalization_matrix: pd.DataFrame = field(default_factory=pd.DataFrame)
    __z_score_matrix: pd.DataFrame = field(default_factory=pd.DataFrame)
    __high_confidence_entrez_gene_ids: list[str] = field(default=list)

    def __post_init__(self):
        for layout in self.layout:
            if layout not in ["paired-end", "single-end", ""]:
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


Density = namedtuple("Density", ["x", "y"])


class zFPKMresult(NamedTuple):
    zfpkm: pd.Series
    density: Density
    mu: float
    std_dev: float
    max_fpkm: float


class ReadMatrixResults(NamedTuple):
    metrics: dict[str, _StudyMetrics]
    entrez_gene_ids: list[str]


NamedMetrics = dict[str, _StudyMetrics]


def k_over_a(k: int, A: float) -> Callable[[npt.NDArray], bool]:
    """
    Return a function that filters rows of an array based on the sum of elements being greater than or equal to A at least k times.

    This code is based on the `kOverA` function found in R's `genefilter` package: https://www.rdocumentation.org/packages/genefilter/versions/1.54.2/topics/kOverA

    :param k: The minimum number of times the sum of elements must be greater than or equal to A.
    :param A: The threshold value.
    :return: A function that accepts a NumPy array to perform the actual filtering
    """

    def filter_func(row: npt.NDArray) -> bool:
        return np.sum(row >= A) >= k

    return filter_func


def genefilter(data: pd.DataFrame | npt.NDArray, filter_func: Callable[[npt.NDArray], bool]) -> npt.NDArray:
    """
    Apply a filter function to the rows of the data and return the filtered array.

    This code is based on the `genefilter` function found in R's `genefilter` package: https://www.rdocumentation.org/packages/genefilter/versions/1.54.2/topics/genefilter

    :param data: The data to filter
    :param filter_func: THe function to filter the data by
    :return: A NumPy array of the filtered data.
    """
    match type(data):
        case pd.DataFrame:
            return data.apply(filter_func, axis=1).values
        case npt.NDArray:
            return np.apply_along_axis(filter_func, axis=1, arr=data)
        case _:
            raise ValueError("Unsupported data type. Must be a Pandas DataFrame or a NumPy array.")


async def _read_counts_matrix(
    *,
    biodbnet: BioDBNet,
    context_name: str,
    counts_matrix_filepath: Path,
    config_filepath: Path,
    gene_info_filepath: Path,
    taxon_id: Taxon,
) -> ReadMatrixResults:
    """
    Reads the counts matrix and returns the results.

    :param biodbnet: The BioDBNet object to use for determining the input gene data type
    :param context_name: The context name being processed. Usually a cell type, but can be any string
    :param counts_matrix_filepath: The file path to the gene count matrix
    :param config_filepath: The file path to the Excel configuration file
    :param gene_info_filepath: The file path to gene information generated by `rnaseq_preprocess.py`
    :param taxon_id: The NCBI Taxon ID
    :return: A dataclass `ReadMatrixResults`
    """
    logger.trace(f"Reading config_filepath at '{config_filepath}'")
    config_df: pd.DataFrame = pd.read_excel(config_filepath, sheet_name=context_name, header=0)
    gene_info: pd.DataFrame = pd.read_csv(gene_info_filepath)
    gene_info = gene_info[gene_info["ensembl_gene_id"] != "-"].reset_index(drop=True)
    gene_sizes: npt.NDArray[np.float32] = np.array(gene_info["size"].values).astype(np.float32)

    match counts_matrix_filepath.suffix:
        case ".csv":
            logger.debug(f"Reading CSV file at '{counts_matrix_filepath}'")
            counts_matrix: pd.DataFrame = pd.read_csv(counts_matrix_filepath, header=0)
        case ".h5ad":
            logger.debug(f"Reading h5ad file at '{counts_matrix_filepath}'")
            adata: sc.AnnData = sc.read_h5ad(counts_matrix_filepath)
            counts_matrix: pd.DataFrame = adata.to_df().T  # Make sample names the columns and gene data the index

            # Coherce the incoming gene data (i.e., Gene Symbols) into Entrez and Ensembl Gene IDs
            cohersion = await biodbnet.db_find(
                values=counts_matrix.index.tolist(), output_db=[Output.GENE_ID, Output.ENSEMBL_GENE_ID], taxon=taxon_id
            )
            cohersion.rename(columns={"Gene ID": "entrez_gene_id", "Ensembl Gene ID": "ensembl_gene_id"}, inplace=True)

            # Set proper column names
            input_type = cohersion["Input Type"][0].replace(" ", "_").lower()
            cohersion.drop(columns=["Input Type"], inplace=True)
            cohersion.columns = pd.Index([input_type] + cohersion.columns.tolist()[1:])

            # Merge new gene data with counts_matrix
            counts_matrix.index.name = input_type
            counts_matrix.reset_index(inplace=True)
            counts_matrix = counts_matrix.merge(cohersion, on=input_type)
            counts_matrix = counts_matrix[counts_matrix["entrez_gene_id"] != "-"]

        case _:
            raise ValueError(f"Unknown file extension '{counts_matrix_filepath.suffix}' for file '{counts_matrix_filepath}'")

    # Only include Entrez and Ensembl Gene IDs that are present in `gene_info`
    counts_matrix["entrez_gene_id"] = counts_matrix["entrez_gene_id"].str.split("//")
    counts_matrix = counts_matrix.explode("entrez_gene_id")

    counts_matrix["entrez_gene_id"] = counts_matrix["entrez_gene_id"].astype(int)
    gene_info["entrez_gene_id"] = gene_info["entrez_gene_id"].astype(int)

    counts_matrix = counts_matrix.merge(gene_info[["entrez_gene_id", "ensembl_gene_id"]], on=["entrez_gene_id", "ensembl_gene_id"], how="inner")
    gene_info = gene_info.merge(counts_matrix[["entrez_gene_id", "ensembl_gene_id"]], on=["entrez_gene_id", "ensembl_gene_id"], how="inner")
    entrez_gene_ids: list[str] = gene_info["entrez_gene_id"].tolist()

    metrics: NamedMetrics = {}
    studies: list[str] = config_df["study"].unique().tolist()
    for study in studies:
        study_sample_names = config_df[config_df["study"] == study]["sample_name"].tolist()
        metrics[study] = _StudyMetrics(
            count_matrix=counts_matrix[counts_matrix.columns.intersection(study_sample_names)],
            fragment_lengths=config_df[config_df["study"] == study]["fragment_length"].values,
            sample_names=study_sample_names,
            layout=config_df[config_df["study"] == study]["layout"].tolist(),
            num_samples=len(study_sample_names),
            entrez_gene_ids=entrez_gene_ids,
            gene_sizes=gene_sizes,
            study=study,
        )
        metrics[study].fragment_lengths[np.isnan(metrics[study].fragment_lengths)] = 0
        metrics[study].count_matrix.index = pd.Index(entrez_gene_ids, name="entrez_gene_id")

    return ReadMatrixResults(metrics=metrics, entrez_gene_ids=gene_info["entrez_gene_id"].tolist())


def calculate_tpm(metrics: NamedMetrics) -> NamedMetrics:
    for sample in metrics.keys():
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
    matrix_values = []
    for study in metrics.keys():
        for sample in range(metrics[study].num_samples):
            layout = metrics[study].layout[sample]

            count_matrix: npt.NDArray = metrics[study].count_matrix.iloc[:, sample].values
            gene_size = metrics[study].gene_sizes

            count_matrix = count_matrix.astype(np.float32)
            gene_size = gene_size.astype(np.float32)

            match layout:
                case "paired-end":  # FPKM
                    mean_fragment_lengths = metrics[study].fragment_lengths[sample]
                    effective_length = [max(0, size - mean_fragment_lengths + 1) for size in gene_size]  # Ensure non-negative value
                    N = count_matrix.sum()
                    fpkm = (count_matrix + 1) * 1e9 / (np.array(effective_length) * N)
                    matrix_values.append(fpkm)
                case "single-end":  # RPKM
                    rate = np.log(count_matrix + 1) - np.log(gene_size)  # Add a pseudocount before log to ensure log(0) does not happen
                    exp_rate = np.exp(rate - np.log(np.sum(count_matrix)) + np.log(1e9))
                    matrix_values.append(exp_rate)
                case _:
                    raise ValueError("Invalid normalization method specified")

        fpkm_matrix = pd.DataFrame(matrix_values).T  # Transpose is needed because values were appended as rows
        fpkm_matrix = fpkm_matrix[~pd.isna(fpkm_matrix)]
        metrics[study].normalization_matrix = fpkm_matrix

        metrics[study].normalization_matrix.columns = metrics[study].count_matrix.columns

    return metrics


def zfpkm_transform(
    fpkm_df: pd.DataFrame,
    bandwidth: int = 0.5,
    peak_parameters: tuple[float, float] = (0.02, 1.0),
) -> tuple[dict[str, zFPKMresult], pd.DataFrame]:
    def calc(col: npt.NDArray, name: str, bandwidth: int, peak_parameters: tuple[float, float]) -> zFPKMresult:
        """
        Log2 Transformations
            - Stabilize the variance in the data to make the distribution more symmetric; this is helpful for Gaussian fitting

        Kernel Density Estimation (kde)
            - Non-parametric method to estimate the probability density function (PDF) of a random variable
            - Estimates the distribution of log2-transformed FPKM values
            - Bandwidth parameter controls the smoothness of the density estimate
            - KDE Explanation
                - A way to smooth a histogram to get a better idea of the underlying distribution of the data
                - Given a set of data points, we want to understand how they are distributed. Histograms can be useful, but are sensitive to bin size and number
                - The KDE places a "kernel", which is a small symmetric function (i.e., Gaussian curve), at each data point
                - The "kernel" acts as a weight, giving more weight to points closer to the center of the kernel and less weight to points farther away
                - Kernel functions are summed along each point on the x-axis
                - A smooth curve is created that represents the estimated density of the data

        Peak Finding
            - Identifies that are above a certain height and separated by a minimum distance
            - Represent potential local maxima in the distribution

        Peak Selection
            - The peak with the highest x-value (from log2-FPKM) is chosen as the mean (mu) of the "inactive" gene distribution
            - The peak representing unexpressed or inactive genes should be at a lower FPKM value compared to the peak representing expressed genes

        Standard Deviation Estimation
             - The mean of log2-FPKM values are greater than the calculated mu
             - Standard deviation is estimated based on the assumption that the right tail of the distribution (which represents expressed genes) can be approximated by a half-normal distribution

         zFPKM Transformation
            - Centers disbribution around 0 and scales it by the standard deviation, making it easier to compare gene expression across different samples
            - Represents the number of standard deviations away from the mean of the "inactive" gene distribution
            - Higher zFPKM values indicate higher expression levels relative to the "inactive" peak
            - A zFPKM value of 0 represents the mean of the "inactive" distribution
            - Research shows that a zFPKM value of -3 or greater can be used as a threshold for calling a gene as "expressed"
                : https://doi.org/10.1186/1471-2164-14-778
        """
        col_log2: npt.NDArray = np.log2(col + 1)
        col_log2 = np.nan_to_num(col_log2, nan=0)

        kde: KernelDensity = KernelDensity(kernel="gaussian", bandwidth=bandwidth).fit(col_log2.reshape(-1, 1))  # type: ignore
        x_range = np.linspace(col_log2.min(), col_log2.max(), 1000)
        density = np.exp(kde.score_samples(x_range.reshape(-1, 1)))
        peaks, _ = find_peaks(density, height=peak_parameters[0], distance=peak_parameters[1])
        peak_positions = x_range[peaks]

        mu = 0
        max_fpkm = 0
        u = 0
        stddev = 1

        if len(peaks) != 0:
            mu = peak_positions.max()
            max_fpkm = density[peaks[np.argmax(peak_positions)]]
            u = col_log2[col_log2 > mu].mean()
            stddev = (u - mu) * np.sqrt(np.pi / 2)
        zfpkm = pd.Series((col_log2 - mu) / stddev, dtype=np.float32, name=name)

        return zFPKMresult(zfpkm=zfpkm, density=Density(x_range, density), mu=mu, std_dev=stddev, max_fpkm=max_fpkm)

    zfpkm_df = pd.DataFrame(data=0, index=fpkm_df.index, columns=fpkm_df.columns)
    results = {}
    for col in fpkm_df.columns:
        results[col] = calc(fpkm_df[col], name=col, bandwidth=bandwidth, peak_parameters=peak_parameters)
        zfpkm_df[col] = results[col].zfpkm

    return results, zfpkm_df


def zfpkm_plot(results, *, plot_xfloor: int = -4, subplot_titles: bool = True):
    """
    Plot the log2(FPKM) density and fitted Gaussian for each sample.

    Args:
        results: A dictionary of intermediate results from zfpkm_transform.
        subplot_titles: Whether to display facet titles (sample names).
        plot_xfloor: Lower limit for the x-axis.
    """

    mega_df = pd.DataFrame()
    for name, result in results.items():
        d = result.density
        mu = result.mu
        stdev = result.std_dev
        max_fpkm = result.max_fpkm

        fitted = norm.pdf(d.x, loc=mu, scale=stdev)
        scale_fitted = fitted * (max_fpkm / fitted.max())

        df = pd.DataFrame({"sample_name": name, "log2fpkm": d.x, "fpkm_density": d.y, "fitted_density_scaled": scale_fitted})
        mega_df = pd.concat([mega_df, df], ignore_index=True)

    mega_df = mega_df.melt(
        id_vars=["log2fpkm", "sample_name"],
        value_vars=["fpkm_density", "fitted_density_scaled"],
        var_name="source",
        value_name="density",
    )

    subplot_titles = list(results.keys()) if subplot_titles else None
    fig = make_subplots(rows=len(results), cols=1, subplot_titles=subplot_titles, vertical_spacing=min(0.05, (1 / (len(results) - 1))))

    for i, (name, group) in enumerate(mega_df.groupby("sample_name"), start=1):
        fig.add_trace(
            trace=go.Scatter(x=group["log2fpkm"], y=group["density"], mode="lines", name=name, legendgroup=name),
            row=i,
            col=1,
        )
        fig.update_xaxes(title_text="log2(FPKM)", range=[plot_xfloor, group["log2fpkm"].max()], row=i, col=1)
        fig.update_yaxes(title_text="[scaled] density", row=i, col=1)
        fig.update_layout(legend_tracegroupgap=0)

    fig.update_layout(height=600 * len(results), width=1000, title_text="zFPKM Plots", showlegend=True)
    fig.write_image("zfpkm_plot.png")


def calculate_z_score(metrics: NamedMetrics) -> NamedMetrics:
    for sample in metrics.keys():
        log_matrix = np.log(metrics[sample].normalization_matrix)
        z_matrix = pd.DataFrame(data=sklearn.preprocessing.scale(log_matrix, axis=1), columns=metrics[sample].sample_names)
        metrics[sample].z_score_matrix = z_matrix
    return metrics


def cpm_filter(*, context_name: str, metrics: NamedMetrics, filtering_options: _FilteringOptions, prep: RNASeqPreparationMethod) -> NamedMetrics:
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
        # This will not impact the final counts per million calculation because the original counts are still 0; thus, (0 / 1) * 1_000_000 = 0
        library_size[library_size == 0] = 1

        output_filepath = config.result_dir / context_name / prep.value / f"CPM_Matrix_{prep.value}_{sample}.csv"
        output_filepath.parent.mkdir(parents=True, exist_ok=True)
        counts_per_million: pd.DataFrame = (counts / library_size) * 1_000_000
        counts_per_million.insert(0, "entrez_gene_ids", pd.Series(entrez_ids))
        logger.debug(f"Writing CPM matrix to {output_filepath}")
        counts_per_million.to_csv(output_filepath, index=False)

        # TODO: Counts per million is adding ~61,500 columns (equal to the number of genes) for some reason. Most likely due to multiplying by 1_000_000, not exactly sure why

        min_samples = round(n_exp * len(counts.columns))
        top_samples = round(n_top * len(counts.columns))
        test_bools = pd.DataFrame({"entrez_gene_ids": entrez_ids})
        for i in range(len(counts_per_million.columns)):
            cutoff = 10e6 / (np.median(np.sum(counts[:, i]))) if cut_off == "default" else 1e6 * cut_off / np.median(np.sum(counts[:, i]))
            test_bools = test_bools.merge(counts_per_million[counts_per_million.iloc[:, i] > cutoff])

    return metrics


def TPM_quantile_filter(*, metrics: NamedMetrics, filtering_options: _FilteringOptions) -> NamedMetrics:
    """
    Apply quantile-based filtering to the TPM matrix for a given sample.
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
        metric.entrez_gene_ids = [gene for gene, keep in zip(entrez_ids, min_genes) if keep]
        metric.gene_sizes = [gene for gene, keep in zip(gene_size, min_genes) if keep]
        metric.count_matrix = metric.count_matrix.iloc[min_genes, :]
        metric.normalization_matrix = metrics[sample].normalization_matrix.iloc[min_genes, :]

        keep_top_genes = [gene for gene, keep in zip(entrez_ids, top_genes) if keep]
        metric.high_confidence_entrez_gene_ids = [gene for gene, keep in zip(entrez_ids, keep_top_genes) if keep]

    metrics = calculate_z_score(metrics)

    return metrics


def zfpkm_filter(*, metrics: NamedMetrics, filtering_options: _FilteringOptions, prep: RNASeqPreparationMethod) -> NamedMetrics:
    n_exp = filtering_options.replicate_ratio
    n_top = filtering_options.high_replicate_ratio
    cut_off = filtering_options.cut_off

    metrics = calculate_fpkm(metrics)

    sample: str
    metric: _StudyMetrics
    for sample, metric in metrics.items():
        fpkm_df = metric.normalization_matrix
        fpkm_df = fpkm_df[fpkm_df.sum(axis=1) > 0]

        minimums = fpkm_df == 0
        results, zfpkm_df = zfpkm_transform(fpkm_df)
        zfpkm_df[minimums] = -4

        # TODO: Plot zFPKM results

        min_samples = round(n_exp * len(zfpkm_df.columns))
        top_samples = round(n_top * len(zfpkm_df.columns))

        min_func = k_over_a(min_samples, cut_off)
        top_func = k_over_a(top_samples, cut_off)

        min_genes: npt.NDArray[bool] = genefilter(zfpkm_df, min_func)
        top_genes: npt.NDArray[bool] = genefilter(zfpkm_df, top_func)

        # Only keep `entrez_gene_ids` that pass `min_genes`
        metric.entrez_gene_ids = [gene for gene, keep in zip(metric.entrez_gene_ids, min_genes) if keep]
        metric.high_confidence_entrez_gene_ids = [gene for gene, keep in zip(metric.entrez_gene_ids, top_genes) if keep]

    return metrics


def umi_filter(*, metrics: NamedMetrics, filtering_options: _FilteringOptions) -> NamedMetrics:
    for metric in metrics.values():
        entrez_ids = metric.entrez_gene_ids
        count_matrix = count_matrix[count_matrix.sum(axis=0) > 0, :]
        minimums = count_matrix[count_matrix == 0]

        zfpkm_results, z_matrix = zfpkm_transform(count_matrix)
        z_matrix[minimums] = -4
        z_matrix.insert(0, "entrez_gene_id", pd.Series(entrez_ids))
        zfpkm_plot(zfpkm_results, plot_xfloor=z_matrix.min().min())

        min_samples = round(filtering_options.replicate_ratio * len(z_matrix.columns))
        top_samples = round(filtering_options.high_replicate_ratio * len(z_matrix.columns))

        min_func = k_over_a(min_samples, filtering_options.cut_off)
        top_func = k_over_a(top_samples, filtering_options.cut_off)

        min_genes: npt.NDArray[bool] = genefilter(z_matrix, min_func)
        top_genes: npt.NDArray[bool] = genefilter(z_matrix, top_func)

        metric.entrez_gene_ids = [gene for gene, keep in zip(entrez_ids, min_genes) if keep]
        metric.high_confidence_entrez_gene_ids = [gene for gene, keep in zip(entrez_ids, top_genes) if keep]

    return metrics


def filter_counts(
    *,
    context_name: str,
    metrics: NamedMetrics,
    technique: FilteringTechnique,
    filtering_options: _FilteringOptions,
    prep: RNASeqPreparationMethod,
) -> NamedMetrics:
    match technique:
        case FilteringTechnique.cpm:
            return cpm_filter(context_name=context_name, metrics=metrics, filtering_options=filtering_options, prep=prep)
        case FilteringTechnique.tpm:
            return TPM_quantile_filter(metrics=metrics, filtering_options=filtering_options)
        case FilteringTechnique.zfpkm:
            return zfpkm_filter(metrics=metrics, filtering_options=filtering_options, prep=prep)
        case FilteringTechnique.umi:
            return umi_filter(metrics=metrics, filtering_options=filtering_options)
        case _:
            raise ValueError(f"Technique must be one of {FilteringTechnique}")


async def save_rnaseq_tests(
    context_name: str,
    counts_matrix_filepath: Path,
    config_filepath: Path,
    gene_info_filepath: Path,
    output_filepath: Path,
    prep: RNASeqPreparationMethod,
    replicate_ratio: float,
    batch_ratio: float,
    high_replicate_ratio: float,
    high_batch_ratio: float,
    technique: FilteringTechnique,
    cut_off: int | float,
):
    biodbnet: BioDBNet = BioDBNet()
    filtering_options = _FilteringOptions(
        replicate_ratio=replicate_ratio,
        batch_ratio=batch_ratio,
        cut_off=cut_off,
        high_replicate_ratio=high_replicate_ratio,
        high_batch_ratio=high_batch_ratio,
    )

    if prep == RNASeqPreparationMethod.SCRNA:
        technique = FilteringTechnique.umi
        logger.warning("Single cell filtration does not normalize and assumes gene counts are counted with Unique Molecular Identifiers (UMIs)")

    read_counts_results: ReadMatrixResults = await _read_counts_matrix(
        biodbnet=biodbnet,
        context_name=context_name,
        counts_matrix_filepath=counts_matrix_filepath,
        config_filepath=config_filepath,
        gene_info_filepath=gene_info_filepath,
    )
    metrics = read_counts_results.metrics
    entrez_gene_ids = read_counts_results.entrez_gene_ids

    metrics = filter_counts(
        context_name=context_name,
        metrics=metrics,
        technique=technique,
        filtering_options=filtering_options,
        prep=prep,
    )

    expressed_genes: list[str] = []
    top_genes: list[str] = []
    for metric in metrics.values():
        expressed_genes.extend(metric.entrez_gene_ids)
        top_genes.extend(metric.high_confidence_entrez_gene_ids)

    expression_frequency = pd.Series(expressed_genes).value_counts()
    expression_df = pd.DataFrame({"entrez_gene_id": expression_frequency.index, "frequency": expression_frequency.values})
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

    boolean_matrix.to_csv(output_filepath, index=False)


if __name__ == "__main__":
    import asyncio

    asyncio.run(
        save_rnaseq_tests(
            context_name="naiveCD4",
            counts_matrix_filepath=Path("/Users/joshl/Projects/COMO/main/data/data_matrices/naiveCD4/gene_counts_matrix_scrna_naiveCD4.h5ad"),
            config_filepath=Path("/Users/joshl/Projects/COMO/main/data/config_sheets/scrnaseq_data_inputs_auto.xlsx"),
            gene_info_filepath=Path("/Users/joshl/Projects/COMO/main/data/gene_info.csv"),
            output_filepath=Path("/Users/joshl/Projects/COMO/main/data/results/naiveCD4/scrna/rnaseq_scrna_naiveCD4.csv"),
            prep=RNASeqPreparationMethod.SCRNA,
            taxon_id=Taxon.HOMO_SAPIENS,
            replicate_ratio=0.5,
            high_replicate_ratio=1.0,
            batch_ratio=0.5,
            high_batch_ratio=1.0,
            technique=FilteringTechnique.umi,
            cut_off=-3,
        )
    )
