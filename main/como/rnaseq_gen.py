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



    for context_name in sheet_names:
        logger.debug(f"Starting '{context_name}'")

        rnaseq_input_filepath = (
            config.data_dir / "data_matrices" / context_name / f"gene_counts_matrix_{prep.value}_{context_name}"
        )
        if prep == RNAPrepMethod.SCRNA:
            rnaseq_input_filepath = rnaseq_input_filepath.with_suffix(".h5ad")
        elif prep in {RNAPrepMethod.TOTAL, RNAPrepMethod.MRNA}:
            rnaseq_input_filepath = rnaseq_input_filepath.with_suffix(".csv")

        if not rnaseq_input_filepath.exists():
            logger.warning(f"Gene counts matrix not found at {rnaseq_input_filepath}, skipping...")
            continue

        gene_info_filepath = config.data_dir / "gene_info.csv"
        rnaseq_output_filepath = (
            config.result_dir / context_name / prep.value / f"rnaseq_{prep.value}_{context_name}.csv"
        )
        rnaseq_output_filepath.parent.mkdir(parents=True, exist_ok=True)

        await save_rnaseq_tests(
            context_name=context_name,
            counts_matrix_filepath=rnaseq_input_filepath,
            config_filepath=config_filepath,
            output_filepath=rnaseq_output_filepath.as_posix(),
            gene_info_filepath=gene_info_filepath,
            prep=prep,
            replicate_ratio=replicate_ratio,
            batch_ratio=batch_ratio,
            high_replicate_ratio=replicate_ratio_high,
            high_batch_ratio=batch_ratio_high,
            technique=technique,
            cut_off=cut_off,
            taxon_id=taxon,
        )
        logger.success(f"Results saved at '{rnaseq_output_filepath}'")


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
