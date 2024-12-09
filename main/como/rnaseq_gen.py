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


    def __post_init__(self):
        self.library_prep = RNAPrepMethod.from_string(str(self.library_prep))
        self.filtering_technique = FilteringTechnique.from_string(str(self.filtering_technique))

        if self.minimum_cutoff is None:
            if self.filtering_technique == FilteringTechnique.tpm:
                self.minimum_cutoff = 25
            elif self.filtering_technique == FilteringTechnique.cpm:
                self.minimum_cutoff = "default"
            elif self.filtering_technique == FilteringTechnique.zfpkm:
                self.minimum_cutoff = -3


async def _handle_context_batch(
    config_filename: str,
    replicate_ratio: float,
    batch_ratio: float,
    replicate_ratio_high: float,
    batch_ratio_high: float,
    technique: FilteringTechnique,
    cut_off: int | float | str,
    prep: RNAPrepMethod,
    taxon: Taxon,
) -> None:
    """Iterate through each context type and create rnaseq expression file.

    :param config_filename: The configuration filename to read
    :param replicate_ratio: The percentage of replicates that a gene must
        appear in for a gene to be marked as "active" in a batch/study
    :param batch_ratio: The percentage of batches that a gene must appear in for a gene to be marked as 'active"
    :param replicate_ratio_high: The percentage of replicates that a gene must
        appear in for a gene to be marked "highly confident" in its expression in a batch/study
    :param batch_ratio_high: The percentage of batches that a gene must
        appear in for a gene to be marked "highly confident" in its expression
    :param technique: The filtering technique to use
    :param cut_off: The cutoff value to use for the provided filtering technique
    :param prep: The library preparation method
    :param taxon: The NCBI Taxon ID
    :return: None
    """
    config = Config()

    config_filepath = config.config_dir / config_filename
    if not config_filepath.exists():
        raise FileNotFoundError(f"Unable to find '{config_filename}' at the path: '{config_filepath}'")
    xl = pd.ExcelFile(config_filepath)
    sheet_names = xl.sheet_names

    logger.info(f"Reading config file: {config_filepath}")

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
