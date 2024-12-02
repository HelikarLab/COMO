from __future__ import annotations

import argparse
import asyncio
from dataclasses import dataclass
from typing import Optional

import pandas as pd
from fast_bioservices import Taxon
from loguru import logger

from como import Config
from como.custom_types import RNASeqPreparationMethod
from como.rnaseq import FilteringTechnique, save_rnaseq_tests


@dataclass
class _Arguments:
    config_file: str
    replicate_ratio: float
    batch_ratio: float
    high_replicate_ratio: float
    high_batch_ratio: float
    filtering_technique: FilteringTechnique
    minimum_cutoff: int | str
    library_prep: RNASeqPreparationMethod
    taxon: Taxon

    def __post_init__(self):
        self.library_prep = RNASeqPreparationMethod.from_string(str(self.library_prep))
        self.filtering_technique = FilteringTechnique.from_string(str(self.filtering_technique))

        if not str(self.taxon).isdigit():
            raise ValueError(f"Expected '--taxon-id' to be an integer; got: {self.taxon}")
        self.taxon = Taxon.from_int(int(self.taxon))

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
    prep: RNASeqPreparationMethod,
    taxon: Taxon,
) -> None:
    """
    Handle iteration through each context type and create rnaseq expression file

    :param config_filename: The configuration filename to read
    :param replicate_ratio: The percentage of replicates that a gene must appear in for a gene to be marked as "active" in a batch/study
    :param batch_ratio: The percentage of batches that a gene must appear in for a gene to be marked as 'active"
    :param replicate_ratio_high: The percentage of replicates that a gene must appear in for a gene to be marked "highly confident" in its expression in a batch/study
    :param batch_ratio_high: The percentage of batches that a gene must appear in for a gene to be marked "highly confident" in its expression
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

        rnaseq_input_filepath = config.data_dir / "data_matrices" / context_name / f"gene_counts_matrix_{prep.value}_{context_name}"
        if prep == RNASeqPreparationMethod.SCRNA:
            rnaseq_input_filepath = rnaseq_input_filepath.with_suffix(".h5ad")
        elif prep in {RNASeqPreparationMethod.TOTAL, RNASeqPreparationMethod.MRNA}:
            rnaseq_input_filepath = rnaseq_input_filepath.with_suffix(".csv")

        if not rnaseq_input_filepath.exists():
            logger.warning(f"Gene counts matrix not found at {rnaseq_input_filepath}, skipping...")
            continue

        gene_info_filepath = config.data_dir / "gene_info.csv"
        rnaseq_output_filepath = config.result_dir / context_name / prep.value / f"rnaseq_{prep.value}_{context_name}.csv"
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
    config_filename: str,
    prep: RNASeqPreparationMethod,
    taxon_id: int | str | Taxon,
    replicate_ratio: float = 0.5,
    high_replicate_ratio: float = 1.0,
    batch_ratio: float = 0.5,
    high_batch_ratio: float = 1.0,
    technique: FilteringTechnique | str = FilteringTechnique.tpm,
    cut_off: Optional[int | float] = None,
) -> None:
    """
    The main entrypoint for processing bulk total, bulk mRNA, and single-cell RNA sequencing data.
    Generate a list of active and high-confidence genes from a gene count matrix.
    Replicates are compared for consensus within the study/batch number according to replicate ratios, then study/batch numbers are checked for consensus according to batch ratios.
    The zFPKM method is outlined here: https://pubmed.ncbi.nlm.nih.gov/24215113/

    :param config_filename: The configuration filename to read
    :param prep: The preparation method
    :param taxon_id: The NCBI Taxon ID
    :param replicate_ratio: The percentage of replicates that a gene must appear in for a gene to be marked as "active" in a batch/study
    :param batch_ratio: The percentage of batches that a gene must appear in for a gene to be marked as 'active"
    :param high_replicate_ratio: The percentage of replicates that a gene must appear in for a gene to be marked "highly confident" in its expression in a batch/study
    :param high_batch_ratio: The percentage of batches that a gene must appear in for a gene to be marked "highly confident" in its expression
    :param technique: The filtering technique to use
    :param cut_off: The cutoff value to use for the provided filtering technique
    :return: None
    """
    if isinstance(technique, str):
        technique = FilteringTechnique(technique.lower())
    if isinstance(taxon_id, (str, int)):
        taxon_id = Taxon.from_string(str(taxon_id))

    if technique.value not in [t.value for t in FilteringTechnique]:
        raise ValueError(f"Technique must be one of {FilteringTechnique}")

    if technique == FilteringTechnique.tpm:
        if cut_off is None:
            cut_off = 25

        if cut_off < 1 or cut_off > 100:
            raise ValueError("Quantile must be between 1 - 100")

    elif technique == FilteringTechnique.cpm:
        if cut_off is not None and cut_off < 0:
            raise ValueError("Cutoff must be greater than 0")

        if cut_off is None:
            cut_off = "default"
    elif technique == FilteringTechnique.zfpkm and cut_off is None:
        cut_off = "default"

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


def _parse_args() -> _Arguments:
    """
    Parse the command line arguments
    :return: The parsed arguments
    """
    parser = argparse.ArgumentParser(
        prog="rnaseq_gen.py",
        description="Generate a list of active and high-confidence genes from a counts matrix using a user defined "
        "at normalization-technique at /work/data/results/<context name>/rnaseq_<context_name>.csv: "
        "https://github.com/HelikarLab/FastqToGeneCounts",
        epilog="For additional help, please post questions/issues in the MADRID GitHub repo at "
        "https://github.com/HelikarLab/MADRID or email babessell@gmail.com",
    )
    parser.add_argument(
        "-c",
        "--config-file",
        type=str,
        required=True,
        dest="config_file",
        help="Name of config .xlsx file in the /work/data/config_files/. Can be generated using "
        "rnaseq_preprocess.py or manually created and imported into the Juypterlab",
    )
    parser.add_argument(
        "-r",
        "--replicate-ratio",
        type=float,
        required=False,
        default=0.5,
        dest="replicate_ratio",
        help="Ratio of replicates required for a gene to be active within that study/batch group "
        "Example: 0.7 means that for a gene to be active, at least 70% of replicates in a group "
        "must pass the cutoff after normalization",
    )
    parser.add_argument(
        "-g",
        "--batch-ratio",
        type=float,
        required=False,
        default=0.5,
        dest="batch_ratio",
        help="Ratio of groups (studies or batches) required for a gene to be active "
        "Example: 0.7 means that for a gene to be active, at least 70% of groups in a study must  "
        "have passed the replicate ratio test",
    )
    parser.add_argument(
        "-rh",
        "--high-replicate-ratio",
        type=float,
        required=False,
        default=1.0,
        dest="high_replicate_ratio",
        help="Ratio of replicates required for a gene to be considered high-confidence. "
        "High-confidence genes ignore consensus with other data-sources, such as proteomics. "
        "Example: 0.9 means that for a gene to be high-confidence, at least 90% of replicates in a group must pass the cutoff after normalization",
    )
    parser.add_argument(
        "-gh",
        "--high-batch-ratio",
        type=float,
        required=False,
        default=1.0,
        dest="high_batch_ratio",
        help="Ratio of groups (studies/batches) required for a gene to be considered high-confidence within that group. "
        "High-confidence genes ignore consensus with other data-sources, like proteomics. "
        "Example: 0.9 means that for a gene to be high-confidence, at least 90% of groups in a study must have passed the replicate ratio test",
    )
    parser.add_argument(
        "--taxon",
        "--taxon-id",
        type=str,
        required=True,
        dest="taxon",
        help="The NCBI Taxonomy ID that is being proessed. '9606' for humans, '10090' for mice.",
    )
    parser.add_argument(
        "-t",
        "--filt-technique",
        type=str,
        required=False,
        default="quantile",
        dest="filtering_technique",
        help="Technique to normalize and filter counts with. Either 'zfpkm', 'quantile', or 'cpm'. More info about each method is discussed in pipeline.ipynb.",
    )
    parser.add_argument(
        "--minimum-cutoff",
        type=int,
        required=False,
        default=None,
        dest="minimum_cutoff",
        help="The minimum cutoff used for the filtration technique. If the filtering technique is zFPKM, the default is -3. If the filtering technique is quantile-tpm, the default is 25. If the filtering technique is flat-cpm, the default is determined dynamically. If the filtering technique is quantile, the default is 25.",
    )
    parser.add_argument(
        "-p",
        "--library-prep",
        required=True,
        choices=["total", "mrna", "scrna"],
        dest="library_prep",
        help="Library preparation used, will separate samples into groups to only compare similarly prepared libraries. For example, mRNA, total-rna, scRNA, etc",
    )
    args = parser.parse_args()
    args.filtering_technique = args.filtering_technique.lower()
    return _Arguments(**vars(args))


if __name__ == "__main__":
    args = _parse_args()
    asyncio.run(
        rnaseq_gen(
            config_filename=args.config_file,
            replicate_ratio=args.replicate_ratio,
            batch_ratio=args.batch_ratio,
            high_replicate_ratio=args.high_replicate_ratio,
            high_batch_ratio=args.high_batch_ratio,
            technique=args.filtering_technique,
            cut_off=args.minimum_cutoff,
            prep=args.library_prep,
            taxon_id=args.taxon,
        )
    )
