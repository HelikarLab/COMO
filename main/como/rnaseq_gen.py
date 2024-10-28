import sys
from pathlib import Path

sys.path.insert(0, Path(__file__).parent.parent.as_posix())


import argparse
from dataclasses import dataclass
from typing import Optional

import pandas as pd
from loguru import logger

from como import Config, RNASeqPreparationMethod
from como.rnaseq import FilteringTechnique, save_rnaseq_tests


@dataclass
class RNASeqGen:
    config_file: str
    replicate_ratio: float
    batch_ratio: float
    high_replicate_ratio: float
    high_batch_ratio: float
    filtering_technique: FilteringTechnique
    minimum_cutoff: int | str
    library_prep: RNASeqPreparationMethod

    def __post_init__(self):
        self.library_prep = RNASeqPreparationMethod.from_string(str(self.library_prep))
        self.filtering_technique = FilteringTechnique.from_string(str(self.filtering_technique))

        if self.minimum_cutoff is None:
            if self.filtering_technique == FilteringTechnique.tpm:
                self.minimum_cutoff = 25
            elif self.filtering_technique == FilteringTechnique.cpm:
                self.minimum_cutoff = "default"
            elif self.filtering_technique == FilteringTechnique.zfpkm:
                self.minimum_cutoff = -3


def _handle_context_batch(
    config_filename,
    replicate_ratio,
    batch_ratio,
    replicate_ratio_high,
    batch_ratio_high,
    technique: FilteringTechnique,
    quantile,
    min_count,
    min_zfpkm,
    prep,
):
    """
    Handle iteration through each context type and create rnaseq expression file by calling rnaseq.R
    """
    r_file_path = Path(__file__).parent / "rscripts" / "rnaseq.R"
    if not r_file_path.exists():
        raise FileNotFoundError(f"Unable to find 'rnaseq.R'! Looking for: {r_file_path}")

    config = Config()

    config_filepath = config.config_dir / config_filename
    xl = pd.ExcelFile(config_filepath)
    sheet_names = xl.sheet_names

    logger.info(f"Reading config file: {config_filepath}")

    for context_name in sheet_names:
        logger.debug(f"Starting '{context_name}'")

        rnaseq_input_filepath = config.data_dir / "data_matrices" / context_name / f"gene_counts_matrix_{prep}_{context_name}.csv"
        if not rnaseq_input_filepath.exists():
            print(f"Gene counts matrix not found at {rnaseq_input_filepath}, skipping...")
            continue

        gene_info_filepath = config.data_dir / "gene_info.csv"
        rnaseq_output_filepath = config.result_dir / context_name / prep / f"rnaseq_{prep}_{context_name}.csv"
        rnaseq_output_filepath.parent.mkdir(parents=True, exist_ok=True)

        logger.info(f"Gene info saved at: {gene_info_filepath}")
        logger.info(f"Output file saved at: {rnaseq_output_filepath}")

        save_rnaseq_tests(
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
            quantile=quantile,
            min_count=min_count,
            min_zfpkm=min_zfpkm,
        )
        logger.success(f"Results saved at '{rnaseq_output_filepath}'")


def rnaseq_gen(
    config_filename: str,
    prep: RNASeqPreparationMethod,
    replicate_ratio: float = 0.5,
    batch_ratio: float = 0.5,
    replicate_ratio_high: float = 1.0,
    batch_ratio_high: float = 1.0,
    technique: FilteringTechnique | str = FilteringTechnique.tpm,
    cut_off: Optional[int] = None,
    prep: Optional[str] = "",
) -> None:
    if isinstance(technique, str):
        technique = FilteringTechnique(technique.lower())

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
    elif technique == FilteringTechnique.zfpkm:
        # if cut_off is not None and (cut_off < -3 or cut_off > -2):
        #     raise ValueError("Cutoff must be between -3 and -2")

        if cut_off is None:
            cut_off = "default"

    prep = prep.replace(" ", "")

    _handle_context_batch(
        config_filename,
        replicate_ratio,
        batch_ratio,
        replicate_ratio_high,
        batch_ratio_high,
        technique,
        cut_off,
        cut_off,
        cut_off,
        prep,
    )


def _parse_args():
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
        "-t",
        "--filt-technique",
        type=str,
        required=False,
        default="quantile-tpm",
        dest="filtering_technique",
        help="Technique to normalize and filter counts with. Either 'zfpkm', 'quantile-tpm' or "
        "'flat-cpm'. More info about each method is discussed in pipeline.ipynb.",
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
        required=False,
        default="",
        dest="library_prep",
        help="Library preparation used, will separate samples into groups to only compare similarly "
        "prepared libraries. For example, mRNA, total-rna, scRNA, etc",
    )
    args = parser.parse_args()
    args.filtering_technique = args.filtering_technique.lower()
    return RNASeqGen(**vars(args))


if __name__ == "__main__":
    """
    Generate a list of active and high-confidence genes from a counts matrix using a user defined
    at normalization-technique at /work/data/results/<context name>/rnaseq_<context_name>.csv
    Currently, can filter raw RNA-seq counts using three normalization techniques. Which are defined in rnaseq.R
    TPM Quantile, where each replicate is normalized with Transcripts-per-million and an upper quantile is taken to
    create a boolean list of active genes for the replicate. Replicates are compared for consensus within the
    study/batch number according to user-defined ratios and then study/batch numbers are checked for consensus
    according to different user defined ratios.   **CITATION NEEDED** **Recomended if user wants more control over the
    size of the model, like a smaller model that allows for only the most expressed reactions, or a larger more
    encompassing one that contains less essential reactions.
    zFPKM method outlined in: https://pubmed.ncbi.nlm.nih.gov/24215113/ can be used. Counts will be normalized using
    zFPKM and genes > -3 will be considered expressed per thier recommendation. Expressed genes will be checked for
    consensus at the replicate and study/batch levels the same as TPM Quantile. **Recommended if user wants to give
    least input over gene essentially determination and use the most standardized method of active gene determination.
    flat cutoff of CPM (counts per million) normalized values, check for consensus the same as other methods.
    """
    args = _parse_args()
    rnaseq_gen(
        config_filename=args.config_file,
        replicate_ratio=args.replicate_ratio,
        batch_ratio=args.batch_ratio,
        replicate_ratio_high=args.high_replicate_ratio,
        batch_ratio_high=args.high_batch_ratio,
        technique=args.filtering_technique,
        cut_off=args.minimum_cutoff,
        prep=args.library_prep,
    )
