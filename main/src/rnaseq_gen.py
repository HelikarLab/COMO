#!/usr/bin/python3

import argparse
import os
import re
import sys
from pathlib import Path

import pandas as pd
import rpy2_api
from project import configs
from rpy2.robjects import pandas2ri

configs = project.Configs()

# enable r to py conversion
pandas2ri.activate()

r_file_path = Path(configs.root_dir, "src", "rscripts", "rnaseq.R")


def load_rnaseq_tests(filename, context_name, lib_type):
    """
    Load rnaseq results returning a dictionary of test (context, context, cell, etc ) names and rnaseq expression data
    """

    def load_dummy_dict():
        savepath = os.path.join(
            configs.root_dir,
            "data",
            "data_matrices",
            "placeholder",
            "placeholder_empty_data.csv",
        )
        dat = pd.read_csv(savepath, index_col="ENTREZ_GENE_ID")
        return "dummy", dat

    if (
        not filename or filename == "None"
    ):  # not using this data type, use empty dummy data matrix
        return load_dummy_dict()

    inquiry_full_path = os.path.join(
        configs.root_dir, "data", "config_sheets", filename
    )
    if not os.path.isfile(
        inquiry_full_path
    ):  # check that config file exist (isn't needed to load, but helps user)
        print(f"Error: Config file not found at {inquiry_full_path}")
        sys.exit()

    if lib_type == "total":  # if using total RNA-seq library prep
        filename = f"rnaseq_total_{context_name}.csv"
    elif lib_type == "mrna":  # if using mRNA-seq library prep
        filename = f"rnaseq_mrna_{context_name}.csv"
    elif lib_type == "scrna":  # if using single-cell RNA-seq
        filename = f"rnaseq_scrna_{context_name}.csv"
    else:
        print(
            f"Unsupported RNA-seq library type: {lib_type}. Must be one of 'total', 'mrna', 'sc'."
        )
        sys.exit()

    fullsavepath = os.path.join(
        configs.root_dir, "data", "results", context_name, lib_type, filename
    )

    if os.path.isfile(fullsavepath):
        data = pd.read_csv(fullsavepath, index_col="ENTREZ_GENE_ID")
        print(f"Read from {fullsavepath}")
        return context_name, data

    else:
        print(
            f"{lib_type} gene expression file for {context_name} was not found at {fullsavepath}. This may be "
            f"intentional. Contexts where {lib_type} data can be found in /work/data/results/{context_name}/ will "
            "still be used if found for other contexts."
        )
        return load_dummy_dict()


def handle_context_batch(
    config_filename,
    replicate_ratio,
    batch_ratio,
    replicate_ratio_high,
    batch_ratio_high,
    technique,
    quantile,
    min_count,
    min_zfpkm,
    prep,
):
    """
    Handle iteration through each context type and create rnaseq expression file by calling rnaseq.R
    """

    rnaseq_config_filepath = os.path.join(
        configs.root_dir, "data", "config_sheets", config_filename
    )
    xl = pd.ExcelFile(rnaseq_config_filepath)
    sheet_names = xl.sheet_names

    print(f"Reading config file: {rnaseq_config_filepath}")

    for context_name in sheet_names:
        print(f"\nStarting '{context_name}'")
        rnaseq_output_file = f"rnaseq_{prep}_{context_name}.csv"
        rnaseq_output_filepath = os.path.join(
            configs.data_dir, "results", context_name, prep, rnaseq_output_file
        )

        rnaseq_input_file = f"gene_counts_matrix_{prep}_{context_name}.csv"
        rnaseq_input_filepath = os.path.join(
            configs.data_dir, "data_matrices", context_name, rnaseq_input_file
        )

        if not os.path.exists(rnaseq_input_filepath):
            print(
                f"Gene counts matrix not found at {rnaseq_input_filepath}, skipping..."
            )
            continue

        gene_info_filepath = os.path.join(configs.data_dir, "gene_info.csv")
        os.makedirs(os.path.dirname(rnaseq_output_filepath), exist_ok=True)

        print(f"Gene info:\t\t{gene_info_filepath}")
        print(f"Count matrix:\t\t{rnaseq_input_filepath}")

        rpy2_api.Rpy2(
            r_file_path=r_file_path,
            counts_matrix_file=rnaseq_input_filepath,
            config_file=rnaseq_config_filepath,
            out_file=rnaseq_output_filepath,
            info_file=gene_info_filepath,
            context_name=context_name,
            prep=prep,
            replicate_ratio=replicate_ratio,
            batch_ratio=batch_ratio,
            replicate_ratio_high=replicate_ratio_high,
            batch_ratio_high=batch_ratio_high,
            technique=technique,
            quantile=quantile,
            min_count=min_count,
            min_zfpkm=min_zfpkm,
        ).call_function("save_rnaseq_tests")

        print(f"Results saved at:\t{rnaseq_output_filepath}")


def main(argv):
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
        dest="config_filename",
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
        dest="replicate_ratio_high",
        help="Ratio of replicates required for a gene to be considered high-confidence. "
        "High-confidence genes ignore consensus with other data-sources like proteomics or "
        "microarray. Example: 0.9 means that for a gene to be high-confidence, at least 90% of "
        "replicates in a group must pass the cutoff after normalization",
    )
    parser.add_argument(
        "-gh",
        "--high-batch-ratio",
        type=float,
        required=False,
        default=1.0,
        dest="batch_ratio_high",
        help="Ratio of groups (studies/batches) required for a gene to be considered high-confidence "
        "within that group. HHigh-confidence genes ignore consensus with other data-sources like "
        "proteomics or microarray. Example: 0.9 means that for a gene to be high-confidence, at "
        "least 90% of groups in a study must have passed the replicate ratio test",
    )
    parser.add_argument(
        "-t",
        "--filt-technique",
        type=str,
        required=False,
        default="quantile-tpm",
        dest="technique",
        help="Technique to normalize and filter counts with. Either 'zfpkm', 'quantile-tpm' or "
        "'flat-cpm'. More info about each method is discussed in pipeline.ipynb.",
    )
    parser.add_argument(
        "-q",
        "--quantile",
        type=int,
        required=False,
        default=25,
        dest="quantile",
        help="Cutoff used for quantile-tpm normalization and filtration technique. Example: 25 means "
        "that genes with TPM > 75% percentile wik=ll be considered active for that replicate.",
    )
    parser.add_argument(
        "-m",
        "--min-count",
        required=False,
        default="default",
        dest="min_count",
        help="Cutoff used for cpm. Minimum number of counts to be considered expressed, alternatively "
        "use 'default' to use method outlined in CITATION NEEDED",
    )
    parser.add_argument(
        "-z",
        "--min-zfpkm",
        required=False,
        default="default",
        dest="min_zfpkm",
        help="Cutoff used for zfpkm. Minimum zfpkm to be considered expressed, according to the paper, "
        "CITATION NEEDED, should be between -3 and -2.",
    )
    parser.add_argument(
        "-p",
        "--library-prep",
        required=False,
        default="",
        dest="prep",
        help="Library preparation used, will separate samples into groups to only compare similarly "
        "prepared libraries. For example, mRNA, total-rna, scRNA, etc",
    )

    args = parser.parse_args(argv)

    config_filename = args.config_filename
    replicate_ratio = args.replicate_ratio
    batch_ratio = args.batch_ratio
    replicate_ratio_high = args.replicate_ratio_high
    batch_ratio_high = args.batch_ratio_high
    technique = args.technique
    quantile = args.quantile
    min_count = args.min_count
    prep = args.prep
    min_zfpkm = args.min_zfpkm

    if re.search("tpm", technique.lower()) or re.search("quantile", technique.lower()):
        technique = "quantile"
    elif re.search("cpm", technique.lower()):
        technique = "cpm"
    elif re.search("zfpkm", technique.lower()):
        technique = "zfpkm"
    else:
        print(
            "Normalization-filtration technique not recognized. Must be 'tpm-quantile', 'cpm', or 'zfpkm'."
        )
        sys.exit()

    if int(quantile) > 100 or int(quantile) < 1:
        print("Quantile must be between 1 - 100")

    prep = prep.replace(" ", "")

    handle_context_batch(
        config_filename,
        replicate_ratio,
        batch_ratio,
        replicate_ratio_high,
        batch_ratio_high,
        technique,
        quantile,
        min_count,
        min_zfpkm,
        prep,
    )
    print("\nDone!")


if __name__ == "__main__":
    main(sys.argv[1:])
