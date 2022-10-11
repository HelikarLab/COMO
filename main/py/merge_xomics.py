#!/usr/bin/python3

import sys
import argparse
import os
import pandas as pd
import json
from collections import Counter
from warnings import warn
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from pathlib import Path

from create_context_specific_model import split_gene_expression_data
import microarray_gen
from project import configs
import proteomics_gen
import rnaseq_gen
import rpy2_api


# enable r to py conversion
pandas2ri.activate()

# import R libraries
ggplot2 = importr("ggplot2")
tidyverse = importr("tidyverse")

# read and translate R functions
# f = open(os.path.join(configs.rootdir, "py", "rscripts", "combine_distributions.R"), "r")
# string = f.read()
# f.close()
# combine_dist_io = SignatureTranslatedAnonymousPackage(string, "combine_dist_io")
r_file_path = Path(configs.rootdir, "py", "rscripts", "combine_distributions.R")


def merge_xomics(
    sheet,
    expression_requirement,
    microarray_file=None,
    proteomics_file=None,
    trnaseq_file=None,
    mrnaseq_file=None,
    scrnaseq_file=None,
    no_hc=False,
    no_na=False,
):
    """
    Merges microarray, rnaseq, and/or proteomics active gene logicals from outputs of their respective "_gen.py"
    scripts.

    :param microarray_file: filename of microarray config file in /main/data/config_sheets/
    :param proteomics_file: filename of proteomics config file in /main/data/config_sheets/
    :param trnaseq_file: filename of Total RNA-seq config file in /main/data/config_sheets/
    :param mrnaseq_file: filename of mRNA-seq config file in /main/data/config_sheets/
    :param scrnaseq_file: filename of single-cell RNA-seq config file in /main/data/config_sheets/
    :param no_hc: True if not adjusting for NA values (happens when gene missing from data source)
    :param no_na: filename of single-cell RNA-seq config file in /main/data/config_sheets/
    :param sheet: sheet name to use, should be context, context, cell type, etc
    :param expression_requirement: integer, minimum number of provided sources with active gene for a it to be in model
    :return: dictionary where keys are contexts, (tissue name, control type etc) and values are expression tables
    """
    print(f"Merging data for {sheet}")
    # load data for each source if it exists. IF not load an empty dummy dataset
    microarray = microarray_gen.load_microarray_tests(filename=microarray_file, context_name=sheet)
    proteomics = proteomics_gen.load_proteomics_tests(filename=proteomics_file, context_name=sheet)
    trnaseq = rnaseq_gen.load_rnaseq_tests(filename=trnaseq_file, context_name=sheet, lib_type="total")  # total RNA-seq
    mrnaseq = rnaseq_gen.load_rnaseq_tests(filename=mrnaseq_file, context_name=sheet, lib_type="mrna")  # mRNA-seq
    scrnaseq = rnaseq_gen.load_rnaseq_tests(filename=scrnaseq_file, context_name=sheet, lib_type="scrna")  # Single-cell RNA-seq

    files_dict = dict()

    exp_list = []
    high_list = []

    if proteomics[0] != "dummy":
        exp_list.append("prote_exp")
        high_list.append("prote_high")
        prote_data = proteomics[1].loc[:, ["expressed", "high"]]
        prote_data.rename(columns={"expressed": "prote_exp", "high": "prote_high"}, inplace=True)
        merge_data = prote_data

    if microarray[0] != "dummy":
        exp_list.append("trans_exp")
        high_list.append("trans_high")
        micro_data = microarray[1].loc[:, ["expressed", "high"]]
        micro_data.rename(columns={"expressed": "trans_exp", "high": "trans_high"}, inplace=True)
        if "merge_data" not in locals():
            merge_data = micro_data
        else:
            merge_data = merge_data.join(micro_data, how="outer")

    if trnaseq[0] != "dummy":
        exp_list.append("trnaseq_exp")
        high_list.append("trnaseq_high")
        trnaseq_data = trnaseq[1].loc[:, ["expressed", "high"]]
        trnaseq_data.rename(columns={"expressed": "trnaseq_exp", "high": "trnaseq_high"}, inplace=True)
        if "merge_data" not in locals():
            merge_data = trnaseq_data
        else:
            merge_data = merge_data.join(trnaseq_data, how="outer")

    if mrnaseq[0] != "dummy":
        exp_list.append("mrnaseq_exp")
        high_list.append("mrnaseq_high")
        mrnaseq_data = mrnaseq[1].loc[:, ["expressed", "high"]]
        mrnaseq_data.rename(columns={"expressed": "mrnaseq_exp", "high": "mrnaseq_high"}, inplace=True)
        if "merge_data" not in locals():
            merge_data = mrnaseq_data
        else:
            merge_data = merge_data.join(mrnaseq_data, how="outer")

    if scrnaseq[0] != "dummy":
        exp_list.append("scrnaseq_exp")
        high_list.append("scrnaseq_high")
        scrnaseq_data = scrnaseq[1].loc[:, ["expressed", "high"]]
        scrnaseq_data.rename(columns={"expressed": "scrnaseq_exp", "high": "scrnaseq_high"}, inplace=True)
        if "merge_data" not in locals():
            merge_data = scrnaseq_data
        else:
            merge_data = merge_data.join(scrnaseq_data, how="outer")

    merge_data = microarray_gen.mergeLogicalTable(merge_data)

    num_sources = len(exp_list)
    merge_data["Active"] = 0
    merge_data["Required"] = 0

    if no_na:  # dont adjust for na values
        merge_data.loc[:, "Required"] = merge_data[exp_list].apply(
            lambda x: expression_requirement
            if (expression_requirement - (num_sources - x.count()) > 0)
            else 1,
            axis=1,
        )
    else:  # subtract one from requirement per NA
        merge_data.loc[:, "Required"] = merge_data[exp_list].apply(
            lambda x: expression_requirement - (num_sources - x.count())
            if (expression_requirement - (num_sources - x.count()) > 0)
            else 1,
            axis=1,
        )

    # count number of sources gene is active in. Set to active in final output if at least adjusted expression reqirmnt
    merge_data["TotalExpressed"] = merge_data[exp_list].sum(axis=1)
    merge_data.loc[merge_data[exp_list].sum(axis=1) >= merge_data["Required"], "Active"] = 1

    if not no_hc:  # set genes that are high-confidence in at least one data source to active
        merge_data.loc[merge_data[high_list].sum(axis=1) > 0, "Active"] = 1

    #merge_data = merge_data.astype(int)
    merge_data = merge_data

    filepath = os.path.join(configs.rootdir, "data", "results", sheet, f"merged_{sheet}.csv")
    merge_data.to_csv(filepath, index_label="ENTREZ_GENE_ID")

    filepath = os.path.join(configs.rootdir, "data", "results", sheet, f"ActiveGenes_{sheet}_Merged.csv")
    merge_data.reset_index(drop=False, inplace=True)

    split_entrez = split_gene_expression_data(merge_data)
    split_entrez.rename(columns={"Gene": "ENTREZ_GENE_ID", "Data": "Active"}, inplace=True)
    split_entrez.to_csv(filepath, index_label="ENTREZ_GENE_ID")
    files_dict[sheet] = filepath

    print(f"{sheet}: save to {filepath}\n")

    return files_dict


def handle_context_batch(
    microarray_file,
    trnaseq_file,
    mrnaseq_file,
    scrnaseq_file,
    proteomics_file,
    tweight,
    mweight,
    sweight,
    pweight,
    expression_requirement,
    adjust_method,
    no_hc,
    no_na,
    custom_df,
    merge_distro,
    keep_gene_score
):
    """
    Handle merging of different data sources for each context type
    """
    sheet_names = []
    for file in [microarray_file, trnaseq_file, mrnaseq_file, scrnaseq_file, proteomics_file]:
        if file is not None:
            config_filepath = os.path.join(configs.rootdir, "data", "config_sheets", file)
            xl = pd.ExcelFile(config_filepath)
            sheet_names += xl.sheet_names

    use_trna = True if trnaseq_file is not None else False
    use_mrna = True if mrnaseq_file is not None else False
    use_scrna = True if scrnaseq_file is not None else False
    use_proteins = True if proteomics_file is not None else False

    counts = Counter(sheet_names)
    sheet_names = sorted(list(set(sheet_names)))
    print(f"Will merge data for: {sheet_names}")
    dict_list = {}

    max_inputs = max(counts.values())
    min_inputs = min(counts.values())

    if merge_distro:
        rpy2_api.Rpy2(
            r_file_path,
            os.path.join(configs.datadir, "results"),
            sheet_names,
            use_mrna,
            use_trna,
            use_scrna,
            use_proteins,
            keep_gene_score,
            tweight,
            mweight,
            sweight,
            pweight
        ).call_function("combine_zscores_main")
        # combine_dist_io.call_function("combine_zscores_main")
        # combine_dist_io.combine_zscores_main(
        #     os.path.join(configs.datadir, "results"),
        #     sheet_names,
        #     use_mrna,
        #     use_trna,
        #     use_scrna,
        #     use_proteins,
        #     keep_gene_score,
        #     tweight,
        #     mweight,
        #     sweight,
        #     pweight
        # )

    files_json = os.path.join(configs.rootdir, "data", "results", "step1_results_files.json")
    for context_name in sheet_names:
        num_sources = counts[context_name]
        if adjust_method == "progressive":
            exp_req = (num_sources - min_inputs) + expression_requirement
        elif adjust_method == "regressive":
            exp_req = expression_requirement - (max_inputs - num_sources)
        elif adjust_method == "flat":
            exp_req = expression_requirement
        else:
            exp_req = int(custom_df.iloc[custom_df["context"] == context_name, "req"].iloc[0])

        print(
            f"Expression requirement of {expression_requirement} adjusted to {exp_req} using {adjust_method} "
            f"adjustment method for {context_name}."
        )

        if exp_req > num_sources:
            warn(
                f"Expression requirement for {context_name} was calculated to be greater than max number of input "
                f"data sources. Will be force changed to {num_sources} to prevent output from having 0 active genes. "
                f"Consider lowering the expression requirement or changing the adjustment method."
            )
            exp_req = num_sources

        if exp_req < 1:  # never allow expression requirement to be less than one
            warn(
                f"Expression requirement for {context_name} was calculated to be lower than 1. Will be force changed "
                f"to 1 to prevent every gene being active."
            )
            exp_req = 1

        files_dict = merge_xomics(
            context_name,
            expression_requirement=exp_req,
            microarray_file=microarray_file,
            proteomics_file=proteomics_file,
            trnaseq_file=trnaseq_file,
            mrnaseq_file=mrnaseq_file,
            scrnaseq_file=scrnaseq_file,
            no_hc=no_hc,
            no_na=no_na,
        )

        dict_list.update(files_dict)

    with open(files_json, "w") as fp:
        json.dump(dict_list, fp)

    return


def main(argv):
    """
    Merge expression tables of multiple sources, microarray, RNA-seq, and/or proteomics into one list
    User can specify the number of sources with an active gene in order for it to be considered active in the model.
    Otherwise, it defaults to the number of sources provided. High-confidence genes from any source will be considered
    active in the model, regardless of agreement with other sources.
    """
    parser = argparse.ArgumentParser(
        prog="merge_xomics.py",
        description="Merge expression tables of multiple sources, microarray, RNA-seq, and/or proteomics into one",
        epilog="For additional help, please post questions/issues in the MADRID GitHub repo at "
        "https://github.com/HelikarLab/MADRID or email babessell@gmail.com",
    )

    parser.add_argument(
        "-d",
        "--merge-distribution",
        action="store_true",
        required=False,
        default=False,
        dest="merge_distro",
        help="Flag to merge zFPKM distributions. Required if using iMAT reconstruction algorithm in "
        "create_context_specific_model.py. Must have run rnaseq_gen.py with 'zFPKM' as "
        "'--technique'. If --proteomics-config-file is given will merge proteomics distributions "
        "with zFPKM distributions using a weighted scheme.",
    )

    parser.add_argument(
        "-k",
        "--keep-gene-scores",
        action="store_true",
        required=False,
        default=True,
        dest="keep_gene_score",
        help="When merging z-score distributions of expression, if using both protein abundance and transcipt zFPKM "
             "flag true if you wish to keep z-score of genes with no protein data, flag false if you wish to discard "
             "and treat as no expression"
    )

    parser.add_argument(
        "-a",
        "--microarray-config-file",
        type=str,
        required=False,
        default=None,
        dest="microarray_file",
        help="Name of microarray config .xlsx file in the /main/data/config_files/.",
    )

    parser.add_argument(
        "-t",
        "--total-rnaseq-config-file",
        type=str,
        required=False,
        default=None,
        dest="trnaseq_file",
        help="Name of total RNA-seq config .xlsx file in the /main/data/config_files/.",
    )

    parser.add_argument(
        "-m",
        "--mrnaseq-config-file",
        type=str,
        required=False,
        default=None,
        dest="mrnaseq_file",
        help="Name of mRNA-seq config .xlsx file in the /main/data/config_files/.",
    )

    parser.add_argument(
        "-s",
        "--scrnaseq-config-file",
        type=str,
        required=False,
        default=None,
        dest="scrnaseq_file",
        help="Name of RNA-seq config .xlsx file in the /main/data/config_files/.",
    )

    parser.add_argument(
        "-p",
        "--proteomics-config-file",
        type=str,
        required=False,
        default=None,
        dest="proteomics_file",
        help="Name of proteomics config .xlsx file in the /main/data/config_files/.",
    )

    parser.add_argument(
        "-e",
        "--expression-requirement",
        type=str,
        required=False,
        default=None,
        dest="expression_requirement",
        help="Number of sources with active gene for it to be considered active even if it is not a "
        "high confidence-gene",
    )

    parser.add_argument(
        "-r",
        "--requirement-adjust",
        type=str,
        required=False,
        default="flat",
        dest="adjust_method",
        help="Technique to adjust expression requirement based on differences in number of provided "
        "data source types.",
    )

    parser.add_argument(
        "-c",
        "--custom-requirement-file",
        required="custom" in argv,  # required if --requriement-adjust is "custom",
        dest="custom_file",
        default="SKIP",
        help="Name of .xlsx file where first column is context names and second column is expression "
        "requirement for that context, in /main/data/",
    )

    parser.add_argument(
        "-hc",
        "--no-hc",
        action="store_true",
        required=False,
        default=False,
        dest="no_hc",
        help="Flag to prevent high-confidence genes forcing a gene to be used in final model "
        "irrespective of other other data sources",
    )

    parser.add_argument(
        "-na",
        "--no-na-adjustment",
        action="store_true",
        required=False,
        default=False,
        dest="no_na",
        help="Flag to prevent genes missing in a data source library, but present in others from "
        "subtracting 1 from the expression requirement per data source that gene is missing in",
    )

    parser.add_argument(
        "-tw",
        "--total-rnaseq-weight",
        required=False,
        default=1,
        type=float,
        dest="tweight",
        help="Total RNA-seq weight for merging zFPKM distribution",
    )

    parser.add_argument(
        "-mw",
        "--mrnaseq-weight",
        required=False,
        default=1,
        type=float,
        dest="mweight",
        help="PolyA enriched (messenger) RNA-seq weight for merging zFPKM distribution",
    )

    parser.add_argument(
        "-sw",
        "--single-cell-rnaseq-weight",
        required=False,
        default=1,
        type=float,
        dest="sweight",
        help="Single-cell RNA-seq weight for merging zFPKM distribution",
    )

    parser.add_argument(
        "-pw",
        "--protein-weight",
        required=False,
        default=2,
        type=float,
        dest="pweight",
        help="Proteomics weight for merging z-score distribution",
    )


    args = parser.parse_args(argv)

    microarray_file = args.microarray_file
    proteomics_file = args.proteomics_file
    trnaseq_file = args.trnaseq_file
    mrnaseq_file = args.mrnaseq_file
    scrnaseq_file = args.scrnaseq_file
    expression_requirement = args.expression_requirement
    adjust_method = args.adjust_method.lower()
    custom_file = args.custom_file
    no_hc = args.no_hc
    no_na = args.no_na
    merge_distro = args.merge_distro
    keep_gene_score = args.keep_gene_score
    tweight = args.tweight
    mweight = args.mweight
    sweight = args.sweight
    pweight = args.pweight

    # read custom expression requirment file if used
    if custom_file != "SKIP":
        custom_filepath = os.path.join(configs.datadir, custom_file)
        custom_df = pd.read_excel(custom_filepath, sheet_name=0)
        custom_df.columns = ["context", "req"]
    else:
        custom_df = pd.DataFrame([])

    def_exp_req = sum(
        [1 for test in [microarray_file,
                        trnaseq_file,
                        mrnaseq_file,
                        scrnaseq_file,
                        proteomics_file]
            if test is None]
    )

    if expression_requirement.lower() == "default":
        expression_requirement = def_exp_req

    else:
        try:
            expression_requirement = int(expression_requirement)
            if int(expression_requirement) < 1:
                print("Expression requirement must be at least 1!")
                sys.exit(1)
        except ValueError:
            print("Expression requirement must be able to be converted to an integer!")
            sys.exit(1)

    if adjust_method not in ["progressive", "regressive", "flat", "custom"]:
        print(
            "Adjust method must be either 'progressive', 'regressive', 'flat', or 'custom'"
        )
        sys.exit(1)

    handle_context_batch(
        microarray_file,
        trnaseq_file,
        mrnaseq_file,
        scrnaseq_file,
        proteomics_file,
        tweight,
        mweight,
        sweight,
        pweight,
        expression_requirement,
        adjust_method,
        no_hc,
        no_na,
        custom_df,
        merge_distro,
        keep_gene_score
    )

    return


if __name__ == "__main__":
    main(sys.argv[1:])
