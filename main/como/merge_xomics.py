#!/usr/bin/python3

import argparse
import json
import re
import sys
from collections import Counter
from enum import Enum

import numpy as np
import pandas as pd
from fast_bioservices import BioDBNet, Input, Output

from como import proteomics_gen, rnaseq_gen, rpy2_api
from como.como_utilities import split_gene_expression_data
from como.project import Config

config = Config()
r_file_path = config.code_dir / "rscripts" / "combine_distributions.R"


class _MergedHeaderNames:
    PROTEOMICS = "prote"
    TRNASEQ = "trnaseq"
    MRNASEQ = "mrnaseq"
    SCRNASEQ = "scrnaseq"


class _ExpressedHeaderNames:
    PROTEOMICS = f"{_MergedHeaderNames.PROTEOMICS}_exp"
    TRNASEQ = f"{_MergedHeaderNames.TRNASEQ}_exp"
    MRNASEQ = f"{_MergedHeaderNames.MRNASEQ}_exp"
    SCRNASEQ = f"{_MergedHeaderNames.SCRNASEQ}_exp"


class _HighExpressionHeaderNames:
    PROTEOMICS = f"{_MergedHeaderNames.PROTEOMICS}_high"
    TRNASEQ = f"{_MergedHeaderNames.TRNASEQ}_high"
    MRNASEQ = f"{_MergedHeaderNames.MRNASEQ}_high"
    SCRNASEQ = f"{_MergedHeaderNames.SCRNASEQ}_high"


# Merge Output
def merge_logical_table(df: pd.DataFrame):
    """
    Merge the Rows of Logical Table belongs to the same ENTREZ_GENE_ID
    :param df:
    :return: pandas dataframe of merged table
    """
    # step 1: get all plural ENTREZ_GENE_IDs in the input table, extract unique IDs

    df.reset_index(drop=False, inplace=True)
    df.dropna(axis=0, subset=["ENTREZ_GENE_ID"], inplace=True)
    df["ENTREZ_GENE_ID"] = df["ENTREZ_GENE_ID"].str.replace(" /// ", "//").astype(str)

    single_entrez_ids: list[str] = df[~df["ENTREZ_GENE_ID"].str.contains("//")]["ENTREZ_GENE_ID"].tolist()
    multiple_entrez_ids: list[str] = df[df["ENTREZ_GENE_ID"].str.contains("//")]["ENTREZ_GENE_ID"].tolist()
    id_list: list[str] = []
    for i in multiple_entrez_ids:
        ids = i.split("//")
        id_list.extend(ids)

        duplicate_rows = pd.DataFrame([])
        for j in ids:
            rows = df.loc[df["ENTREZ_GENE_ID"] == i].copy()
            rows["ENTREZ_GENE_ID"] = j
            duplicate_rows = pd.concat([duplicate_rows, rows], axis=0)

        df = pd.concat([df, pd.DataFrame(duplicate_rows)], axis=0, ignore_index=True)
        df.drop(df[df["ENTREZ_GENE_ID"] == i].index, inplace=True)

    full_entrez_id_sets: set[str] = set()
    entrez_dups_list: list[list[str]] = []
    multi_entrez_index = list(range(len(multiple_entrez_ids)))

    for i in range(len(multiple_entrez_ids)):
        if i not in multi_entrez_index:
            continue

        set1 = set(multiple_entrez_ids[i].split("//"))
        multi_entrez_index.remove(i)

        for j in multi_entrez_index:
            set2 = set(multiple_entrez_ids[j].split("//"))
            intersect = set1.intersection(set2)
            if bool(intersect):
                set1 = set1.union(set2)
                multi_entrez_index.remove(j)

        sortlist = list(set1)
        sortlist.sort(key=int)
        new_entrez_id = " /// ".join(sortlist)
        full_entrez_id_sets.add(new_entrez_id)

    for i in full_entrez_id_sets:
        entrez_dups_list.append(i.split(" /// "))
    entrez_dups_dict = dict(zip(full_entrez_id_sets, entrez_dups_list))

    for merged_entrez_id, entrez_dups_list in entrez_dups_dict.items():
        df["ENTREZ_GENE_ID"].replace(to_replace=entrez_dups_list, value=merged_entrez_id, inplace=True)

    df.set_index("ENTREZ_GENE_ID", inplace=True)
    df = df.fillna(-1).groupby(level=0).max()
    df.replace(-1, np.nan, inplace=True)

    # TODO: Test if this is working properly
    """
    There seems to be an error when running Step 2.1 in the pipeline.ipynb file
    The commented-out return statement tries to return the df_output dataframe values as integers, but NaN values exist
        Because of this, it is unable to do so.
    If we change this to simply output the database, the line "np.where(posratio >= top_proportion . . ." (line ~162)
        Fails because it is comparing floats and strings

    I am unsure what to do in this situation
    """
    # return df_output.astype(int)
    return df


def get_transcriptmoic_details(merged_df: pd.DataFrame) -> pd.DataFrame:
    """
    This function will get the following details of transcriptomic data:
    - Gene Symbol
    - Gene Name
    - ENTREZ_GENE_ID

    The resulting dataframe will have its columns created in the order listed above
    It will return a pandas dataframe with this information

    :param merged_df: A dataframe containing all the transcriptomic and proteomic data after determining which genes are active
    :return: A dataframe with the above-listed columns
    """
    # import numpy as np
    # merged_df["prote_exp"] = np.random.choice([0, 1], size=len(merged_df))
    # merged_df["prote_high"] = pd.NA
    # for row in merged_df.itertuples():
    #     if row.prote_exp == 1:
    #         merged_df.at[row.Index, "TotalExpressed"] += 1

    # merged_df["TotalExpressed"] += 1 if merged_df["prote_exp"] == 1 else 0

    # If _ExpressedHeaderNames.PROTEOMICS.value is in the dataframe, lower the required expression by 1
    # We are only trying to get details for transcriptomic data
    if _ExpressedHeaderNames.PROTEOMICS in merged_df.columns:
        # Get the number of sources required for a gene to be marked "expressed"
        required_expression = merged_df["Required"].iloc[0]

        # Subtract 1 from merged_df["TotalExpressed"] if the current value is greater than or equal to 1
        # This is done to take into account the removal of proteomic expression
        merged_df["TotalExpressed"] = merged_df["TotalExpressed"].apply(lambda x: x - 1 if x >= 1 else x)

        # Subtract required_expression by 1 if it is greater than 1
        if required_expression > 1:
            required_expression -= 1

        # Create a new dataframe without [_ExpressedHeaderNames.PROTEOMICS.value, _HighExpressionHeaderNames.PROTEOMICS.value] columns
        transcriptomic_df: pd.DataFrame = merged_df.drop(
            columns=[
                _ExpressedHeaderNames.PROTEOMICS,
                _HighExpressionHeaderNames.PROTEOMICS,
            ],
            inplace=False,
        )

        # Must recalculate TotalExpressed because proteomic data was removed
        # If the TotalExpressed column is less than the Required column, set active to 1, otherwise set it to 0
        transcriptomic_df.loc[
            transcriptomic_df["TotalExpressed"] >= transcriptomic_df["Required"],
            "Active",
        ] = 1

    else:
        transcriptomic_df: pd.DataFrame = merged_df.copy()

    biodbnet = BioDBNet()
    gene_details: pd.DataFrame = biodbnet.db2db(
        input_values=transcriptomic_df.index.astype(str).values.tolist(),
        input_db=Input.GENE_ID,
        output_db=[
            Output.GENE_SYMBOL,
            Output.ENSEMBL_GENE_INFO,
            Output.GENE_INFO,
        ],
    )
    gene_details["entrez_gene_id"] = gene_details.index
    gene_details.reset_index(drop=True, inplace=True)

    # Apply regex to search for "[Description: XXXXXX]" and retrieve the XXXXXX
    # It excludes the square brackets and "Description: ", and only returns the description
    # descriptions: list[str] = [
    gene_details["description"] = [
        i.group(1) if isinstance(i, re.Match) else "No Description Available"
        for i in gene_details["Ensembl Gene Info"].apply(lambda x: re.search(r"\[Description: (.*)\]", x))
    ]

    gene_details["gene_info_type"] = [
        i.group(1) if isinstance(i, re.Match) else "None" for i in gene_details["Gene Info"].apply(lambda x: re.search(r"\[Gene Type: (.*)\]", x))
    ]
    gene_details["ensembl_info_type"] = [
        i.group(1) if isinstance(i, re.Match) else "None"
        for i in gene_details["Ensembl Gene Info"].apply(lambda x: re.search(r"\[Gene Type: (.*)\]", x))
    ]

    gene_type: list[str] = []
    row: pd.DataFrame
    for row in gene_details.itertuples():
        if row.gene_info_type != "None":
            gene_type.append(row.gene_info_type)
        elif row.ensembl_info_type != "None":
            gene_type.append(row.ensembl_info_type)
        else:
            gene_type.append("No Gene Type Available")
    gene_details["gene_type"] = gene_type

    # Drop gene_info_type and ensembl_info_type columns
    gene_details.drop(
        columns=[
            "Ensembl Gene Info",
            "Gene Info",
            "gene_info_type",
            "ensembl_info_type",
        ],
        inplace=True,
    )
    gene_details.rename(columns={"entrez_gene_id": "Entrez Gene ID"}, inplace=True)

    return gene_details


def _merge_xomics(
    context_name: str,
    expression_requirement,
    proteomics_file=None,
    trnaseq_file=None,
    mrnaseq_file=None,
    scrnaseq_file=None,
    no_hc=False,
    no_na=False,
):
    """
    Merges rnaseq and/or proteomics active gene logicals from outputs of their respective "_gen.py"
    scripts.

    :param proteomics_file: filename of proteomics config file in /main/data/config_sheets/
    :param trnaseq_file: filename of Total RNA-seq config file in /main/data/config_sheets/
    :param mrnaseq_file: filename of mRNA-seq config file in /main/data/config_sheets/
    :param scrnaseq_file: filename of single-cell RNA-seq config file in /main/data/config_sheets/
    :param no_hc: True if not adjusting for NA values (happens when gene missing from data source)
    :param no_na: filename of single-cell RNA-seq config file in /main/data/config_sheets/
    :param context_name: sheet name to use, should be context, context, cell type, etc
    :param expression_requirement: integer, minimum number of provided sources with active gene for a it to be in model
    :return: dictionary where keys are contexts, (tissue name, control type etc) and values are expression tables
    """
    config = Config()
    print(f"Merging data for {context_name}")
    # load data for each source if it exists. IF not load an empty dummy dataset
    trnaseq = rnaseq_gen.load_rnaseq_tests(filename=trnaseq_file, context_name=context_name, lib_type="total")
    mrnaseq = rnaseq_gen.load_rnaseq_tests(filename=mrnaseq_file, context_name=context_name, lib_type="mrna")
    scrnaseq = rnaseq_gen.load_rnaseq_tests(filename=scrnaseq_file, context_name=context_name, lib_type="scrna")
    proteomics = proteomics_gen.load_proteomics_tests(filename=proteomics_file, context_name=context_name)

    files_dict = dict()

    expression_list = []
    high_confidence_list = []
    merge_data = None

    if trnaseq[0] != "dummy":
        expression_list.append(_ExpressedHeaderNames.TRNASEQ)
        high_confidence_list.append(_HighExpressionHeaderNames.TRNASEQ)
        trnaseq_data = trnaseq[1].loc[:, ["expressed", "high"]]
        trnaseq_data.rename(
            columns={
                "expressed": _ExpressedHeaderNames.TRNASEQ,
                "high": _HighExpressionHeaderNames.TRNASEQ,
            },
            inplace=True,
        )
        merge_data = trnaseq_data

    if mrnaseq[0] != "dummy":
        expression_list.append(_ExpressedHeaderNames.MRNASEQ)
        high_confidence_list.append(_HighExpressionHeaderNames.MRNASEQ)
        mrnaseq_data = mrnaseq[1].loc[:, ["expressed", "high"]]
        mrnaseq_data.rename(
            columns={
                "expressed": _ExpressedHeaderNames.MRNASEQ,
                "high": _HighExpressionHeaderNames.MRNASEQ,
            },
            inplace=True,
        )
        merge_data = mrnaseq_data if merge_data is None else merge_data.join(mrnaseq_data, how="outer")

    if scrnaseq[0] != "dummy":
        expression_list.append(_ExpressedHeaderNames.SCRNASEQ)
        high_confidence_list.append(_HighExpressionHeaderNames.SCRNASEQ)
        scrnaseq_data = scrnaseq[1].loc[:, ["expressed", "high"]]
        scrnaseq_data.rename(
            columns={
                "expressed": _ExpressedHeaderNames.SCRNASEQ,
                "high": _HighExpressionHeaderNames.SCRNASEQ,
            },
            inplace=True,
        )
        merge_data = scrnaseq_data if merge_data is None else merge_data.join(scrnaseq_data, how="outer")

    if proteomics[0] != "dummy":
        expression_list.append(_ExpressedHeaderNames.PROTEOMICS)
        high_confidence_list.append(_HighExpressionHeaderNames.PROTEOMICS)
        prote_data = proteomics[1].loc[:, ["expressed", "high"]]
        prote_data.rename(
            columns={
                "expressed": _ExpressedHeaderNames.PROTEOMICS,
                "high": _HighExpressionHeaderNames.PROTEOMICS,
            },
            inplace=True,
        )
        merge_data = prote_data if merge_data is None else merge_data.join(prote_data, how="outer")

    merge_data = merge_logical_table(merge_data)

    num_sources = len(expression_list)
    merge_data["Active"] = 0
    merge_data["Required"] = 0

    if no_na:  # dont adjust for na values
        merge_data.loc[:, "Required"] = merge_data[expression_list].apply(
            lambda x: expression_requirement if (expression_requirement - (num_sources - x.count()) > 0) else 1,
            axis=1,
        )
    else:  # subtract one from requirement per NA
        merge_data.loc[:, "Required"] = merge_data[expression_list].apply(
            lambda x: expression_requirement - (num_sources - x.count()) if (expression_requirement - (num_sources - x.count()) > 0) else 1,
            axis=1,
        )

    # count number of sources gene is active in. Set to active in final output if at least adjusted expression reqirmnt
    merge_data["TotalExpressed"] = merge_data[expression_list].sum(axis=1)
    merge_data.loc[merge_data["TotalExpressed"] >= merge_data["Required"], "Active"] = 1

    if not no_hc:  # set genes that are high-confidence in at least one data source to active
        merge_data.loc[merge_data[high_confidence_list].sum(axis=1) > 0, "Active"] = 1

    # merge_data = merge_data.astype(int)
    merge_data = merge_data

    filepath = config.result_dir / context_name / f"merged_{context_name}.csv"
    merge_data.to_csv(filepath, index_label="ENTREZ_GENE_ID")

    filepath = config.result_dir / context_name / f"ActiveGenes_{context_name}_Merged.csv"
    merge_data.reset_index(drop=False, inplace=True)

    split_entrez = split_gene_expression_data(merge_data)
    split_entrez.rename(columns={"Gene": "ENTREZ_GENE_ID", "Data": "Active"}, inplace=True)
    split_entrez.to_csv(filepath, index_label="ENTREZ_GENE_ID")
    files_dict[context_name] = filepath.as_posix()  # Must return as string to be able to serialize

    transcriptomic_details = get_transcriptmoic_details(merge_data)
    transcriptomic_details_filepath = filepath.parent / f"TranscriptomicDetails_{context_name}.csv"
    transcriptomic_details.to_csv(transcriptomic_details_filepath, index=False)

    print(f"{context_name}: Save to {filepath}\n")

    return files_dict


def handle_context_batch(
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
    keep_gene_score,
):
    """
    Handle merging of different data sources for each context type
    """
    config = Config()
    sheet_names = []
    for file in [
        trnaseq_file,
        mrnaseq_file,
        scrnaseq_file,
        proteomics_file,
    ]:
        if file is not None:
            config_filepath = config.config_dir / file
            xl = pd.ExcelFile(config_filepath, engine="openpyxl")
            sheet_names += xl.sheet_names

    use_trna = True if trnaseq_file is not None else False
    use_mrna = True if mrnaseq_file is not None else False
    use_scrna = True if scrnaseq_file is not None else False
    use_proteins = True if proteomics_file is not None else False

    counts = Counter(sheet_names)
    sheet_names = sorted(list(set(sheet_names)))
    print("The data provided for each context listed will be merged. Data BETWEEN contexts will not be merged, only WITHIN a context")
    for i in sheet_names:
        # Print the sheet names in a list, like so
        # name1, name2, and name3
        print(i, end="")
        if counts[i] > 1:
            print(f" ({counts[i]}x)", end="")
        print(", ", end="")
    print("\b\b")

    dict_list = {}

    max_inputs = max(counts.values())
    min_inputs = min(counts.values())

    if merge_distro:
        print(f"Using {merge_distro} distribution for merging")
        rpy2_api.Rpy2(
            r_file_path,
            config.result_dir.as_posix(),
            sheet_names,
            use_mrna,
            use_trna,
            use_scrna,
            use_proteins,
            keep_gene_score,
            tweight,
            mweight,
            sweight,
            pweight,
        ).call_function("combine_zscores_main")

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

        print(f"Expression requirement of {expression_requirement} adjusted to {exp_req} using {adjust_method} adjustment method for {context_name}.")

        if exp_req > num_sources:
            print(
                f"WARNING: Expression requirement for {context_name} was calculated to be greater than max number of input data sources."
                f"Will be force changed to {num_sources} to prevent output from having 0 active genes. "
                f"Consider lowering the expression requirement or changing the adjustment method."
            )
            exp_req = num_sources

        if exp_req < 1:  # never allow expression requirement to be less than one
            print(
                "WARNING: Expression requirement for {context_name} was calculated to be less than 1. "
                "Will be force changed to 1 to prevent output from having 0 active genes. "
            )
            exp_req = 1

        files_dict = _merge_xomics(
            context_name,
            expression_requirement=exp_req,
            proteomics_file=proteomics_file,
            trnaseq_file=trnaseq_file,
            mrnaseq_file=mrnaseq_file,
            scrnaseq_file=scrnaseq_file,
            no_hc=no_hc,
            no_na=no_na,
        )

        dict_list.update(files_dict)

    files_json = config.result_dir / "step1_results_files.json"
    files_json.parent.mkdir(parents=True, exist_ok=True)
    with open(files_json.as_posix(), "w") as fp:
        json.dump(dict_list, fp)

    return


class AdjustMethod(Enum):
    PROGRESSIVE = "progressive"
    REGRESSIVE = "regressive"
    FLAT = "flat"
    CUSTOM = "custom"


def merge_xomics(
    trnaseq_file: str = None,
    mrnaseq_file: str = None,
    scrnaseq_file: str = None,
    proteomics_file: str = None,
    tweight: float = 1,
    mweight: float = 1,
    sweight: float = 1,
    pweight: float = 2,
    expression_requirement: int = None,
    adjust_method: AdjustMethod = AdjustMethod.FLAT,
    no_hc: bool = False,
    no_na: bool = False,
    custom_file: str = None,
    merge_distro: bool = False,
    keep_gene_score: bool = True,
):
    config = Config()
    # read custom expression requirment file if used
    if custom_file is not None:
        custom_filepath = config.data_dir / custom_file
        custom_df = pd.read_excel(custom_filepath, sheet_name=0)
        custom_df.columns = ["context", "req"]
    else:
        custom_df = pd.DataFrame([])

    if expression_requirement is None:
        expression_requirement = sum(
            1
            for test in [
                trnaseq_file,
                mrnaseq_file,
                scrnaseq_file,
                proteomics_file,
            ]
            if test is not None
        )
    elif expression_requirement < 1:
        raise ValueError("Expression requirement must be at least 1!")

    if adjust_method not in AdjustMethod:
        raise ValueError("Adjust method must be either 'progressive', 'regressive', 'flat', or 'custom'")

    handle_context_batch(
        trnaseq_file,
        mrnaseq_file,
        scrnaseq_file,
        proteomics_file,
        tweight,
        mweight,
        sweight,
        pweight,
        expression_requirement,
        adjust_method.value,
        no_hc,
        no_na,
        custom_df,
        merge_distro,
        keep_gene_score,
    )


def main():
    """
    Merge expression tables of multiple sources, (RNA-seq, proteomics) into one list
    User can specify the number of sources with an active gene in order for it to be considered active in the model.
    Otherwise, it defaults to the number of sources provided. High-confidence genes from any source will be considered
    active in the model, regardless of agreement with other sources.
    """
    parser = argparse.ArgumentParser(
        prog="merge_xomics.py",
        description="Merge expression tables of multiple sources (RNA-seq, proteomics) into one",
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
        "and treat as no expression",
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
        help="Number of sources with active gene for it to be considered active even if it is not a " "high confidence-gene",
    )

    parser.add_argument(
        "-r",
        "--requirement-adjust",
        type=str,
        required=False,
        default="flat",
        dest="adjust_method",
        help="Technique to adjust expression requirement based on differences in number of provided " "data source types.",
    )

    parser.add_argument(
        "-c",
        "--custom-requirement-file",
        required="custom" in sys.argv,  # required if --requriement-adjust is "custom",
        dest="custom_file",
        default="SKIP",
        help="Name of .xlsx file where first column is context names and second column is expression " "requirement for that context, in /main/data/",
    )

    parser.add_argument(
        "-hc",
        "--no-hc",
        action="store_true",
        required=False,
        default=False,
        dest="no_hc",
        help="Flag to prevent high-confidence genes forcing a gene to be used in final model " "irrespective of other other data sources",
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

    args = parser.parse_args()

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
    config = Config()

    # read custom expression requirment file if used
    if custom_file != "SKIP":
        custom_filepath = config.data_dir / custom_file
        custom_df = pd.read_excel(custom_filepath, sheet_name=0)
        custom_df.columns = ["context", "req"]
    else:
        custom_df = pd.DataFrame([])

    def_exp_req = sum(
        [
            1
            for test in [
                trnaseq_file,
                mrnaseq_file,
                scrnaseq_file,
                proteomics_file,
            ]
            if test is None
        ]
    )

    if expression_requirement.lower() == "default":
        expression_requirement = def_exp_req

    else:
        try:
            expression_requirement = int(expression_requirement)
            if int(expression_requirement) < 1:
                raise ValueError("Expression requirement must be at least 1!")
        except ValueError:
            raise ValueError("Expression requirement must be able to be converted to an integer!")

    if adjust_method not in ["progressive", "regressive", "flat", "custom"]:
        raise ValueError("Adjust method must be either 'progressive', 'regressive', 'flat', or 'custom'")

    handle_context_batch(
        trnaseq_file=trnaseq_file,
        mrnaseq_file=mrnaseq_file,
        scrnaseq_file=scrnaseq_file,
        proteomics_file=proteomics_file,
        tweight=tweight,
        mweight=mweight,
        sweight=sweight,
        pweight=pweight,
        expression_requirement=expression_requirement,
        adjust_method=adjust_method,
        no_hc=no_hc,
        no_na=no_na,
        custom_df=custom_df,
        merge_distro=merge_distro,
        keep_gene_score=keep_gene_score,
    )

    print("\nDone!")


if __name__ == "__main__":
    main()
