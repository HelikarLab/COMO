from __future__ import annotations

import argparse
import json
import re
import sys
from collections import Counter
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

import numpy as np
import pandas as pd
from fast_bioservices import BioDBNet, Input, Output
from loguru import logger

from como import proteomics_gen, return_placeholder_data
from como.combine_distributions import _combine_zscores
from como.project import Config
from como.types import RNAPrepMethod
from como.utils import split_gene_expression_data


class _MergedHeaderNames:
    TRNASEQ = "trnaseq"
    MRNASEQ = "mrnaseq"
    SCRNASEQ = "scrnaseq"
    PROTEOMICS = "prote"


class _ExpressedHeaderNames:
    TRNASEQ = f"{_MergedHeaderNames.TRNASEQ}_exp"
    MRNASEQ = f"{_MergedHeaderNames.MRNASEQ}_exp"
    SCRNASEQ = f"{_MergedHeaderNames.SCRNASEQ}_exp"
    PROTEOMICS = f"{_MergedHeaderNames.PROTEOMICS}_exp"


class _HighExpressionHeaderNames:
    TRNASEQ = f"{_MergedHeaderNames.TRNASEQ}_high"
    MRNASEQ = f"{_MergedHeaderNames.MRNASEQ}_high"
    SCRNASEQ = f"{_MergedHeaderNames.SCRNASEQ}_high"
    PROTEOMICS = f"{_MergedHeaderNames.PROTEOMICS}_high"


class AdjustmentMethod(Enum):
    """Adjustment method for expression requirement based on differences in number of provided data source types."""

    PROGRESSIVE = "progressive"
    REGRESSIVE = "regressive"
    FLAT = "flat"
    CUSTOM = "custom"

    @classmethod
    def from_string(cls, value: str) -> AdjustmentMethod:
        """Convert a string to an AdjustmentMethod enum."""
        if value.lower() not in [t.value for t in cls]:
            raise ValueError(f"Adjustment method must be one of {cls}; got: {value}")
        return cls(value)


def _load_rnaseq_tests(filename, context_name, prep_method: RNAPrepMethod) -> tuple[str, pd.DataFrame]:
    """Load rnaseq results.

    Returns a dictionary of test (context, context, cell, etc ) names and rnaseq expression data
    """
    config = Config()

    def load_dummy_dict():
        df = return_placeholder_data()
        return "dummy", df

    if not filename or filename == "None":  # not using this data type, use empty dummy data matrix
        return load_dummy_dict()

    inquiry_full_path = config.data_dir / "config_sheets" / filename
    if not inquiry_full_path.exists():  # check that config file exist (isn't needed to load, but helps user)
        raise FileNotFoundError(f"Error: Config file not found at {inquiry_full_path}")

    match prep_method:
        case RNAPrepMethod.TOTAL:
            filename = f"rnaseq_total_{context_name}.csv"
        case RNAPrepMethod.MRNA:
            filename = f"rnaseq_mrna_{context_name}.csv"
        case RNAPrepMethod.SCRNA:
            filename = f"rnaseq_scrna_{context_name}.csv"
        case _:
            raise ValueError(
                f"Unsupported RNA-seq library type: {prep_method.value}. "
                f"Must be an option defined in 'RNASeqPreparationMethod'."
            )

    save_filepath = config.result_dir / context_name / prep_method.value / filename
    if save_filepath.exists():
        data = pd.read_csv(save_filepath, index_col="entrez_gene_id")
        return context_name, data

    else:
        logger.warning(
            f"'{prep_method.value}' gene expression file for '{context_name}' was not found at '{save_filepath}'. "
            f"If this is not intentional, please fix the filename to match '{save_filepath}'."
        )
        return load_dummy_dict()


# Merge Output
def _merge_logical_table(df: pd.DataFrame):
    """Merge rows of Logical Table belonging to the same entrez_gene_id.

    :param df: pandas dataframe of logical table
    :return: pandas dataframe of merged table
    """
    # step 1: get all plural ENTREZ_GENE_IDs in the input table, extract unique IDs

    df.reset_index(drop=False, inplace=True)
    df.dropna(axis=0, subset=["entrez_gene_id"], inplace=True)
    df["entrez_gene_id"] = df["entrez_gene_id"].astype(str).str.replace(" /// ", "//").astype(str)

    id_list: list[str] = df.loc[
        ~df["entrez_gene_id"].str.contains("//"), "entrez_gene_id"
    ].tolist()  # Collect "single" ids, like "123"
    multiple_entrez_ids: list[str] = df.loc[
        df["entrez_gene_id"].str.contains("//"), "entrez_gene_id"
    ].tolist()  # Collect "double" ids, like "123//456"

    for i in multiple_entrez_ids:
        ids = i.split("//")
        id_list.extend(ids)

        duplicate_rows = pd.DataFrame([])
        for j in ids:
            rows = df.loc[df["entrez_gene_id"] == i].copy()
            rows["entrez_gene_id"] = j
            duplicate_rows = pd.concat([duplicate_rows, rows], axis=0)

        df = pd.concat([df, pd.DataFrame(duplicate_rows)], axis=0, ignore_index=True)
        df.drop(df[df["entrez_gene_id"] == i].index, inplace=True)

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

    entrez_dups_list.extend(i.split(" /// ") for i in full_entrez_id_sets)
    entrez_dups_dict = dict(zip(full_entrez_id_sets, entrez_dups_list))

    for merged_entrez_id, entrez_dups_list in entrez_dups_dict.items():
        df["entrez_gene_id"].replace(to_replace=entrez_dups_list, value=merged_entrez_id, inplace=True)

    df.set_index("entrez_gene_id", inplace=True)
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


async def _get_transcriptmoic_details(merged_df: pd.DataFrame) -> pd.DataFrame:
    """Get details of transcriptomic data.

    This function will get the following details of transcriptomic data:
    - Gene Symbol
    - Gene Name
    - entrez_gene_id

    The resulting dataframe will have its columns created in the order listed above
    It will return a pandas dataframe with this information

    :param merged_df: A dataframe containing all active transcriptomic and proteomic genes
    :return: A dataframe with the above-listed columns
    """
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
    gene_details: pd.DataFrame = await biodbnet.async_db2db(
        values=transcriptomic_df.index.astype(str).values.tolist(),
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
        i.group(1) if isinstance(i, re.Match) else "None"
        for i in gene_details["Gene Info"].apply(lambda x: re.search(r"\[Gene Type: (.*)\]", x))
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


async def _merge_xomics(
    context_name: str,
    expression_requirement,
    proteomics_file=None,
    trnaseq_file=None,
    mrnaseq_file=None,
    scrnaseq_file=None,
    no_hc=False,
    no_na=False,
):
    """Merge rnaseq and/or proteomics active genes.

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
    logger.info(f"Merging data for {context_name}")
    # load data for each source if it exists. IF not load an empty dummy dataset
    trnaseq = _load_rnaseq_tests(filename=trnaseq_file, context_name=context_name, prep_method=RNAPrepMethod.TOTAL)
    mrnaseq = _load_rnaseq_tests(filename=mrnaseq_file, context_name=context_name, prep_method=RNAPrepMethod.MRNA)
    scrnaseq = _load_rnaseq_tests(filename=scrnaseq_file, context_name=context_name, prep_method=RNAPrepMethod.SCRNA)
    proteomics = proteomics_gen.load_proteomics_tests(filename=proteomics_file, context_name=context_name)

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

    if merge_data is None:
        logger.critical(
            f"No data is available for the '{context_name}' context. If this is intentional, ignore this error."
        )
        return {}
    merge_data = _merge_logical_table(merge_data)

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
            lambda x: expression_requirement - (num_sources - x.count())
            if (expression_requirement - (num_sources - x.count()) > 0)
            else 1,
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
    merge_data.to_csv(filepath, index_label="entrez_gene_id")

    filepath = config.result_dir / context_name / f"ActiveGenes_{context_name}_Merged.csv"
    merge_data.reset_index(drop=False, inplace=True)

    split_entrez = split_gene_expression_data(merge_data)
    split_entrez.to_csv(filepath, index_label="entrez_gene_id")
    files_dict = {context_name: filepath.as_posix()}

    transcriptomic_details = await _get_transcriptmoic_details(merge_data)
    transcriptomic_details_filepath = filepath.parent / f"TranscriptomicDetails_{context_name}.csv"
    transcriptomic_details.to_csv(transcriptomic_details_filepath, index=False)

    logger.success(f"{context_name}: Save to {filepath}\n")

    return files_dict


async def _handle_context_batch(  # noqa: C901
    trnaseq_file: Path | None,
    mrnaseq_file: Path | None,
    scrnaseq_file: Path | None,
    proteomics_file: Path | None,
    tweight: float,
    mweight: float,
    sweight: float,
    pweight: float,
    expression_requirement: int,
    adjust_method: AdjustmentMethod,
    no_hc: bool,
    no_na: bool,
    merge_zfpkm_distribution: bool,
    keep_gene_score: bool,
):
    """Merge different data sources for each context type."""
    if all(file is None for file in [trnaseq_file, mrnaseq_file, scrnaseq_file, proteomics_file]):
        raise ValueError("No configuration file was passed!")

    config = Config()
    sheet_names = []
    for file in [trnaseq_file, mrnaseq_file, scrnaseq_file, proteomics_file]:
        if file is not None:
            config_filepath = config.config_dir / file
            try:
                xl = pd.ExcelFile(config_filepath, engine="openpyxl")
            except Exception as e:
                raise ValueError(f"Unable to read file '{config_filepath}'") from e
            sheet_names += xl.sheet_names

    use_trna = trnaseq_file is not None
    use_mrna = mrnaseq_file is not None
    use_scrna = scrnaseq_file is not None
    use_proteins = proteomics_file is not None

    counts = Counter(sheet_names)
    sheet_names = sorted(set(sheet_names))
    logger.info("Beginning to merge data within contexts")
    dict_list = {}

    max_inputs = max(counts.values())
    min_inputs = min(counts.values())

    if merge_zfpkm_distribution:
        logger.debug("Using zFPKM distribution for merging")
        _combine_zscores(
            working_dir=config.result_dir.as_posix(),
            context_names=sheet_names,
            global_use_mrna=use_mrna,
            global_use_trna=use_trna,
            global_use_scrna=use_scrna,
            global_use_proteins=use_proteins,
            keep_gene_scores=keep_gene_score,
            global_trna_weight=tweight,
            global_mrna_weight=mweight,
            global_scrna_weight=sweight,
            global_protein_weight=pweight,
        )

    for context_name in sheet_names:
        num_sources = counts[context_name]
        match adjust_method:
            case AdjustmentMethod.PROGRESSIVE:
                adjusted_expression_requirement = (num_sources - min_inputs) + expression_requirement
            case AdjustmentMethod.REGRESSIVE:
                adjusted_expression_requirement = expression_requirement - (max_inputs - num_sources)
            case AdjustmentMethod.FLAT:
                adjusted_expression_requirement = expression_requirement
            case _:
                adjusted_expression_requirement = int(
                    custom_df.iloc[custom_df["context"] == context_name, "req"].iloc[0]
                )

        if adjusted_expression_requirement != expression_requirement:
            logger.debug(
                f"Expression requirement of '{expression_requirement}' adjusted to "
                f"'{adjusted_expression_requirement}' using '{adjust_method.value}' adjustment method "
                f"for '{context_name}'."
            )

        if adjusted_expression_requirement > num_sources:
            logger.warning(
                f"Expression requirement for {context_name} was calculated to be greater "
                f"than max number of input data sources. "
                f"Will be force changed to {num_sources} to prevent output from having 0 active genes. "
                f"Consider lowering the expression requirement or changing the adjustment method."
            )
            adjusted_expression_requirement = num_sources

        if adjusted_expression_requirement < 1:  # never allow expression requirement to be less than one
            logger.warning(
                f"Expression requirement for {context_name} was calculated to be less than 1. "
                "Will be changed to 1 to prevent output from having 0 active genes. "
            )
            adjusted_expression_requirement = 1

        files_dict = await _merge_xomics(
            context_name,
            expression_requirement=adjusted_expression_requirement,
            proteomics_file=proteomics_file,
            trnaseq_file=trnaseq_file,
            mrnaseq_file=mrnaseq_file,
            scrnaseq_file=scrnaseq_file,
            no_hc=no_hc,
            no_na=no_na,
        )

        dict_list |= files_dict

    files_json = config.result_dir / "step1_results_files.json"
    files_json.parent.mkdir(parents=True, exist_ok=True)
    with files_json.open("w") as f:
        json.dump(dict_list, f)  # type: ignore

    return


async def merge_xomics(
    trnaseq_filepath: str | Path | None = None,
    mrnaseq_filepath: str | Path | None = None,
    scrnaseq_filepath: str | Path | None = None,
    proteomics_filepath: str | Path | None = None,
    trna_weight: float = 1,
    mrna_weight: float = 1,
    scrna_weight: float = 1,
    proteomics_weight: float = 2,
    expression_requirement: int | None = None,
    adjust_method: AdjustmentMethod = AdjustmentMethod.FLAT,
    no_high_confidence: bool = False,
    no_na: bool = False,
    merge_zfpkm_distribution: bool = False,
    keep_transcriptomics_score: bool = True,
):
    """Merge expression tables of multiple sources (RNA-seq, proteomics) into one."""
    if expression_requirement is None:
        expression_requirement = sum(
            test is not None
            for test in [
                trnaseq_filepath,
                mrnaseq_filepath,
                scrnaseq_filepath,
                proteomics_filepath,
            ]
        )
    elif expression_requirement < 1:
        raise ValueError("Expression requirement must be at least 1!")

    trnaseq_filepath = Path(trnaseq_filepath) if trnaseq_filepath else None
    mrnaseq_filepath = Path(mrnaseq_filepath) if mrnaseq_filepath else None
    scrnaseq_filepath = Path(scrnaseq_filepath) if scrnaseq_filepath else None
    proteomics_filepath = Path(proteomics_filepath) if proteomics_filepath else None

    await _handle_context_batch(
        trnaseq_filepath,
        mrnaseq_filepath,
        scrnaseq_filepath,
        proteomics_filepath,
        trna_weight,
        mrna_weight,
        scrna_weight,
        proteomics_weight,
        expression_requirement,
        adjust_method,
        no_high_confidence,
        no_na,
        merge_zfpkm_distribution,
        keep_transcriptomics_score,
    )
    )
    )
    )

    )
    )
    )
    )
    )
