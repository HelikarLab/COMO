from __future__ import annotations

import asyncio
import sys
from io import TextIOWrapper
from pathlib import Path

import numpy as np
import pandas as pd
from fast_bioservices.biothings.mygene import MyGene
from loguru import logger

from como.combine_distributions import (
    _begin_combining_distributions,
)
from como.data_types import (
    AdjustmentMethod,
    LogLevel,
    RNAType,
    SourceTypes,
    _BatchEntry,
    _BatchNames,
    _InputMatrices,
    _OutputCombinedSourceFilepath,
    _SourceWeights,
)
from como.project import Config
from como.utils import _log_and_raise_error, _read_file, _set_up_logging, get_missing_gene_data, return_placeholder_data


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


# TODO: If function is no longer needed, remove?
async def _load_rnaseq_tests(filename, context_name, prep_method: RNAType) -> tuple[str, pd.DataFrame]:
    """Load rnaseq results."""
    logger.debug(f"Loading data for context '{context_name}' using preparation method '{prep_method.value}'")
    config = Config()

    def load_dummy_dict():
        df = return_placeholder_data()
        return "dummy", df

    if not filename or filename == "None":  # not using this data type, use empty dummy data matrix
        return load_dummy_dict()

    inquiry_full_path = config.data_dir / "config_sheets" / filename
    if not inquiry_full_path.exists():
        _log_and_raise_error(
            f"Config file not found at {inquiry_full_path}",
            error=FileNotFoundError,
            level=LogLevel.ERROR,
        )

    match prep_method:
        case RNAType.TRNA:
            filename = f"{RNAType.TRNA.value}_{context_name}.csv"
        case RNAType.MRNA:
            filename = f"{RNAType.MRNA.value}_{context_name}.csv"
        case RNAType.SCRNA:
            filename = f"{RNAType.SCRNA.value}_{context_name}.csv"
        case _:
            _log_and_raise_error(
                f"Unsupported RNA-seq library type: {prep_method.value}. Must be one of {', '.join(RNAType)}.",
                error=ValueError,
                level=LogLevel.ERROR,
            )

    save_filepath = config.result_dir / context_name / prep_method.value / filename
    if save_filepath.exists():
        logger.debug(f"Loading RNA-seq data from: {save_filepath}")
        data = pd.read_csv(save_filepath, index_col="entrez_gene_id")
        logger.success(f"Successfully loaded RNA-seq data from: {save_filepath}")
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


async def _get_transcriptmoic_details(merged_df: pd.DataFrame, taxon_id: int) -> pd.DataFrame:
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
    transcriptomic_df: pd.DataFrame = merged_df.copy()
    if _ExpressedHeaderNames.PROTEOMICS in merged_df.columns:
        # Get the number of sources required for a gene to be marked "expressed"
        required_expression = merged_df["required"].iloc[0]

        # Subtract 1 from merged_df["TotalExpressed"] if the current value is greater than or equal to 1
        # This is done to take into account the removal of proteomic expression
        merged_df["total_expressed"] = merged_df["total_expressed"].apply(lambda x: x - 1 if x >= 1 else x)

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
            transcriptomic_df["total_expressed"] >= transcriptomic_df["required"],
            "active",
        ] = 1

    my_gene = MyGene()
    gene_details: pd.DataFrame = pd.DataFrame(
        data=pd.NA,
        columns=["entrez_gene_id", "gene_symbol", "description", "gene_type"],
        index=list(range(len(transcriptomic_df))),
    )
    for i, detail in enumerate(
        await my_gene.query(
            items=transcriptomic_df["entrez_gene_id"].tolist(),
            taxon=taxon_id,
            scopes="entrezgene",
        )
    ):
        gene_details.at[i, "entrez_gene_id"] = detail["entrezgene"]
        gene_details.at[i, "gene_symbol"] = detail["symbol"]
        gene_details.at[i, "description"] = detail["name"]
        gene_details.at[i, "gene_type"] = detail["type_of_gene"]
    return gene_details


async def _merge_xomics(
    context_name: str,
    expression_requirement: int,
    trna_boolean_matrix: pd.DataFrame | None,
    mrna_boolean_matrix: pd.DataFrame | None,
    scrna_boolean_matrix: pd.DataFrame | None,
    proteomic_boolean_matrix: pd.DataFrame | None,
    output_merged_filepath: Path,
    output_gene_activity_filepath: Path,
    output_transcriptomic_details_filepath: Path,
    taxon_id: int,
    force_activate_high_confidence: bool = True,
    adjust_for_missing_sources: bool = False,
):
    expression_list: list[str] = []
    high_confidence_list: list[str] = []
    merge_data = pd.DataFrame()

    for matrix, expressed_sourcetype, high_expressed_sourcetype in (
        (trna_boolean_matrix, _ExpressedHeaderNames.TRNASEQ, _HighExpressionHeaderNames.TRNASEQ),
        (mrna_boolean_matrix, _ExpressedHeaderNames.MRNASEQ, _HighExpressionHeaderNames.MRNASEQ),
        (scrna_boolean_matrix, _ExpressedHeaderNames.SCRNASEQ, _HighExpressionHeaderNames.SCRNASEQ),
        (proteomic_boolean_matrix, _ExpressedHeaderNames.PROTEOMICS, _HighExpressionHeaderNames.PROTEOMICS),
    ):
        if matrix is None:
            continue

        matrix: pd.DataFrame  # re-define type to assist in type hinting for IDEs
        expression_list.append(expressed_sourcetype)
        high_confidence_list.append(high_expressed_sourcetype)
        matrix.rename(columns={"expressed": expressed_sourcetype, "high": high_expressed_sourcetype}, inplace=True)
        merge_data = matrix if merge_data.empty else merge_data.merge(matrix, on="entrez_gene_id", how="outer")

    if merge_data.empty:
        logger.critical(
            f"No data is available for the '{context_name}' context. If this is intentional, ignore this error."
        )
        return {}

    merge_data = _merge_logical_table(merge_data)
    num_sources = len(expression_list)
    merge_data["active"] = 0
    merge_data["required"] = 0

    if adjust_for_missing_sources:  # Subtract 1 from requirement per missing source
        merge_data.loc[:, "required"] = merge_data[expression_list].apply(
            lambda x: expression_requirement - (num_sources - x.count())
            if (expression_requirement - (num_sources - x.count()) > 0)
            else 1,
            axis=1,
        )
    else:  # Do not adjust for missing sources
        merge_data.loc[:, "required"] = merge_data[expression_list].apply(
            lambda x: expression_requirement if (expression_requirement - (num_sources - x.count()) > 0) else 1, axis=1
        )

    # Count the number of sources each gene is active in
    # set to active in final output if we meet the adjusted expression requirement
    merge_data["total_expressed"] = merge_data[expression_list].sum(axis=1)
    merge_data.loc[merge_data["total_expressed"] >= merge_data["required"], "active"] = 1

    if force_activate_high_confidence:  # If a gene is high-confidence in at least 1 data source, set it to active
        merge_data.loc[merge_data[high_confidence_list].sum(axis=1) > 0, "active"] = 1

    merge_data.to_csv(output_merged_filepath, index=False)
    transcriptomic_details = await _get_transcriptmoic_details(merge_data, taxon_id=taxon_id)
    transcriptomic_details.to_csv(output_transcriptomic_details_filepath, index=False)

    return {context_name: output_gene_activity_filepath.as_posix()}


async def _process(
    *,
    context_name: str,
    trna_matrix: pd.DataFrame | None,
    mrna_matrix: pd.DataFrame | None,
    scrna_matrix: pd.DataFrame | None,
    proteomic_matrix: pd.DataFrame | None,
    trna_boolean_matrix: pd.DataFrame | None,
    mrna_boolean_matrix: pd.DataFrame | None,
    scrna_boolean_matrix: pd.DataFrame | None,
    proteomic_boolean_matrix: pd.DataFrame | None,
    trna_batches: dict[int, list[str]],
    mrna_batches: dict[int, list[str]],
    scrna_batches: dict[int, list[str]],
    proteomic_batches: dict[int, list[str]],
    trna_weight: int,
    mrna_weight: int,
    scrna_weight: int,
    proteomic_weight: int,
    taxon_id: int,
    minimum_source_expression: int,
    expression_requirement: int,
    weighted_z_floor: int,
    weighted_z_ceiling: int,
    adjust_method: AdjustmentMethod,
    merge_zfpkm_distribution: bool,
    force_activate_high_confidence: bool,
    adjust_for_missing_sources: bool,
    output_merge_activity_filepath: Path,
    output_transcriptomic_details_filepath: Path,
    output_trna_activity_filepath: Path | None,
    output_mrna_activity_filepath: Path | None,
    output_scrna_activity_filepath: Path | None,
    output_proteomic_activity_filepath: Path | None,
    output_final_model_scores_filepath: Path,
    output_figure_dirpath: Path | None,
):
    """Merge different data sources for each context type."""
    num_sources = sum(1 for source in [trna_matrix, mrna_matrix, scrna_matrix, proteomic_matrix] if source is not None)

    if merge_zfpkm_distribution:
        _combine_zscores(
            context_name=context_name,
            input_matrices=_InputMatrices(
                trna=trna_matrix,
                mrna=mrna_matrix,
                scrna=scrna_matrix,
                proteomics=proteomic_matrix,
            ),
            batch_names=_BatchNames(
                trna=[_BatchEntry(batch_num=n, sample_names=s) for n, s in trna_batches.items()],
                mrna=[_BatchEntry(batch_num=n, sample_names=s) for n, s in mrna_batches.items()],
                scrna=[_BatchEntry(batch_num=n, sample_names=s) for n, s in scrna_batches.items()],
                proteomics=[_BatchEntry(batch_num=n, sample_names=s) for n, s in proteomic_batches.items()],
            ),
            source_weights=_SourceWeights(
                trna=trna_weight,
                mrna=mrna_weight,
                scrna=scrna_weight,
                proteomics=proteomic_weight,
            ),
            output_filepaths=_OutputCombinedSourceFilepath(
                trna=output_trna_activity_filepath,
                mrna=output_mrna_activity_filepath,
                scrna=output_scrna_activity_filepath,
                proteomics=output_proteomic_activity_filepath,
            ),
            output_figure_dirpath=output_figure_dirpath,
            output_final_model_scores=output_final_model_scores_filepath,
            weighted_z_floor=weighted_z_floor,
            weighted_z_ceiling=weighted_z_ceiling,
        )

    # the more data sources available, the higher the expression requirement for the gene
    if adjust_method == AdjustmentMethod.PROGRESSIVE:
        adjusted_expression_requirement = (num_sources - minimum_source_expression) + expression_requirement
    # the more data sources available, the lower the expression requirement for the gene
    elif adjust_method == AdjustmentMethod.REGRESSIVE:
        # we use a hardcoded 4 here because that is the maximum number of contexts available
        # (trna, mrna, scrna, and proteomics is 4 sources)
        adjusted_expression_requirement = expression_requirement - (4 - num_sources)
    elif adjust_method == AdjustmentMethod.FLAT:
        adjusted_expression_requirement = expression_requirement

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

    await _merge_xomics(
        context_name=context_name,
        expression_requirement=adjusted_expression_requirement,
        trna_boolean_matrix=trna_boolean_matrix,
        mrna_boolean_matrix=mrna_boolean_matrix,
        scrna_boolean_matrix=scrna_boolean_matrix,
        proteomic_boolean_matrix=proteomic_boolean_matrix,
        output_merged_filepath=output_merge_activity_filepath,
        output_gene_activity_filepath=output_final_model_scores_filepath,
        output_transcriptomic_details_filepath=output_transcriptomic_details_filepath,
        taxon_id=taxon_id,
        force_activate_high_confidence=force_activate_high_confidence,
        adjust_for_missing_sources=adjust_for_missing_sources,
    )


async def merge_xomics(  # noqa: C901
    context_name: str,
    trna_matrix_or_filepath: Path | pd.DataFrame | None,
    mrna_matrix_or_filepath: Path | pd.DataFrame | None,
    scrna_matrix_or_filepath: Path | pd.DataFrame | None,
    proteomic_matrix_or_filepath: Path | pd.DataFrame | None,
    trna_boolean_matrix_or_filepath: Path | pd.DataFrame | None,
    mrna_boolean_matrix_or_filepath: Path | pd.DataFrame | None,
    scrna_boolean_matrix_or_filepath: Path | pd.DataFrame | None,
    proteomic_boolean_matrix_or_filepath: Path | pd.DataFrame | None,
    trna_batches: dict[int, list[str]],
    mrna_batches: dict[int, list[str]],
    scrna_batches: dict[int, list[str]],
    proteomic_batches: dict[int, list[str]],
    output_merge_activity_filepath: Path,
    output_transcriptomic_details_filepath: Path,
    output_trna_activity_filepath: Path | None,
    output_mrna_activity_filepath: Path | None,
    output_scrna_activity_filepath: Path | None,
    output_proteomic_activity_filepath: Path | None,
    output_final_model_scores_filepath: Path,
    output_figure_dirpath: Path | None,
    taxon_id: int,
    trna_weight: int = 1,
    mrna_weight: int = 1,
    scrna_weight: int = 1,
    proteomic_weight: int = 2,
    minimum_source_expression: int = 1,
    expression_requirement: int | None = None,
    adjust_method: AdjustmentMethod = AdjustmentMethod.FLAT,
    force_activate_high_confidence: bool = False,
    adjust_for_na: bool = False,
    merge_zfpkm_distribution: bool = False,
    weighted_z_floor: int = -6,
    weighted_z_ceiling: int = 6,
):
    """Merge expression tables of multiple sources (RNA-seq, proteomics) into one."""
    if expression_requirement < 1:
        raise ValueError("Expression requirement must be at least 1!")

    if expression_requirement is None:
        expression_requirement = sum(
            test is not None
            for test in [
                trna_matrix_or_filepath,
                mrna_matrix_or_filepath,
                scrna_matrix_or_filepath,
                proteomic_matrix_or_filepath,
            ]
        )

    if all(
        file is None
        for file in [
            trna_matrix_or_filepath,
            mrna_matrix_or_filepath,
            scrna_matrix_or_filepath,
            proteomic_matrix_or_filepath,
        ]
    ):
        raise ValueError("No data was passed!")

    if adjust_method not in AdjustmentMethod:
        raise ValueError(f"Adjustment method must be one of {AdjustmentMethod}; got: {adjust_method}")

    output_final_model_scores_filepath.parent.mkdir(parents=True, exist_ok=True)
    if output_merge_activity_filepath:
        output_merge_activity_filepath.parent.mkdir(parents=True, exist_ok=True)
    if output_transcriptomic_details_filepath:
        output_transcriptomic_details_filepath.parent.mkdir(parents=True, exist_ok=True)
    if output_trna_activity_filepath:
        output_trna_activity_filepath.parent.mkdir(parents=True, exist_ok=True)
    if output_mrna_activity_filepath:
        output_mrna_activity_filepath.parent.mkdir(parents=True, exist_ok=True)
    if output_scrna_activity_filepath:
        output_scrna_activity_filepath.parent.mkdir(parents=True, exist_ok=True)
    if output_proteomic_activity_filepath:
        output_proteomic_activity_filepath.parent.mkdir(parents=True, exist_ok=True)
    if output_figure_dirpath:
        output_figure_dirpath.mkdir(parents=True, exist_ok=True)

    # fmt: off
    trna_matrix: pd.DataFrame | None = (
        pd.read_csv(trna_matrix_or_filepath) if isinstance(trna_matrix_or_filepath, Path)
        else trna_matrix_or_filepath if isinstance(trna_matrix_or_filepath, pd.DataFrame)
        else None
    )
    mrna_matrix: pd.DataFrame | None = (
        pd.read_csv(mrna_matrix_or_filepath) if isinstance(mrna_matrix_or_filepath, Path)
        else mrna_matrix_or_filepath if isinstance(mrna_matrix_or_filepath, pd.DataFrame)
        else None
    )
    scrna_matrix: pd.DataFrame | None = (
        pd.read_csv(scrna_matrix_or_filepath) if isinstance(scrna_matrix_or_filepath, Path)
        else scrna_matrix_or_filepath if isinstance(scrna_matrix_or_filepath, pd.DataFrame)
        else None
    )
    proteomic_matrix: pd.DataFrame | None = (
        pd.read_csv(proteomic_matrix_or_filepath) if isinstance(proteomic_matrix_or_filepath, Path)
        else proteomic_matrix_or_filepath if isinstance(proteomic_matrix_or_filepath, pd.DataFrame)
        else None
    )

    trna_boolean_matrix: pd.DataFrame | None = (
        pd.read_csv(trna_boolean_matrix_or_filepath) if isinstance(trna_boolean_matrix_or_filepath, Path)
        else trna_boolean_matrix_or_filepath if isinstance(trna_boolean_matrix_or_filepath, pd.DataFrame)
        else None
    )
    mrna_boolean_matrix: pd.DataFrame | None = (
        pd.read_csv(mrna_boolean_matrix_or_filepath) if isinstance(mrna_boolean_matrix_or_filepath, Path)
        else mrna_boolean_matrix_or_filepath if isinstance(mrna_boolean_matrix_or_filepath, pd.DataFrame)
        else None
    )
    scrna_boolean_matrix: pd.DataFrame | None = (
        pd.read_csv(scrna_boolean_matrix_or_filepath) if isinstance(scrna_boolean_matrix_or_filepath, Path)
        else scrna_boolean_matrix_or_filepath if isinstance(scrna_boolean_matrix_or_filepath, pd.DataFrame)
        else None
    )
    proteomic_boolean_matrix: pd.DataFrame | None = (
        pd.read_csv(proteomic_boolean_matrix_or_filepath) if isinstance(proteomic_boolean_matrix_or_filepath, Path)
        else proteomic_boolean_matrix_or_filepath if isinstance(proteomic_boolean_matrix_or_filepath, pd.DataFrame)
        else None
    )
    # fmt: on

    await _process(
        context_name=context_name,
        trna_matrix=trna_matrix,
        mrna_matrix=mrna_matrix,
        scrna_matrix=scrna_matrix,
        proteomic_matrix=proteomic_matrix,
        trna_boolean_matrix=trna_boolean_matrix,
        mrna_boolean_matrix=mrna_boolean_matrix,
        scrna_boolean_matrix=scrna_boolean_matrix,
        proteomic_boolean_matrix=proteomic_boolean_matrix,
        trna_batches=trna_batches,
        mrna_batches=mrna_batches,
        scrna_batches=scrna_batches,
        proteomic_batches=proteomic_batches,
        trna_weight=trna_weight,
        mrna_weight=mrna_weight,
        scrna_weight=scrna_weight,
        proteomic_weight=proteomic_weight,
        taxon_id=taxon_id,
        minimum_source_expression=minimum_source_expression,
        expression_requirement=expression_requirement,
        weighted_z_floor=weighted_z_floor,
        weighted_z_ceiling=weighted_z_ceiling,
        adjust_method=adjust_method,
        merge_zfpkm_distribution=merge_zfpkm_distribution,
        force_activate_high_confidence=force_activate_high_confidence,
        adjust_for_missing_sources=adjust_for_na,
        output_merge_activity_filepath=output_merge_activity_filepath,
        output_transcriptomic_details_filepath=output_transcriptomic_details_filepath,
        output_trna_activity_filepath=output_trna_activity_filepath,
        output_mrna_activity_filepath=output_mrna_activity_filepath,
        output_scrna_activity_filepath=output_scrna_activity_filepath,
        output_proteomic_activity_filepath=output_proteomic_activity_filepath,
        output_final_model_scores_filepath=output_final_model_scores_filepath,
        output_figure_dirpath=output_figure_dirpath,
    )
