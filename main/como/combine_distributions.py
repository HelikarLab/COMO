from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from loguru import logger

from como.data_types import (
    GeneIdentifier,
    SourceTypes,
    _BatchEntry,
    _BatchNames,
    _CombineOmicsInput,
    _InputMatrices,
    _OutputCombinedSourceFilepath,
    _SourceWeights,
)
from como.pipelines.identifier import get_remaining_identifiers
from como.utils import num_columns


def _combine_z_distribution_for_batch(
    context_name: str,
    batch: _BatchEntry,
    matrix: pd.DataFrame,
    source: SourceTypes,
    output_combined_matrix_filepath: Path,
    output_figure_dirpath: Path,
) -> pd.DataFrame:
    """Combine z-score distributions across samples for a single batch.

    :param context_name: Name of the context (e.g., tissue or condition).
    :param batch: Batch entry containing batch number and sample names.
    :param matrix: DataFrame with 'ensembl_gene_id' and sample columns.
    :param source: Source type (e.g., trna, mrna, scrna, proteomics).
    :param output_combined_matrix_filepath: Path to save the combined z-score matrix.
    :param output_figure_dirpath: Path to save the z-score distribution figure.

    :returns: A pandas dataframe of the weighted z-distributions
    """
    output_combined_matrix_filepath.parent.mkdir(parents=True, exist_ok=True)
    output_figure_dirpath.mkdir(parents=True, exist_ok=True)

    logger.trace(
        f"Combining z-score distributions: batch #{batch.batch_num}, "
        f"samples: {len(batch.sample_names)}, "
        f"source: '{source.value}', "
        f"context: '{context_name}'"
    )
    if num_columns(matrix) < 2:
        logger.trace(
            f"A single sample exists for batch '{batch.batch_num}'. "
            f"Returning as-is because no additional combining can be done"
        )
        with_batch_num = matrix.copy()
        with_batch_num.columns = [batch.batch_num]
        return with_batch_num

    values = matrix.copy().replace([np.inf, -np.inf], np.nan)
    values_array = values.to_numpy(dtype=float)
    mask = ~np.isnan(values_array)
    weights = np.ones(values.shape[1], dtype=float)
    masked_values = np.where(mask, values_array * weights, 0.0)
    weighted_sum = masked_values.sum(axis=1)

    weights_squared = weights**2
    sum_weights_squared = (mask * weights_squared).sum(axis=1)
    denominator = np.sqrt(sum_weights_squared)
    combined_z = np.divide(
        weighted_sum,
        denominator,
        out=np.full_like(weighted_sum, np.nan, dtype=float),
        where=denominator != 0.0,
    )
    df = pd.DataFrame({"combine_z": combined_z}, index=values.index)

    # merge_df = pd.concat([matrix, pd.Series(weighted_matrix, name="combined")], axis=1)
    # stack_df = pd.melt(
    #     merge_df,
    #     id_vars=["ensembl_gene_id"],
    #     # Get all 'data' columns
    #     value_vars=[col for col in merge_df.columns if col not in GeneIdentifier._member_map_],
    #     var_name="source",
    #     value_name="zscore",
    # )
    # if len(stack_df["source"].unique()) > 10:
    #     stack_df = stack_df[stack_df["source"] == "combined"]
    # graph_zscore_distribution(
    #     stack_df,
    #     title=f"Combined Z-score Distribution for {context_name} - batch #{batch.batch_num}",
    #     output_filepath=output_figure_dirpath
    #     / f"{context_name}_{source.value}_batch{batch.batch_num}_combined_zscore_distribution.pdf",
    # )

    df.columns = [batch.batch_num]
    df.index = df.index.astype(int)
    df = df.sort_index()
    df.to_csv(output_combined_matrix_filepath, index=True)
    return df


def _combine_z_distribution_for_source(
    merged_source_data: pd.DataFrame,
    context_name: str,
    replicates_per_batch: Sequence[int],
    output_combined_matrix_filepath: Path,
    output_figure_filepath: Path,
    weighted_z_floor: int = -6,
    weighted_z_ceiling: int = 6,
) -> pd.DataFrame:
    """Combine z-score distributions across batches for a single source.

    Args:
        merged_source_data : DataFrame with columns: ["ensembl_gene_id", <batch1>, <batch2>, ...].
            Each batch column contains the within-batch combined z for that gene.
        context_name: tissue or condition name
        replicates_per_batch: vector of replicate counts per batch, aligned exactly to the order of
            batch columns in `merged_source_data.iloc[:, 1:]`.
        output_combined_matrix_filepath: directory to write combined z-score figures
        output_figure_filepath: filepath to write z-score figure
        weighted_z_floor: lower bound to clip z-scores
        weighted_z_ceiling : upper bound to clip z-scores

    Returns:
          A pandas dataframe of the weighted z-distributions
    """
    # If only one batch column exists, return as-is (rename like R path-through).
    if num_columns(merged_source_data) <= 1:
        logger.warning("A single batch exists for this source; returning matrix as-is.")
        out_df = merged_source_data.copy()
        out_df.columns = ["combine_z"]
        return out_df

    output_combined_matrix_filepath.parent.mkdir(parents=True, exist_ok=True)
    output_figure_filepath.parent.mkdir(parents=True, exist_ok=True)

    # alidate alignment between columns and replicate vector
    batch_column_names: list[str] = list(merged_source_data.columns)
    expected_num_batches = len(batch_column_names)
    if expected_num_batches != len(replicates_per_batch):
        raise ValueError(
            f"[{context_name}] Mismatch between number of batch columns "
            f"({expected_num_batches}: {batch_column_names}) and length of "
            f"replicates_per_batch ({len(replicates_per_batch)}: {list(replicates_per_batch)})."
        )

    logger.trace(
        f"[{context_name}] Combining {expected_num_batches} batches with"
        f" replicate counts: {dict(zip(batch_column_names, replicates_per_batch, strict=True))}"
    )

    # extract values and build masks (shape is in gene x batch)
    gene_by_batch_values: np.ndarray = merged_source_data.to_numpy(dtype=float)
    present_mask: np.ndarray = ~np.isnan(gene_by_batch_values)

    # per-batch weights: use equal weights (Stouffer's Z method)
    # Replicate counts are not statistically valid weights for z-scores because
    #   z-scores are already standardized (mean=0, std=1 under null hypothesis)
    #   Weighting by sample size assumes information about precision that was lost during standardization
    replicates_per_batch_array: np.ndarray = np.asarray(replicates_per_batch, dtype=float)
    base_weights_per_batch = np.ones_like(replicates_per_batch_array, dtype=float)  # shape: (B,)

    # drop NA columns per gene by zeroing their weights
    masked_weights_per_batch: np.ndarray = np.where(present_mask, base_weights_per_batch, 0.0)  # (G,B)
    masked_values: np.ndarray = np.where(present_mask, gene_by_batch_values, 0.0)  # (G,B)

    # weighted z combine per gene: sum(w*x)/sqrt(sum(w^2)), with NA rows → NaN
    numerator_per_gene: np.ndarray = np.sum(masked_weights_per_batch * masked_values, axis=1)  # (G,)
    denominator_per_gene: np.ndarray = np.sqrt(np.sum(masked_weights_per_batch**2, axis=1))  # (G,)
    combined_z_per_gene: np.ndarray = np.divide(
        numerator_per_gene,
        denominator_per_gene,
        out=np.full_like(numerator_per_gene, np.nan, dtype=float),
        where=denominator_per_gene != 0.0,
    )

    combined_z_per_gene = np.clip(combined_z_per_gene, weighted_z_floor, weighted_z_ceiling)

    logger.trace(f"[{context_name}] Finished combining batch z-distributions for source.")

    df = pd.DataFrame({"combine_z": combined_z_per_gene}, index=merged_source_data.index)
    df.index = df.index.astype(int)
    df = df.sort_index()
    return df


def _combine_z_distribution_for_context(
    context: str,
    zscore_results: list[_CombineOmicsInput],
    output_graph_filepath: Path,
) -> pd.DataFrame:
    if not zscore_results:
        logger.warning("No zscore results exist, returning empty dataframe")
        return pd.DataFrame({"ensembl_gene_id": [], "combine_z": []})

    z_matrices: list[pd.DataFrame] = []
    for res in zscore_results:
        matrix = res.z_score_matrix.copy()
        if len(matrix.columns) > 1:
            raise ValueError(
                f"Expected a single column for combined z-score dataframe for data '{res.type.value.lower()}'. "
                f"Got '{len(matrix.columns)}' columns"
            )

        matrix.columns = [res.type.value.lower()]
        matrix = matrix.loc[matrix.index != "-"]

        # if the matrix has duplicate ids (i.e., multiple entrez ids map to a single ensembl id, etc.)
        #   take the max value of them to remove the duplicates
        if not matrix.index.is_unique:
            matrix = matrix.groupby(level=0).max()
        z_matrices.append(matrix)

    z_matrix = pd.concat(z_matrices, axis="columns", join="outer", ignore_index=False)
    if num_columns(z_matrix) <= 1:
        logger.trace(
            f"Only 1 source exists for '{context}', returning dataframe as-is becuase no data exists to combine"
        )
        z_matrix.columns = ["combine_z"]
        return z_matrix

    values = z_matrix.values
    weights = np.array([r.weight for r in zscore_results])
    mask = ~np.isnan(values)
    masked_values = np.where(mask, values, 0)
    masked_weights = np.where(mask, weights, 0)

    numerator = np.sum(masked_weights * masked_values, axis=1)
    denominator = np.sqrt(np.sum(masked_weights**2, axis=1))
    combined_z_matrix = np.divide(
        numerator,
        denominator,
        out=np.full_like(numerator, np.nan, dtype=float),
        where=denominator != 0.0,
    )
    combined_z_matrix_df = pd.DataFrame(
        {
            "entrez_gene_id": z_matrix.index,
            "combine_z": combined_z_matrix,
        }
    )

    # combined_df = pd.DataFrame(
    #     {
    #         "ensembl_gene_id": z_matrix.index,
    #         "zscore": combined_z_matrix,
    #         "source": "combined",
    #     }
    # )
    # stack_df = pd.melt(
    #     z_matrix,
    #     id_vars=["ensembl_gene_id"],
    #     value_vars=z_matrix.columns[1:],
    #     var_name="source",
    #     value_name="zscore",
    # )
    # stack_df = pd.concat([stack_df, combined_df])
    # graph_zscore_distribution(
    #     df=stack_df,
    #     title=f"Combined Z-score Distribution for {context}",
    #     output_filepath=output_graph_filepath,
    # )
    combined_z_matrix_df["entrez_gene_id"] = combined_z_matrix_df["entrez_gene_id"].astype(int)
    combined_z_matrix_df = combined_z_matrix_df.sort_values(by="entrez_gene_id")
    return combined_z_matrix_df


def _begin_combining_distributions(
    context_name: str,
    taxon: int,
    input_matrices: _InputMatrices,
    batch_names: _BatchNames,
    source_weights: _SourceWeights,
    output_filepaths: _OutputCombinedSourceFilepath,
    output_figure_dirpath: Path,
    output_final_model_scores: Path,
    weighted_z_floor: int = -6,
    weighted_z_ceiling: int = 6,
):
    logger.info(f"Starting to combine z-scores for context '{context_name}'")
    output_figure_dirpath.mkdir(parents=True, exist_ok=True)

    z_score_results: list[_CombineOmicsInput] = []

    matrix: pd.DataFrame | None
    for source, matrix in input_matrices:
        if matrix is None:
            logger.trace(f"Source '{source.value}' is None, skipping")
            continue
        if source not in SourceTypes:
            logger.critical(f"Invalid source; got '{source.value}', expected 'trna', 'mrna', 'scrna', or 'proteomics'.")
            raise ValueError("Invalid source")

        batch_results: list[pd.DataFrame] = []
        for batch in batch_names[source.value]:
            batch: _BatchEntry
            if isinstance(matrix, pd.DataFrame):
                matrix_subset = matrix[[GeneIdentifier.entrez_gene_id.value, *batch.sample_names]]
                matrix_subset = matrix_subset.set_index(keys=[GeneIdentifier.entrez_gene_id.value], drop=True)
                matrix_subset = matrix_subset.drop(columns=["gene_symbol", "ensembl_gene_id"], errors="ignore")
            elif isinstance(matrix, sc.AnnData):
                conversion = get_remaining_identifiers(ids=matrix.var_names.tolist(), taxon=taxon)
                conversion.reset_index(drop=False, inplace=True)
                matrix = matrix.to_df().T
                matrix.reset_index(inplace=True, drop=False, names=["gene_symbol"])
                matrix = matrix.merge(conversion, left_on="gene_symbol", right_on="gene_symbol", how="left")
                matrix_subset: pd.DataFrame = matrix[[GeneIdentifier.entrez_gene_id.value, *batch.sample_names]]
                matrix_subset: pd.DataFrame = matrix_subset.set_index(
                    keys=[GeneIdentifier.entrez_gene_id.value], drop=True
                )
            else:
                raise TypeError(
                    f"Unexpected matrix type for source '{source.value}'; "
                    f"expected 'pandas.DataFrame' or 'scanpy.AnnData': {type(matrix)}"
                )

            output_fp = (
                output_filepaths[source.value].parent
                / f"{source.value}_batch{batch.batch_num}_combined_z_distrobution_{context_name}.csv"
            )
            batch_results.append(
                _combine_z_distribution_for_batch(
                    context_name=context_name,
                    batch=batch,
                    matrix=matrix_subset,
                    source=source,
                    output_combined_matrix_filepath=output_fp,
                    output_figure_dirpath=output_figure_dirpath,
                )
            )

        index_name: str = (
            "ensembl_gene_id"
            if all(df.index.name == "ensembl_gene_id" for df in batch_results)
            else "entrez_gene_id"
            if all(df.index.name == "entrez_gene_id" for df in batch_results)
            else "gene_symbol"
            if all(df.index.name == "gene_symbol" for df in batch_results)
            else ""
        )
        if not index_name:
            raise ValueError(
                f"Unable to find common gene identifier across batches for source "
                f"'{source.value}' in context '{context_name}'"
            )
        merged_batch_results = pd.concat(batch_results, axis="columns")
        merged_batch_results.index.name = index_name

        replicates_by_batch_num_map: dict[str, int] = {
            str(batch.batch_num): batch.num_samples for batch in batch_names[source.value]
        }
        batch_column_names: list[str] = list(merged_batch_results.columns)
        replicates_per_batch: list[int] = []
        for col in batch_column_names:
            key = str(col)
            if key not in replicates_by_batch_num_map:
                raise KeyError(
                    f"Missing replicate count for batch column '{col}' in source '{source.value}'. "
                    f"Known: {list(replicates_by_batch_num_map.keys())}"
                )
            replicates_per_batch.append(replicates_by_batch_num_map[key])

        merged_source_results: pd.DataFrame = _combine_z_distribution_for_source(
            merged_source_data=merged_batch_results,
            context_name=context_name,
            replicates_per_batch=replicates_per_batch,
            output_combined_matrix_filepath=(
                output_filepaths[source.value].parent
                / f"{context_name}_{source.value}_combined_zscore_distribution.csv"
            ),
            output_figure_filepath=(
                output_figure_dirpath / f"{context_name}_{source.value}_combined_zscore_distribution.pdf"
            ),
            weighted_z_floor=weighted_z_floor,
            weighted_z_ceiling=weighted_z_ceiling,
        )
        z_score_results.append(
            _CombineOmicsInput(
                z_score_matrix=merged_source_results,
                type=source,
                weight=source_weights[source.value],
            )
        )
        merged_source_results.to_csv(output_filepaths[source.value], index=True)
        logger.success(
            f"Wrote z-scores for source '{source.value}' in context '{context_name}' "
            f"to '{output_filepaths[source.value]}'"
        )

    logger.trace(f"Combining z-score distributions for all sources in context '{context_name}'")
    merged_context_results = _combine_z_distribution_for_context(
        context=context_name,
        zscore_results=z_score_results,
        output_graph_filepath=output_figure_dirpath / f"{context_name}_combined_omics_distribution.pdf",
    )
    merged_context_results.to_csv(output_final_model_scores, index=len(merged_context_results.columns) <= 1)
    logger.success(f"Wrote combined z-scores for context '{context_name}' to {output_final_model_scores}")
