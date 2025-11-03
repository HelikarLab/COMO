from __future__ import annotations

from collections.abc import Sequence
from pathlib import Path
from typing import cast

import numpy as np
import pandas as pd
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
from como.utils import LogLevel, _log_and_raise_error, _num_columns


def _combine_z_distribution_for_batch(
    context_name: str,
    batch: _BatchEntry,
    matrix: pd.DataFrame,
    source: SourceTypes,
    output_combined_matrix_filepath: Path,
    output_figure_dirpath: Path,
    weighted_z_floor: int,
    weighted_z_ceiling: int,
) -> pd.DataFrame:
    """Combine z-score distributions across samples for a single batch.

    Args:
        context_name: Name of the context (e.g., tissue or condition).
        batch: Batch entry containing batch number and sample names.
        matrix: DataFrame with 'ensembl_gene_id' and sample columns.
        source: Source type (e.g., trna, mrna, scrna, proteomics).
        output_combined_matrix_filepath: Path to save the combined z-score matrix.
        output_figure_dirpath: Path to save the z-score distribution figure.
        weighted_z_floor: Minimum z-score value after combining.
        weighted_z_ceiling: Maximum z-score value after combining.

    Returns:
            A pandas dataframe of the weighted z-distributions
    """
    output_combined_matrix_filepath.parent.mkdir(parents=True, exist_ok=True)
    output_figure_dirpath.mkdir(parents=True, exist_ok=True)

    logger.trace(
        f"Combining z-score distributions: batch #{batch.batch_num}, "
        f"samples: {len(batch.sample_names)}, "
        f"source: '{source.value}', "
        f"context: '{context_name}'"
    )
    if _num_columns(matrix) < 2:
        logger.trace(f"A single sample exists for batch '{batch.batch_num}'. Returning as-is because no additional combining can be done")
        return matrix

    values = matrix.values
    weighted_matrix = np.clip(
        a=np.sum(values, axis=1) / np.sqrt(values.shape[1]),  # calculate a weighted matrix
        a_min=weighted_z_floor,
        a_max=weighted_z_ceiling,
    ).astype(float)

    weighted_matrix = pd.DataFrame({"combine_z": weighted_matrix}, index=matrix.index)

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

    weighted_matrix.columns = [batch.batch_num]
    weighted_matrix.to_csv(output_combined_matrix_filepath, index=True)
    return weighted_matrix


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
        merged_source_data : DataFrame with columns: ["ensembl_gene_id", <batch1>, <batch2>, ...]. Each batch column contains the within-batch combined z for that gene.
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
    if _num_columns(merged_source_data) <= 1:
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
        f"[{context_name}] Combining {expected_num_batches} batches with replicate counts: {dict(zip(batch_column_names, replicates_per_batch))}"
    )

    # extract values and build masks (shape is in gene x batch)
    gene_by_batch_values: np.ndarray = merged_source_data.to_numpy(dtype=float)
    present_mask: np.ndarray = ~np.isnan(gene_by_batch_values)

    # per-batch weights from replicate counts, normalized over all batches
    replicates_per_batch_array: np.ndarray = np.asarray(replicates_per_batch, dtype=float)
    total_replicates: float = replicates_per_batch_array.sum()
    if total_replicates <= 0.0:
        raise ValueError(f"[{context_name}] Sum of replicates_per_batch must be > 0; got {total_replicates}.")
    base_weights_per_batch: np.ndarray = replicates_per_batch_array / total_replicates  # shape: (B,)

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

    return pd.DataFrame({"combine_z": combined_z_per_gene}, index=merged_source_data.index)


def _combine_z_distribution_for_context(
    context: str,
    zscore_results: list[_CombineOmicsInput],
    output_graph_filepath: Path,
    weighted_z_floor: int = -6,
    weighted_z_ceiling: int = 6,
):
    if not zscore_results:
        logger.warning("No zscore results exist, returning empty dataframe")
        return pd.DataFrame({"ensembl_gene_id": [], "combine_z": []})

    z_matrices: list[pd.DataFrame] = []
    for res in zscore_results:
        matrix = res.z_score_matrix.copy()
        if len(matrix.columns) > 1:
            _log_and_raise_error(
                f"Expected a single column for combined z-score dataframe for data '{res.type.value.lower()}'. Got '{len(matrix.columns)}' columns",
                error=ValueError,
                level=LogLevel.ERROR,
            )

        matrix.columns = [res.type.value.lower()]
        matrix = cast(pd.DataFrame, matrix.loc[matrix.index != "-"])

        # if the matrix has duplicate ids (i.e., multiple entrez ids map to a single ensembl id, etc.)
        #   take the mean value of them to remove the duplicates
        if not matrix.index.is_unique:
            # matrix = cast(pd.DataFrame, matrix.groupby(level=0).mean())
            matrix = cast(pd.DataFrame, matrix.groupby(level=0).max())
        z_matrices.append(matrix)

    z_matrix = pd.concat(z_matrices, axis="columns", join="outer", ignore_index=False)
    if _num_columns(z_matrix) <= 1:
        logger.trace(f"Only 1 source exists for '{context}', returning dataframe as-is becuase no data exists to combine")
        z_matrix.columns = ["combine_z"]
        return z_matrix

    values = z_matrix.iloc[:, 1:].values
    weights = np.array([r.weight for r in zscore_results])
    mask = ~np.isnan(values)
    masked_values = np.where(mask, values, 0)
    masked_weights = np.where(mask, weights, 0)

    normalized_weights = masked_weights / np.sum(masked_weights, axis=1, keepdims=True)
    numerator = np.sum(normalized_weights * masked_values, axis=1)
    denominator = np.sqrt(np.sum(normalized_weights**2, axis=1))
    combined_z_matrix = numerator / denominator
    combined_z_matrix = np.clip(combined_z_matrix, weighted_z_floor, weighted_z_ceiling, dtype=float)
    combined_z_matrix_df = pd.DataFrame(
        {
            "ensembl_gene_id": z_matrix.index,
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
    return combined_z_matrix_df


async def _begin_combining_distributions(
    context_name: str,
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
    for source, matrix in input_matrices:
        if matrix is None:
            logger.trace(f"Source '{source.value}' is None, skipping")
            continue
        if source not in SourceTypes:
            logger.critical(f"Invalid source; got '{source.value}', expected 'trna', 'mrna', 'scrna', or 'proteomics'.")
            raise ValueError("Invalid source")

        batch_results = await asyncio.gather(
            *[
                _combine_z_distribution_for_batch(
                    context_name=context_name,
                    batch=batch,
                    matrix=matrix[[GeneIdentifier.ENSEMBL_GENE_ID.value, *batch.sample_names]],
                    source=source,
                    output_combined_matrix_filepath=(
                        output_filepaths[source.value].parent / f"{context_name}_{source.value}_batch{batch.batch_num}_combined_z_distribution.csv"
                    ),
                    output_figure_dirpath=output_figure_dirpath,
                    weighted_z_floor=weighted_z_floor,
                    weighted_z_ceiling=weighted_z_ceiling,
                )
                for batch in batch_names[source.value]
            ]
        )

        merged_batch_results = pd.DataFrame()
        for df in batch_results:
            merged_batch_results = df if merged_batch_results.empty else merged_batch_results.merge(df, on="ensembl_gene_id", how="outer")

        merged_source_results: pd.DataFrame = _combine_z_distribution_for_source(
            merged_source_data=merged_batch_results,
            context_name=context_name,
            replicates_per_batch=replicates_per_batch,
            output_combined_matrix_filepath=(
                output_filepaths[source.value].parent / f"{context_name}_{source.value}_combined_zscore_distribution.csv"
            ),
            output_figure_filepath=(output_figure_dirpath / f"{context_name}_{source.value}_combined_zscore_distribution.pdf"),
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
        merged_source_results.to_csv(output_filepaths[source.value], index=False)
        logger.success(f"Wrote z-scores for source '{source.value}' in context '{context_name}' to '{output_filepaths[source.value]}'")

    logger.trace(f"Combining z-score distributions for all sources in context '{context_name}'")
    merged_context_results = _combine_z_distribution_for_context(
        context=context_name,
        zscore_results=z_score_results,
        output_graph_filepath=output_figure_dirpath / f"{context_name}_combined_omics_distribution.pdf",
    )
    merged_context_results.to_csv(output_final_model_scores, index=False)
    logger.success(f"Finished combining z-scores for context '{context_name}'")
