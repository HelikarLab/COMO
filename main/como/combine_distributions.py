from __future__ import annotations

import asyncio
from pathlib import Path

import aiofiles.os
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
from como.graph import z_score_distribution as graph_zscore_distribution
from como.utils import (
    _num_columns,
)


async def _combine_z_distribution_for_batch(
    context_name: str,
    batch: _BatchEntry,
    matrix: pd.DataFrame,
    source: SourceTypes,
    output_combined_matrix_filepath: Path,
    output_figure_dirpath: Path,
    weighted_z_floor: int,
    weighted_z_ceiling: int,
) -> pd.DataFrame:
    await asyncio.gather(
        *[
            aiofiles.os.makedirs(output_combined_matrix_filepath.parent.as_posix(), exist_ok=True),
            aiofiles.os.makedirs(output_figure_dirpath.as_posix(), exist_ok=True),
        ]
    )

    logger.trace(
        f"Combining z-score distributions: batch #{batch.batch_num}, "
        f"samples: {len(batch.sample_names)}, "
        f"source: '{source.value}', "
        f"context: '{context_name}'"
    )
    if _num_columns(matrix) < 2:
        logger.trace(
            f"A single sample exists for batch '{batch.batch_num}'. "
            f"Returning as-is because no additional combining can be done"
        )
        return matrix

    values = matrix.iloc[:, 1:].values
    weighted_matrix = np.sum(values, axis=1) / np.sqrt(values.shape[1])
    weighted_matrix = np.clip(weighted_matrix, weighted_z_floor, weighted_z_ceiling).astype(np.int8)

    merge_df = pd.concat([matrix, pd.Series(weighted_matrix, name="combined")], axis=1)
    weighted_matrix = pd.DataFrame(
        {
            "ensembl_gene_id": matrix["ensembl_gene_id"],
            "combine_z": weighted_matrix,
        },
    )
    stack_df = pd.melt(
        merge_df,
        id_vars=["ensembl_gene_id"],
        # Get all columns except ensembl_gene_id
        value_vars=[col for col in merge_df.columns if col not in GeneIdentifier._member_map_],
        var_name="source",
        value_name="zscore",
    )
    if len(stack_df["source"].unique()) > 10:
        stack_df = stack_df[stack_df["source"] == "combined"]

    graph_zscore_distribution(
        stack_df,
        title=f"Combined Z-score Distribution for {context_name} - batch #{batch.batch_num}",
        output_filepath=output_figure_dirpath
        / f"{context_name}_{source.value}_batch{batch.batch_num}_combined_zscore_distribution.pdf",
    )

    weighted_matrix.columns = ["ensembl_gene_id", batch.batch_num]
    weighted_matrix.to_csv(output_combined_matrix_filepath, index=False)
    return weighted_matrix


async def _combine_z_distribution_for_source(
    merged_source_data: pd.DataFrame,
    context_name: str,
    num_replicates: int,
    output_combined_matrix_filepath: Path,
    output_figure_filepath: Path,
    weighted_z_floor: int = -6,
    weighted_z_ceiling: int = 6,
):
    if _num_columns(merged_source_data) <= 2:
        logger.warning("A single source exists, returning matrix as-is because no additional combining can be done")
        merged_source_data.columns = ["ensembl_gene_id", "combine_z"]
        return merged_source_data

    print(f"Making directory for '{output_combined_matrix_filepath.parent}'")
    print(f"Making directory for '{output_figure_filepath.parent}'")
    await asyncio.gather(
        *[
            aiofiles.os.makedirs(output_combined_matrix_filepath.parent.as_posix(), exist_ok=True),
            aiofiles.os.makedirs(output_figure_filepath.parent.as_posix(), exist_ok=True),
        ]
    )

    logger.trace(f"Found {_num_columns(merged_source_data) - 1} samples for context '{context_name}' to combine")
    values = merged_source_data.iloc[:, 1:].values
    mask = ~np.isnan(values)
    masked_values = np.where(mask, values, 0)  # Replace NaN with 0
    masked_num_replicates = np.where(mask, num_replicates, 0)

    weights = masked_num_replicates / np.sum(masked_num_replicates, axis=1, keepdims=True)
    numerator = np.sum(weights * masked_values, axis=1)
    denominator = np.sqrt(np.sum(weights**2, axis=1))
    weighted_matrix = numerator / denominator
    weighted_matrix = np.clip(weighted_matrix, weighted_z_floor, weighted_z_ceiling)
    logger.trace("Finished combining z-distribution")
    merge_df = pd.concat([merged_source_data, pd.Series(weighted_matrix, name="combined")], axis=1)
    weighted_matrix = pd.DataFrame(
        {
            "ensembl_gene_id": merged_source_data["ensembl_gene_id"],
            "combine_z": weighted_matrix,
        }
    )
    return weighted_matrix


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

    z_matrices = [
        res.z_score_matrix.set_index("ensembl_gene_id").rename(
            columns={col: res.type.value for col in res.z_score_matrix.columns[1:]}
        )
        for res in zscore_results
    ]
    z_matrix = pd.concat(z_matrices, axis=1, join="outer").reset_index()
    if _num_columns(z_matrix) <= 1:
        logger.trace(
            f"Only 1 source exists for '{context}', returning dataframe as-is becuase no data exists to combine"
        )
        if _num_rows(z_matrix) > 2
        else z_matrix.iloc[:, 1:].values
    ).ravel()
    merge_df = pd.concat([z_matrix, pd.Series(combined_z_matrix, name="combined")], axis=1)
    combined_z_matrix = pd.DataFrame(
    values = z_matrix.iloc[:, 1:].values
    weights = np.array([r.weight for r in zscore_results])
    mask = ~np.isnan(values)
    masked_values = np.where(mask, values, 0)
    masked_weights = np.where(mask, weights, 0)

    normalized_weights = masked_weights / np.sum(masked_weights, axis=1, keepdims=True)
    numerator = np.sum(normalized_weights * masked_values, axis=1)
    denominator = np.sqrt(np.sum(normalized_weights**2, axis=1))
    combined_z_matrix = numerator / denominator
    combined_z_matrix = np.clip(combined_z_matrix, weighted_z_floor, weighted_z_ceiling)
    combined_z_matrix_df = pd.DataFrame(
        {
            "ensembl_gene_id": z_matrix["ensembl_gene_id"],
            "combine_z": combined_z_matrix,
        }
    )

    stack_df = pd.concat(
        [
            pd.DataFrame(
                {
                    "ensembl_gene_id": merge_df["ensembl_gene_id"],
                    "zscore": merge_df[col].astype(float),
                    "source": col,
                }
            )
            for col in merge_df.columns[1:]
        ]
    )

    graph_zscore_distribution(
        df=stack_df,
        title=f"Combined Omics Z-score Distribution for {context}",
        output_png_filepath=output_png_filepath,
    )
    return combined_z_matrix


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
                        output_filepaths[source.value].parent
                        / f"{context_name}_{source.value}_batch{batch.batch_num}_combined_z_distribution_.csv"
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
            merged_batch_results = (
                df if merged_batch_results.empty else merged_batch_results.merge(df, on="ensembl_gene_id", how="outer")
            )
            zscore_results.append(_CombineOmicsInput(z_score_matrix=combine_batches_zscore, type=source, weight=weight))  # type: ignore
            combine_batches_zscore.to_csv(output_filepath, index=False)

    combined_z_omics = _combine_omics_zdistros(
        context=context_name,
        zscore_results=zscore_results,
        output_png_filepath=output_figure_dirpath / f"{context_name}_combined_omics_distribution.png",
    )
    combined_z_omics.to_csv(output_final_model_scores, index=False)
