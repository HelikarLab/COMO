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

class _IterableDataclass:
    def __iter__(self) -> Iterator[str]:
        return iter(getattr(self, field.name) for field in fields(self))  # type: ignore


class _BatchEntry(NamedTuple):
    batch_num: int
    sample_names: list[str]


@dataclass
class _InputMatrices(_IterableDataclass):
    trna: Path | pd.DataFrame | None
    mrna: Path | pd.DataFrame | None
    scrna: Path | pd.DataFrame | None
    proteomics: Path | pd.DataFrame | None


@dataclass
class _BatchNames(_IterableDataclass):
    trna: list[_BatchEntry]
    mrna: list[_BatchEntry]
    scrna: list[_BatchEntry]
    proteomics: list[_BatchEntry]


@dataclass
class _SourceWeights(_IterableDataclass):
    trna: int
    mrna: int
    scrna: int
    proteomics: int


@dataclass
class _OutputCombinedSourceFilepath(_IterableDataclass):
    trna: Path | None
    mrna: Path | None
    scrna: Path | None
    proteomics: Path | None


@dataclass
class _CombineOmicsInput:
    z_score_matrix: pd.DataFrame
    type: Literal["totalrna", "mrna", "scrna", "proteomics"]
    weight: int


def _combine_batch_zdistro(
    matrix: pd.DataFrame,
    context_name: str,
    batch_num: int,
    output_png_filepath: Path,
    weighted_z_floor: int = -6,
    weighted_z_ceiling: int = 6,
) -> pd.DataFrame:
    def weighted_z(x: npt.NDArray[float], floor: int, ceiling: int) -> npt.NDArray[float]:
        result = np.sum(x) / np.sqrt(len(x))
        return np.clip(result, floor, ceiling)

    if _num_rows(matrix) < 2:
        return matrix

    weighted_matrix = np.apply_along_axis(
        weighted_z,
        axis=1,
        arr=matrix.iloc[:, 1:].values,
        floor=weighted_z_floor,
        ceiling=weighted_z_ceiling,
    )
    merge_df = pd.concat([matrix, pd.Series(weighted_matrix, name="combined")], axis=1)
    weighted_matrix = pd.DataFrame(
        {
            "ensembl_gene_id": matrix["ensembl_gene_id"].astype(str),
            "combine_z": weighted_matrix,
        },
    )

    stack_df = pd.concat(
        [
            pd.DataFrame(
                {"ensembl_gene_id": merge_df["ensembl_gene_id"], "zscore": merge_df[col].astype(float), "source": col}
            )
            for col in merge_df.columns[1:]
        ]
    )
    if len(stack_df["source"].unique()) > 10:
        stack_df = stack_df[stack_df["source"] == "combined"]

    graph_zscore_distribution(
        stack_df,
        title=f"Combined Z-score Distribution for {context_name} - batch #{batch_num}",
        output_png_filepath=output_png_filepath,
    )
    return weighted_matrix


def _combine_context_zdistro(
    matrix: pd.DataFrame,
    context_name: str,
    batch_num: int,
    num_replicates: int,
    output_png_filepath: Path,
    weighted_z_floor: int = -6,
    weighted_z_ceiling: int = 6,
):
    def weighted_z(
        x: npt.NDArray[float],
        n_reps: int,
        floor: int,
        ceiling: int,
    ) -> npt.NDArray[float]:
        na_values = np.where(np.isnan(x))[0]
        if len(na_values) > 0:
            x = np.delete(x, na_values)
            n_reps = np.delete(n_reps, na_values)
        weights = n_reps / np.sum(n_reps)
        numerator = np.sum(weights * x)
        denominator = np.sqrt(np.sum(weights**2))
        result = numerator / denominator
        return np.clip(result, floor, ceiling)

    if _num_rows(matrix) < 2:
        matrix.columns = ["entrez_gene_id", "combine_z"]
        return matrix

    weighted_matrix = np.apply_along_axis(
        weighted_z,
        axis=1,
        arr=matrix.iloc[:, 1:].values,
        n_reps=num_replicates,
        floor=weighted_z_floor,
        ceiling=weighted_z_ceiling,
    )
    merge_df = pd.concat([matrix, pd.Series(weighted_matrix, name="combined")], axis=1)
    weighted_matrix = pd.DataFrame(
        {"ensembl_gene_id": matrix["ensembl_gene_id"].astype(str), "combine_z": weighted_matrix}
    )
    stack_df = pd.concat(
        [
            pd.DataFrame(
                {"ensembl_gene_id": merge_df["ensembl_gene_id"], "zscore": merge_df[col].astype(float), "source": col}
                for col in merge_df.columns[1:]
            )
        ]
    )
    graph_zscore_distribution(
        df=stack_df,
        title=f"Combined Z-score Distribution for {context_name} - batch #{batch_num}",
        output_png_filepath=output_png_filepath,
    )
    return weighted_matrix


def _combine_omics_zdistros(
    context: str,
    zscore_results: list[_CombineOmicsInput],
    output_png_filepath: Path,
    weighted_z_floor: int = -6,
    weighted_z_ceiling: int = 6,
):
    def weighted_z(
        x: npt.NDArray[float],
        weights,
        floor: int,
        ceiling: int,
    ):
        na_values = np.where(np.isnan(x))[0]
        if len(na_values) > 0:
            x = np.delete(x, na_values)
            weights = np.delete(weights, na_values)
        weights = weights / np.sum(weights)
        numerator = np.sum(weights * x)
        denominator = np.sqrt(np.sum(weights**2))
        result = numerator / denominator
        return np.clip(result, floor, ceiling)

    z_matrix = pd.DataFrame()
    for result in zscore_results:
        result.z_score_matrix.columns = ["ensembl_gene_id", result.type]
        z_matrix = (
            result.z_score_matrix
            if z_matrix.empty
            else pd.merge(z_matrix, result.z_score_matrix, on="ensembl_gene_id", how="outer")
        )

    combined_z_matrix = (
        np.apply_along_axis(
            weighted_z,
            axis=1,
            arr=z_matrix.iloc[:, 1:].values,
            weights=[r.weight for r in zscore_results],
            floor=weighted_z_floor,
            ceiling=weighted_z_ceiling,
        )
        if _num_rows(z_matrix) > 2
        else z_matrix.iloc[:, 1:].values
    ).ravel()
    merge_df = pd.concat([z_matrix, pd.Series(combined_z_matrix, name="combined")], axis=1)
    combined_z_matrix = pd.DataFrame(
        {
            "ensembl_gene_id": z_matrix["ensembl_gene_id"].astype(str),
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


def _combine_zscores(
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
    output_figure_dirpath.mkdir(parents=True, exist_ok=True)
    source_name: list[str] = ["totalrna", "mrna", "scrna", "proteomics"]
    zscore_results: list[_CombineOmicsInput] = []
    for matrix, source in zip(input_matrices, source_name):
        matrix: pd.DataFrame | None = pd.read_csv(matrix) if isinstance(matrix, Path) else matrix

        if matrix is not None:
            if source == "totalrna":
                batch_data = batch_names.trna
                weight = source_weights.trna
                output_filepath = output_filepaths.trna
            elif source == "mrna":
                batch_data = batch_names.mrna
                weight = source_weights.mrna
                output_filepath = output_filepaths.mrna
            elif source == "scrna":
                batch_data = batch_names.scrna
                weight = source_weights.scrna
                output_filepath = output_filepaths.scrna
            elif source == "proteomics":
                batch_data = batch_names.proteomics
                weight = source_weights.proteomics
                output_filepath = output_filepaths.proteomics
            else:
                raise ValueError(f"Invalid source; got '{source}', expected one of '{','.join(source_name)}'")

            replicate_count: list[int] = []
            merge_z_data = pd.DataFrame()

            batch: _BatchEntry
            for batch in batch_data:
                replicate_count.append(len(batch.sample_names))

                batch_df: pd.DataFrame = matrix[["ensembl_gene_id", *batch.sample_names]]
                # graph.z_score_distribution(
                #     batch_df,
                #     title=f"Z-Score Distribution for {context_name} - batch #{batch.batch_num} - {source}",
                #     output_png_filepath=output_figure_dirpath
                #     / f"{source}_batch{batch.batch_num}_zscore_distribution.png",
                # )
                combine_z_matrix: pd.DataFrame = _combine_batch_zdistro(
                    matrix=batch_df,
                    context_name=context_name,
                    batch_num=batch.batch_num,
                    output_png_filepath=(
                        output_figure_dirpath
                        / f"combined_{source}_{context_name}_batch{batch.batch_num}_distribution.png"
                    ),
                    weighted_z_floor=weighted_z_floor,
                    weighted_z_ceiling=weighted_z_ceiling,
                )
                combine_z_matrix.columns = ["ensembl_gene_id", batch.batch_num]
                merge_z_data = (
                    combine_z_matrix
                    if merge_z_data.empty
                    else pd.merge(merge_z_data, combine_z_matrix, on="ensembl_gene_id", how="outer")
                )
            combine_batches_zscore = _combine_context_zdistro(
                matrix=merge_z_data,
                context_name=context_name,
                batch_num=batch.batch_num,
                num_replicates=sum(replicate_count),
                output_png_filepath=output_figure_dirpath / f"totalrna_{context_name}_combined_distribution.png",
                weighted_z_floor=weighted_z_floor,
                weighted_z_ceiling=weighted_z_ceiling,
            )
            zscore_results.append(_CombineOmicsInput(z_score_matrix=combine_batches_zscore, type=source, weight=weight))  # type: ignore
            combine_batches_zscore.to_csv(output_filepath, index=False)

    combined_z_omics = _combine_omics_zdistros(
        context=context_name,
        zscore_results=zscore_results,
        output_png_filepath=output_figure_dirpath / f"{context_name}_combined_omics_distribution.png",
    )
    combined_z_omics.to_csv(output_final_model_scores, index=False)
