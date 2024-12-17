from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass, fields
from pathlib import Path
from typing import Literal, NamedTuple

import numpy as np
import numpy.typing as npt
import pandas as pd

from como import graph
from como.utils import _num_rows


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



def _combine_batch_zdistro(wd, context, batch, zmat):
    plot_name_png = wd / "figures" / f"plot_{context}_{batch}_combine_distro.png"

    def weighted_z(x):
        floor_score = -6
        ceil_score = 6
        x = np.array(x, dtype=float)
        result = np.sum(x) / np.sqrt(len(x))
        return np.clip(result, floor_score, ceil_score)

    combine_z = zmat
    if zmat.shape[1] > 2:
        combine_z = np.apply_along_axis(weighted_z, axis=1, arr=zmat.iloc[:, 1:].values)
        merge_df = pd.concat([zmat, pd.Series(combine_z, name="combined")], axis=1)
        combine_z = pd.DataFrame({"entrez_gene_id": zmat["entrez_gene_id"].astype(str), "combine_z": combine_z})

        stack_df = pd.concat(
            [
                pd.DataFrame(
                    {"entrez_gene_id": merge_df["entrez_gene_id"], "zscore": merge_df[col].astype(float), "source": col}
                )
                for col in merge_df.columns[1:]
            ]
        )

        # Simplified plot for many sources (optional)
        if len(stack_df["source"].unique()) > 10:
            stack_df = stack_df[stack_df["source"] == "combined"]

        fig = px.histogram(
            stack_df,
            x="zscore",
            color="source",
            nbins=100,
            marginal="rug",
            title=f"Combined Z-score Distribution for {context} - {batch}",
        )
        fig.update_layout(xaxis_title="Z-score", yaxis_title="Frequency", font={"family": "sans-serif", "size": 12})
        fig.write_image(plot_name_png)

    return combine_z


def _combine_context_zdistro(wd, context, n_reps, zmat):
    plot_name_png = wd / "figures" / f"plot_{context}_combine_batches_distro.png"

    def weighted_z(x, n_reps):
        floor_score = -6
        ceil_score = 6
        x = np.array(x, dtype=float)
        nas = np.where(np.isnan(x))[0]
        if len(nas) > 0:
            x = np.delete(x, nas)
            n_reps = np.delete(n_reps, nas)
        weights = n_reps / np.sum(n_reps)
        numer = np.sum(weights * x)
        denom = np.sqrt(np.sum(weights**2))
        result = numer / denom
        return np.clip(result, floor_score, ceil_score)

    if zmat.shape[1] > 2:
        combine_z = np.apply_along_axis(weighted_z, axis=1, arr=zmat.iloc[:, 1:].values, n_reps=n_reps)
        merge_df = pd.concat([zmat, pd.Series(combine_z, name="combined")], axis=1)
        combine_z = pd.DataFrame({"entrez_gene_id": zmat["entrez_gene_id"].astype(str), "combine_z": combine_z})

        stack_df = pd.concat(
            [
                pd.DataFrame(
                    {"entrez_gene_id": merge_df["entrez_gene_id"], "zscore": merge_df[col].astype(float), "source": col}
                )
                for col in merge_df.columns[1:]
            ]
        )

        fig = px.histogram(
            stack_df,
            x="zscore",
            color="source",
            nbins=100,  # Adjust as needed
            marginal="rug",
            title=f"Combined Batches Z-score Distribution for {context}",
        )

        fig.update_layout(xaxis_title="Z-score", yaxis_title="Frequency", font=dict(family="sans-serif", size=12))  # noqa: C408

        fig.write_image(plot_name_png)

    else:
        combine_z = zmat
        combine_z.columns = ["entrez_gene_id", "combine_z"]

    return combine_z


def _combine_omics_zdistros(
    wd,
    context,
    comb_batches_z_trna,
    comb_batches_z_mrna,
    comb_batches_z_scrna,
    comb_batches_z_prot,
    trna_weight: int,
    mrna_weight: int,
    scrna_weight: int,
    proteomics_weight: int,
    keep_gene_scores=True,
):
    fig_path = wd / context / "figures"
    if not fig_path.exists():
        fig_path.mkdir(parents=True)
    plot_name_png = fig_path / f"plot_{context}_combine_omics_distro.png"

    weights = []
    names = []
    dfs = []
    counter = 0
    if comb_batches_z_trna is not None:
        counter += 1
        weights.append(trna_weight)
        names.append("total")
        dfs.append(comb_batches_z_trna)
    if comb_batches_z_mrna is not None:
        counter += 1
        weights.append(mrna_weight)
        names.append("polyA")
        dfs.append(comb_batches_z_mrna)
    if comb_batches_z_scrna is not None:
        counter += 1
        weights.append(scrna_weight)
        names.append("singleCell")
        dfs.append(comb_batches_z_scrna)
    if comb_batches_z_prot is not None:
        counter += 1
        weights.append(proteomics_weight)
        names.append("proteome")
        dfs.append(comb_batches_z_prot)

    def weighted_z(x, weights):
        floor_score = -6
        ceil_score = 10
        x = np.array(x, dtype=float)
        nas = np.where(np.isnan(x))[0]
        if len(nas) > 0:
            x = np.delete(x, nas)
            weights = np.delete(weights, nas)
        weights = weights / np.sum(weights)
        numer = np.sum(weights * x)
        denom = np.sqrt(np.sum(weights**2))
        result = numer / denom
        return np.clip(result, floor_score, ceil_score)

    zmat = dfs[0].copy()
    zmat.columns = ["entrez_gene_id", names[0]]
    for i in range(1, counter):
        add_df = dfs[i]
        add_df.columns = ["entrez_gene_id", names[i]]
        zmat = pd.merge(zmat, add_df, on="entrez_gene_id", how="outer")

    if zmat.shape[1] > 2:
        combine_z = np.apply_along_axis(weighted_z, axis=1, arr=zmat.iloc[:, 1:].values, weights=weights)
    else:
        combine_z = zmat.iloc[:, 1:].values

    combine_z = combine_z.ravel()
    merge_df = pd.concat([zmat, pd.Series(combine_z, name="combined")], axis=1)
    combine_z = pd.DataFrame({"entrez_gene_id": zmat["entrez_gene_id"].astype(str), "combine_z": combine_z})

    stack_df = pd.concat(
        [
            pd.DataFrame(
                {"entrez_gene_id": merge_df["entrez_gene_id"], "zscore": merge_df[col].astype(float), "source": col}
            )
            for col in merge_df.columns[1:]
        ]
    )

    fig = px.histogram(
        stack_df,
        x="zscore",
        color="source",
        nbins=100,  # Adjust as needed
        marginal="rug",
        title=f"Combined Omics Z-score Distribution for {context}",
    )

    fig.update_layout(xaxis_title="Z-score", yaxis_title="Frequency", font={"family": "sans-serif", "size": 12})

    fig.write_image(plot_name_png)

    return combine_z


def _combine_zscores(  # noqa: C901
    working_dir,
    context_names,
    global_use_mrna,
    global_use_trna,
    global_use_scrna,
    global_use_proteins,
    keep_gene_scores,
    global_trna_weight,
    global_mrna_weight,
    global_scrna_weight,
    global_protein_weight,
):
    working_dir = Path(working_dir)
    figure_output_dir = working_dir / "figures"
    figure_output_dir.mkdir(parents=True, exist_ok=True)

    global_trna_batches = _parse_contexts_zfpkm(working_dir, context_names, "total")
    global_mrna_batches = _parse_contexts_zfpkm(working_dir, context_names, "mrna")
    global_scrna_batches = _parse_contexts_zumi(working_dir, context_names, "scrna")
    global_protein_batches = _parse_contexts_zscore_prot(working_dir, context_names)

    for context in context_names:
        context_use_trna = global_use_trna
        context_use_mrna = global_use_mrna
        context_use_scrna = global_use_scrna
        context_use_proteins = global_use_proteins

        context_trna_batch = global_trna_batches.get(context, [])
        if len(context_trna_batch) == 0 and global_use_trna:
            context_use_trna = False
            logger.warning(f"No total RNA-seq zFPKM Matrix files found for {context}. Will not use for this context.")

        context_mrna_batch = global_mrna_batches.get(context, [])
        if len(context_mrna_batch) == 0 and global_use_mrna:
            context_use_mrna = False
            logger.warning(f"No polyA RNA-seq zFPKM Matrix files found for {context}. Will not use for this context.")

        context_scrna_batch = global_scrna_batches.get(context, [])
        if len(context_scrna_batch) == 0 and global_use_scrna:
            context_use_scrna = False
            logger.warning(f"No SC RNA-seq zFPKM Matrix files found for {context}. Will not use for this context.")

        context_protein_batch = global_protein_batches.get(context, [])
        if len(context_protein_batch) == 0 and global_use_proteins:
            context_use_proteins = False
            logger.warning(f"No proteomics z-score Matrix files found for {context}. Will not use for this context.")

        comb_batches_z_trna = None
        if context_use_trna:
            logger.info("Will merge total RNA-seq distributions")
            trna_workdir = working_dir / context / "total"
            num_reps = []
            merge_z_data = pd.DataFrame()

            for batch in context_trna_batch:
                res = _merge_batch(trna_workdir, context, batch)
                z_matrix = res[0]
                num_reps.extend(res[1])
                combine_z_matrix = _combine_batch_zdistro(trna_workdir, context, batch, z_matrix)
                combine_z_matrix.columns = ["entrez_gene_id", batch]
                merge_z_data = (
                    combine_z_matrix
                    if merge_z_data.empty
                    else pd.merge(merge_z_data, combine_z_matrix, on="entrez_gene_id", how="outer")
                )

            comb_batches_z_trna = _combine_context_zdistro(trna_workdir, context, num_reps, merge_z_data)
            filename = trna_workdir / f"combined_zFPKM_{context}.csv"
            comb_batches_z_trna.to_csv(filename, index=False)

            if not context_use_proteins and not context_use_mrna and not context_use_scrna:
                filename = trna_workdir / f"model_scores_{context}.csv"
                comb_batches_z_trna.to_csv(filename, index=False)

        comb_batches_z_mrna = None
        if context_use_mrna:
            logger.info("Will merge polyA enriched RNA-seq distributions")
            mrna_workdir = working_dir / context / "mrna"
            num_reps = []
            merge_z_data = pd.DataFrame()
            for batch in context_mrna_batch:
                res = _merge_batch(mrna_workdir, context, batch)
                z_matrix = res[0]
                num_reps.extend(res[1])
                combine_z_matrix = _combine_batch_zdistro(mrna_workdir, context, batch, z_matrix)
                combine_z_matrix.columns = ["entrez_gene_id", batch]
                merge_z_data = (
                    combine_z_matrix
                    if merge_z_data.empty
                    else pd.merge(merge_z_data, combine_z_matrix, on="entrez_gene_id", how="outer")
                )

            comb_batches_z_mrna = _combine_context_zdistro(mrna_workdir, context, num_reps, merge_z_data)
            filename = mrna_workdir / f"combined_zFPKM_{context}.csv"
            comb_batches_z_mrna.to_csv(filename, index=False)

            if not context_use_proteins and not context_use_trna and not context_use_scrna:
                filename = mrna_workdir / f"model_scores_{context}.csv"
                comb_batches_z_mrna.to_csv(filename, index=False)

        comb_batches_z_scrna = None
        if context_use_scrna:
            logger.info(f"Will merge single-cell RNA-seq distributions for {context}")
            scrna_workdir = working_dir / context / "scrna"
            num_reps = []
            merge_z_data = pd.DataFrame()
            for batch in context_scrna_batch:
                res = _merge_batch(scrna_workdir, context, batch)
                z_matrix = res[0]
                num_reps.extend(res[1])
                combine_z_matrix = _combine_batch_zdistro(scrna_workdir, context, batch, z_matrix)
                combine_z_matrix.columns = ["entrez_gene_id", batch]
                merge_z_data = (
                    combine_z_matrix
                    if merge_z_data.empty
                    else pd.merge(merge_z_data, combine_z_matrix, on="entrez_gene_id", how="outer")
                )

            comb_batches_z_scrna = _combine_context_zdistro(scrna_workdir, context, num_reps, merge_z_data)
            filename = scrna_workdir / f"combined_zFPKM_{context}.csv"
            comb_batches_z_scrna.to_csv(filename, index=False)

            if not context_use_proteins and not context_use_trna and not context_use_mrna:
                filename = scrna_workdir / f"model_scores_{context}.csv"
                comb_batches_z_scrna.to_csv(filename, index=False)

        comb_batches_z_prot = None
        if context_use_proteins:
            logger.info("Will merge protein abundance distributions")
            protein_workdir = working_dir / context / "proteomics"
            num_reps = []
            merge_z_data = pd.DataFrame()
            for batch in context_protein_batch:
                res = _merge_batch(protein_workdir, context, batch)
                z_matrix = res[0]
                num_reps.extend(res[1])
                combine_z_matrix = _combine_batch_zdistro(protein_workdir, context, batch, z_matrix)
                combine_z_matrix.columns = ["entrez_gene_id", batch]
                merge_z_data = (
                    combine_z_matrix
                    if merge_z_data.empty
                    else pd.merge(merge_z_data, combine_z_matrix, on="entrez_gene_id", how="outer")
                )

            comb_batches_z_prot = _combine_context_zdistro(protein_workdir, context, num_reps, merge_z_data)
            filename = protein_workdir / f"combined_zscore_proteinAbundance_{context}.csv"
            comb_batches_z_prot.to_csv(filename, index=False)

            if not context_use_mrna and not context_use_trna and not context_use_scrna:
                filename = protein_workdir / f"model_scores_{context}.csv"
                comb_batches_z_prot.to_csv(filename, index=False)

        if (
            comb_batches_z_mrna is None
            and comb_batches_z_trna is None
            and comb_batches_z_scrna is None
            and comb_batches_z_prot is None
        ):
            logger.critical(
                f"The context '{context}' was found in the configuration file but no data was found on disk!"
            )
            continue

        comb_omics_z = _combine_omics_zdistros(
            wd=working_dir,
            context=context,
            comb_batches_z_trna=comb_batches_z_trna,
            comb_batches_z_mrna=comb_batches_z_mrna,
            comb_batches_z_scrna=comb_batches_z_scrna,
            comb_batches_z_prot=comb_batches_z_prot,
            trna_weight=global_trna_weight,
            mrna_weight=global_mrna_weight,
            scrna_weight=global_scrna_weight,
            proteomics_weight=global_protein_weight,
        )

        filename = working_dir / context / f"model_scores_{context}.csv"
        comb_omics_z.to_csv(filename, index=False)
