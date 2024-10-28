from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
from loguru import logger


def _get_batch_name(x):
    return Path(x).stem


def _parse_contexts_zfpkm(wd, contexts, prep):
    batches = {}
    for context in contexts:
        files = (wd / context / prep).glob(f"zFPKM_Matrix_{prep}_*.csv")
        batches[context] = [_get_batch_name(f) for f in files]
    return batches


def _parse_contexts_zumi(wd, contexts, prep):
    batches = {}
    for context in contexts:
        files = (wd / context / prep).glob(f"zUMI_Matrix_{prep}_*.csv")
        batches[context] = [_get_batch_name(f) for f in files]
    return batches


def _parse_contexts_zscore_prot(wd, contexts):
    batches = {}
    for context in contexts:
        files = (wd / context / "proteomics").glob("protein_zscore_Matrix_*.csv")
        batches[context] = [_get_batch_name(f) for f in files]
    return batches


def _merge_batch(wd, context, batch):
    files = list(wd.glob(f"*{batch}*"))
    nrep = []
    if not files:
        raise ValueError(f"No files found for {context} - {batch}")

    for f in files:
        zmat = pd.read_csv(f)
        zmat.columns = pd.Index([c.lower() for c in zmat.columns])
        zmat["entrez_gene_id"] = zmat["entrez_gene_id"].str.split("//").str[0]

        zmat = zmat.astype({col: float for col in zmat.columns if col != "entrez_gene_id"})
        zmat = zmat.astype({"entrez_gene_id": str})
        zmat = zmat.groupby("entrez_gene_id").max().reset_index()
        zmat = zmat.dropna()

        nrep.append(zmat.shape[1] - 1)
        entrez_gene = zmat["entrez_gene_id"]
        rep_names = zmat.columns
        zmat = pd.concat([zmat[col] for col in zmat.columns[1:]], axis=1)
        zmat = pd.concat([entrez_gene, zmat], axis=1)
        zmat = zmat.dropna()
        zmat.columns = rep_names

        stack_df = pd.concat(
            [pd.DataFrame({"entrez_gene_id": zmat["entrez_gene_id"], "zscore": zmat[col].astype(float), "source": col}) for col in zmat.columns[1:]]
        )

        plot_name_png = wd / "figures" / f"plot_{context}_{Path(f).stem}.png"

        fig = px.histogram(
            stack_df,
            x="zscore",
            color="source",
            nbins=100,  # Adjust as needed
            marginal="rug",
            title=f"Z-score Distribution for {context} - {Path(f).stem}",
        )

        fig.update_layout(xaxis_title="Z-score", yaxis_title="Frequency", font=dict(family="sans-serif", size=12))

        # Simplified plot for many sources (optional)
        if len(stack_df["source"].unique()) > 10:
            fig.update_layout(showlegend=False)

        fig.write_image(plot_name_png)

    return zmat, nrep


def _combine_batch_zdistro(wd, context, batch, zmat):
    plot_name_png = wd / "figures" / f"plot_{context}_{batch}_combine_distro.png"

    def weighted_z(x):
        floor_score = -6
        ceil_score = 6
        x = np.array(x, dtype=float)
        result = np.sum(x) / np.sqrt(len(x))
        return np.clip(result, floor_score, ceil_score)

    if zmat.shape[1] > 2:
        combine_z = np.apply_along_axis(weighted_z, axis=1, arr=zmat.iloc[:, 1:].values)
        merge_df = pd.concat([zmat, pd.Series(combine_z, name="combined")], axis=1)
        combine_z = pd.DataFrame({"entrez_gene_id": zmat["entrez_gene_id"].astype(str), "combine_z": combine_z})

        stack_df = pd.concat(
            [
                pd.DataFrame({"entrez_gene_id": merge_df["entrez_gene_id"], "zscore": merge_df[col].astype(float), "source": col})
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
            nbins=100,  # Adjust as needed
            marginal="rug",
            title=f"Combined Z-score Distribution for {context} - {batch}",
        )

        fig.update_layout(xaxis_title="Z-score", yaxis_title="Frequency", font=dict(family="sans-serif", size=12))

        fig.write_image(plot_name_png)

    else:
        combine_z = zmat

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
                pd.DataFrame({"entrez_gene_id": merge_df["entrez_gene_id"], "zscore": merge_df[col].astype(float), "source": col})
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

        fig.update_layout(xaxis_title="Z-score", yaxis_title="Frequency", font=dict(family="sans-serif", size=12))

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
    tweight,
    mweight,
    sweight,
    pweight,
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
    if tweight > 0:
        counter += 1
        weights.append(tweight)
        names.append("total")
        dfs.append(comb_batches_z_trna)
    if mweight > 0:
        counter += 1
        weights.append(mweight)
        names.append("polyA")
        dfs.append(comb_batches_z_mrna)
    if sweight > 0:
        counter += 1
        weights.append(sweight)
        names.append("singleCell")
        dfs.append(comb_batches_z_scrna)
    if pweight > 0:
        counter += 1
        weights.append(pweight)
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

    merge_df = pd.concat([zmat, pd.Series(combine_z, name="combined")], axis=1)
    combine_z = pd.DataFrame({"entrez_gene_id": zmat["entrez_gene_id"].astype(str), "combine_z": combine_z})

    stack_df = pd.concat(
        [
            pd.DataFrame({"entrez_gene_id": merge_df["entrez_gene_id"], "zscore": merge_df[col].astype(float), "source": col})
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

    fig.update_layout(xaxis_title="Z-score", yaxis_title="Frequency", font=dict(family="sans-serif", size=12))

    fig.write_image(plot_name_png)

    return combine_z


def combine_zscores_main(
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
    if not figure_output_dir.exists():
        figure_output_dir.mkdir()

    global_trna_batches = _parse_contexts_zfpkm(working_dir, context_names, "total")
    global_mrna_batches = _parse_contexts_zfpkm(working_dir, context_names, "mrna")
    global_scrna_batches = _parse_contexts_zumi(working_dir, context_names, "scrna")
    global_protein_batches = _parse_contexts_zscore_prot(working_dir, context_names)

    for context in context_names:
        context_use_trna = global_use_trna
        context_use_mrna = global_use_mrna
        context_use_scrna = global_use_scrna
        context_use_proteins = global_use_proteins

        context_trna_weight = global_trna_weight
        context_mrna_weight = global_mrna_weight
        context_scrna_weight = global_scrna_weight
        context_protein_weight = global_protein_weight

        context_trna_batch = global_trna_batches.get(context, [])
        context_mrna_batch = global_mrna_batches.get(context, [])
        context_scrna_batch = global_scrna_batches.get(context, [])
        context_protein_batch = global_protein_batches.get(context, [])

        if len(context_trna_batch) == 0 and global_use_trna:
            context_use_trna = False
            logger.warning(f"No total RNA-seq zFPKM Matrix files found for {context}. Will not use for this context.")

        if len(context_mrna_batch) == 0 and global_use_mrna:
            context_use_mrna = False
            logger.warning(f"No polyA RNA-seq zFPKM Matrix files found for {context}. Will not use for this context.")

        if len(context_scrna_batch) == 0 and global_use_scrna:
            context_use_scrna = False
            logger.warning(f"No SC RNA-seq zFPKM Matrix files found for {context}. Will not use for this context.")

        if len(context_protein_batch) == 0 and global_use_proteins:
            context_use_proteins = False
            logger.warning(f"No proteomics z-score Matrix files found for {context}. Will not use for this context.")

        if context_use_trna:
            logger.info("Will merge total RNA-seq distributions")
            trna_workdir = working_dir / context / "total"
            num_reps = []
            count = 0
            merge_z = pd.DataFrame()  # Initialize an empty DataFrame
            for batch in context_trna_batch:
                res = merge_batch(trna_workdir, context, batch)
                zmat = res[0]
                num_reps.extend(res[1])
                comb_z = combine_batch_zdistro(trna_workdir, context, batch, zmat)
                comb_z.columns = ["ENTREZ_GENE_ID", batch]
                if merge_z.empty:
                    merge_z = comb_z
                combine_z_matrix = _combine_batch_zdistro(trna_workdir, context, batch, z_matrix)
                combine_z_matrix.columns = ["entrez_gene_id", batch]
                if merge_z_data.empty:
                    merge_z_data = combine_z_matrix
                else:
                    merge_z_data = pd.merge(merge_z_data, combine_z_matrix, on="entrez_gene_id", how="outer")
                count += 1

            comb_batches_z_trna = _combine_context_zdistro(trna_workdir, context, num_reps, merge_z_data)
            filename = trna_workdir / f"combined_zFPKM_{context}.csv"
            comb_batches_z_trna.to_csv(filename, index=False)

            if not context_use_proteins and not context_use_mrna and not context_use_scrna:
                filename = working_dir / context / "total" / f"model_scores_{context}.csv"
                comb_batches_z_trna.to_csv(filename, index=False)

        else:
            comb_batches_z_trna = None

        if context_use_mrna:
            logger.info("Will merge polyA enriched RNA-seq distributions")
            mrna_workdir = working_dir / context / "mrna"
            num_reps = []
            count = 0
            merge_z = pd.DataFrame()  # Initialize an empty DataFrame
            for batch in context_mrna_batch:
                res = _merge_batch(mrna_workdir, context, batch)
                z_matrix = res[0]
                num_reps.extend(res[1])
                combine_z_matrix = _combine_batch_zdistro(mrna_workdir, context, batch, z_matrix)
                combine_z_matrix.columns = ["entrez_gene_id", batch]
                if merge_z_data.empty:
                    merge_z_data = combine_z_matrix
                else:
                    merge_z = pd.merge(merge_z, comb_z, on="ENTREZ_GENE_ID", how="outer")
                count += 1

            comb_batches_z_mrna = _combine_context_zdistro(mrna_workdir, context, num_reps, merge_z_data)
            filename = mrna_workdir / f"combined_zFPKM_{context}.csv"
            comb_batches_z_mrna.to_csv(filename, index=False)

            if not context_use_proteins and not context_use_trna and not context_use_scrna:
                filename = mrna_workdir / f"model_scores_{context}.csv"
                comb_batches_z_mrna.to_csv(filename, index=False)

        else:
            comb_batches_z_mrna = None

        if context_use_scrna:
            logger.info("Will merge single-cell RNA-seq distributions")
            scrna_workdir = working_dir / context / "scrna"
            num_reps = []
            count = 0
            merge_z = pd.DataFrame()  # Initialize an empty DataFrame
            for batch in context_scrna_batch:
                res = merge_batch(scrna_workdir, context, batch)
                zmat = res[0]
                num_reps.extend(res[1])
                combine_z_matrix = _combine_batch_zdistro(scrna_workdir, context, batch, z_matrix)
                combine_z_matrix.columns = ["entrez_gene_id", batch]
                if merge_z_data.empty:
                    merge_z_data = combine_z_matrix
                else:
                    merge_z_data = pd.merge(merge_z_data, combine_z_matrix, on="entrez_gene_id", how="outer")
                count += 1

            comb_batches_z_scrna = _combine_context_zdistro(scrna_workdir, context, num_reps, merge_z_data)
            filename = scrna_workdir / f"combined_zFPKM_{context}.csv"
            comb_batches_z_scrna.to_csv(filename, index=False)

            if not context_use_proteins and not context_use_trna and not context_use_mrna:
                filename = scrna_workdir / f"model_scores_{context}.csv"
                comb_batches_z_scrna.to_csv(filename, index=False)

        else:
            comb_batches_z_scrna = None

        if context_use_proteins:
            logger.info("Will merge protein abundance distributions")
            protein_workdir = working_dir / context / "proteomics"
            num_reps = []
            count = 0
            merge_z = pd.DataFrame()  # Initialize an empty DataFrame
            for batch in context_protein_batch:
                res = _merge_batch(protein_workdir, context, batch)
                z_matrix = res[0]
                num_reps.extend(res[1])
                combine_z_matrix = _combine_batch_zdistro(protein_workdir, context, batch, z_matrix)
                combine_z_matrix.columns = ["entrez_gene_id", batch]
                if merge_z_data.empty:
                    merge_z_data = combine_z_matrix
                else:
                    merge_z_data = pd.merge(merge_z_data, combine_z_matrix, on="entrez_gene_id", how="outer")
                count += 1

            comb_batches_z_prot = _combine_context_zdistro(protein_workdir, context, num_reps, merge_z_data)
            filename = protein_workdir / f"combined_zscore_proteinAbundance_{context}.csv"
            comb_batches_z_prot.to_csv(filename, index=False)

            if not context_use_mrna and not context_use_trna and not context_use_scrna:
                filename = protein_workdir / f"model_scores_{context}.csv"
                comb_batches_z_prot.to_csv(filename, index=False)

        else:
            comb_batches_z_prot = None

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
