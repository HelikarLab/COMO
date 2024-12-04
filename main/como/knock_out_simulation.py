# ruff: noqa

import argparse
import os
import re
import sys
from concurrent.futures import Future, ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Union

import cobra
import numpy as np
import pandas as pd
from fast_bioservices import BioDBNet, Input, Output
from project import Config

configs = Config()


@dataclass
class KnockoutResults:
    model: cobra.Model
    gene_ind2genes: set[str]
    genes_with_metabolic_effects: list[str]
    flux_solution: pd.DataFrame
    flux_solution_ratios: pd.DataFrame
    flux_solution_diffs: pd.DataFrame


def _perform_knockout(
    model: cobra.Model,
    gene_id: str,
    reference_solution,
) -> tuple[str, pd.Series]:
    """This function will perform a single gene knockout. It will be used in multiprocessing"""
    with model:
        gene: cobra.Gene = model.genes.get_by_id(gene_id)
        gene.knock_out()
        optimized_model: cobra.Solution = cobra.flux_analysis.moma(model, solution=reference_solution, linear=False)
    return gene_id, optimized_model.fluxes


def knock_out_simulation(
    model: cobra.Model,
    inhibitors_filepath: Path,
    drug_db: pd.DataFrame,
    reference_flux_filepath: Union[str, Path, None],
    test_all: bool,
    pars_flag: bool,
) -> KnockoutResults:
    reference_solution: cobra.Solution
    if reference_flux_filepath is not None:
        reference_flux_filepath: Path = Path(reference_flux_filepath)
        if not reference_flux_filepath.exists():
            raise FileNotFoundError(f"Reference flux file not found at {reference_flux_filepath.as_posix()}")
        reference_flux_df: pd.DataFrame = pd.read_csv(reference_flux_filepath)
        if "rxn" not in reference_flux_df.columns or "flux" not in reference_flux_df.columns:
            raise KeyError("Reference flux file must be a CSV file with the columns 'rxn' and 'flux' with the same number of rows as the number of reactions in the given context-specific model!")  # fmt: skip
        reference_flux_df.set_index("rxn", inplace=True)
        reference_flux = reference_flux_df["flux"].squeeze()
        reference_solution = cobra.core.solution.Solution(model.objective, "OPTIMAL", reference_flux)  # fmt: skip
    else:
        reference_solution = cobra.flux_analysis.pfba(model) if pars_flag else model.optimize()

    drug_target_genes: pd.DataFrame
    if inhibitors_filepath.exists():
        print(f"Inhibitors file found at: {inhibitors_filepath}")
        drug_target_genes = pd.read_csv(inhibitors_filepath, sep="\t")
        # dt_genes.rename(columns={0: "Gene ID"}, inplace=True)
        drug_target_genes["Gene ID"] = drug_target_genes["Gene ID"].astype(str)
    else:
        # only keep inhibitors
        drug_db = drug_db[drug_db["moa"].str.lower().str.contains("inhibitor")]
        drug_target_genes = pd.DataFrame(columns=["Gene ID"])
        drug_target_genes["Gene ID"] = drug_db["Gene ID"].astype(str)
        drug_target_genes.replace("-", np.nan, inplace=True)
        drug_target_genes.dropna(axis=0, inplace=True)
        drug_target_genes.to_csv(inhibitors_filepath, header=True, sep="\t", index=False)
        print(f"Inhibitors file written to: {inhibitors_filepath}")

    gene_ind2genes = set(x.id for x in model.genes)
    dt_model = list(set(drug_target_genes["Gene ID"].tolist()).intersection(gene_ind2genes))
    print(f"{len(gene_ind2genes)} genes in model, {len(dt_model)} can be targeted by inhibitors")

    wild_type_model = cobra.flux_analysis.moma(model, solution=reference_solution).to_frame()
    wild_type_model[abs(wild_type_model) < 1e-6] = 0.0

    genes_with_metabolic_effects = []
    for id_ in dt_model:
        gene: cobra.Gene = model.genes.get_by_id(id_)
        for rxn in gene.reactions:
            gene_reaction_rule = rxn.gene_reaction_rule
            gene_ids = re.findall(r"\d+", gene_reaction_rule)
            for gene_id in gene_ids:
                boolval = "False" if gene_id == id_ else str(model.genes.get_by_id(gene_id).functional)
                gene_reaction_rule = gene_reaction_rule.replace(gene_id, boolval, 1)
            if not eval(gene_reaction_rule) or test_all:
                genes_with_metabolic_effects.append(id_)
                break
    print(f"Found {len(genes_with_metabolic_effects)} genes with potentially-significant metabolic impacts")  # fmt: skip

    futures: list[Future[tuple[str, pd.Series]]] = []
    with ProcessPoolExecutor(max_workers=os.cpu_count()) as executor:
        for i, id_ in enumerate(genes_with_metabolic_effects, start=1):
            future: Future = executor.submit(_perform_knockout, model, id_, reference_solution)
            futures.append(future)

        gene_id: str
        knock_out_flux: pd.Series
        flux_solution: pd.DataFrame = pd.DataFrame()
        for result in as_completed(futures):
            gene_id, knock_out_flux = result.result()
            flux_solution[gene_id] = knock_out_flux

    flux_solution[abs(flux_solution) < 1e-6] = 0.0
    flux_solution_ratios = flux_solution.div(wild_type_model["fluxes"], axis=0)
    flux_solution_diffs = flux_solution.sub(wild_type_model["fluxes"], axis=0)

    return KnockoutResults(
        model=model,
        gene_ind2genes=gene_ind2genes,
        genes_with_metabolic_effects=genes_with_metabolic_effects,
        flux_solution=flux_solution,
        flux_solution_ratios=flux_solution_ratios,
        flux_solution_diffs=flux_solution_diffs,
    )


def create_gene_pairs(
    datadir: Union[str, Path],
    model: cobra.Model,
    gene_ind2genes: set[str],
    flux_solution: pd.DataFrame,
    flux_solution_ratios: pd.DataFrame,
    flux_solution_diffs: pd.DataFrame,
    has_effects_gene: list[str],
    disease_genes_filename: Path,
):
    disease_genes_df: pd.DataFrame = pd.read_csv(str(os.path.join(datadir, disease_genes_filename)))
    if len(disease_genes_df.columns) != 1:
        raise ValueError(f"Expected 1 column in {disease_genes_filename}, got {len(disease_genes_df.columns)}")

    disease_genes_df.columns = ["Gene ID"]
    disease_genes_df["Gene ID"] = disease_genes_df["Gene ID"].astype(str)
    metabolic_disease_genes = set(disease_genes_df["Gene ID"]).intersection(gene_ind2genes)

    gene_df = pd.DataFrame(columns=["Gene ID", "Reaction ID"])
    for id_ in metabolic_disease_genes:
        model_gene: cobra.Gene = model.genes.get_by_id(id_)
        gene_reactions: list[str] = [rxn.id for rxn in model_gene.reactions]
        gene_df = pd.concat([gene_df, pd.DataFrame({"Gene ID": id_, "Reaction ID": gene_reactions})], ignore_index=True)
    gene_df.set_index("Reaction ID", drop=True, inplace=True)

    dag_rxn_flux_ratio: pd.DataFrame = flux_solution_ratios.loc[gene_df.index.tolist()]
    dag_rxn_flux_diffs: pd.DataFrame = flux_solution_diffs.loc[gene_df.index.tolist()]
    dag_rxn_flux_value: pd.DataFrame = flux_solution.loc[gene_df.index.tolist()]
    gene_mat_out: list[pd.DataFrame] = []

    for id_ in has_effects_gene:
        pegene = pd.DataFrame()
        pegene["Gene ID"] = gene_df["Gene ID"]
        pegene["rxn_flux_ratio"] = dag_rxn_flux_ratio[id_]
        pegene["Gene"] = id_

        rxn_flux_diffs = dag_rxn_flux_diffs[id_]
        rxn_flux_value = dag_rxn_flux_value[id_]
        pegene = pegene.loc[(~pegene["rxn_flux_ratio"].isna()) & (abs(rxn_flux_diffs) + abs(rxn_flux_value) > 1e-8)]
        pegene.index.name = "reaction"
        gene_mat_out.append(pegene)

    gene_pairs = pd.concat(gene_mat_out, ignore_index=True)
    return gene_pairs


def score_gene_pairs(gene_pairs, filename, input_reg):
    p_model_genes = gene_pairs.Gene.unique()
    d_score = pd.DataFrame([], columns=["score"])
    for p_gene in p_model_genes:
        data_p = gene_pairs.loc[gene_pairs["Gene"] == p_gene].copy()
        total_aff = data_p["Gene ID"].unique().size
        n_aff_down = data_p.loc[abs(data_p["rxn_flux_ratio"]) < 0.9, "Gene ID"].unique().size
        n_aff_up = data_p.loc[abs(data_p["rxn_flux_ratio"]) > 1.1, "Gene ID"].unique().size
        if input_reg == "up":
            d_s = (n_aff_down - n_aff_up) / total_aff
        else:
            d_s = (n_aff_up - n_aff_down) / total_aff

        d_score.at[p_gene, "score"] = d_s

    d_score.index.name = "Gene ID"
    d_score.to_csv(configs.data_dir / filename)
    return d_score


def score_gene_pairs_diff(gene_pairs, file_full_path):
    p_model_genes = gene_pairs.Gene.unique()
    d_score = pd.DataFrame([], columns=["score"])
    for p_gene in p_model_genes:
        data_p = gene_pairs.loc[gene_pairs["Gene"] == p_gene].copy()
        total_aff = data_p["Gene ID"].unique().size
        n_aff_down = data_p.loc[data_p["rxn_flux_ratio"] < -1e-8, "Gene ID"].unique().size
        n_aff_up = data_p.loc[data_p["rxn_flux_ratio"] > 1e-8, "Gene ID"].unique().size
        d_s = (n_aff_down - n_aff_up) / total_aff
        d_score.at[p_gene, "score"] = d_s

    d_score.index.name = "Gene"
    d_score.to_csv(file_full_path)
    return d_score


def repurposing_hub_preproc(drug_info_filepath: Path, biodbnet: BioDBNet):
    drug_info_df: pd.DataFrame = pd.read_csv(drug_info_filepath, sep="\t")
    drug_info_df["target"] = drug_info_df["target"].str.split("|").explode().reset_index(drop=True)
    drug_info_df = (
        drug_info_df.drop(columns=["disease_area", "indication"])
        .rename(columns={"pert_iname": "name", "clinical_phase": "phase"})
        .dropna(subset=["target", "moa"])
    )
    # for index, row in drug_db.iterrows():
    #     if pd.isnull(row["target"]):
    #         continue
    #     for target in row["target"].split("|"):
    #         drug_db_new = pd.concat(
    #             [
    #                 drug_db_new,
    #                 pd.DataFrame(
    #                     [
    #                         {
    #                             "Name": row["pert_iname"],
    #                             "MOA": row["moa"],
    #                             "Target": target.strip(),
    #                             "Phase": row["clinical_phase"],
    #                         }
    #                     ]
    #                 ),
    #             ],
    #             ignore_index=True,
    #         )
    # drug_db_new.reset_index(inplace=True)
    entrez_ids = biodbnet.db2db(
        input_values=drug_info_df["target"].tolist(),
        input_db=Input.GENE_SYMBOL,
        output_db=Output.GENE_ID,
    )
    entrez_ids.rename(columns={"Gene Symbol": "target"}, inplace=True)
    drug_info_df = pd.merge(drug_info_df, entrez_ids, on="target")
    # entrez_ids.reset_index(drop=False, inplace=True)
    # drug_db_new["ENTREZ_GENE_ID"] = entrez_ids["Gene ID"]
    # drug_db_new = drug_db_new[["Name", "MOA", "Target", "ENTREZ_GENE_ID", "Phase"]]
    return drug_info_df


def drug_repurposing(drug_db: pd.DataFrame, perturbation_score: pd.DataFrame, biodbnet: BioDBNet):
    perturbation_score["Gene ID"] = perturbation_score["Gene ID"].astype(str)

    conversion = biodbnet.db2db(
        input_values=perturbation_score["Gene ID"].tolist(),
        input_db=Input.GENE_ID,
        output_db=[Output.GENE_SYMBOL],
    )

    perturbation_score = pd.merge(perturbation_score, conversion, on="Gene ID", how="left")
    # d_score.set_index("Gene", inplace=True)
    # d_score["Gene Symbol"] = d_score_gene_sym["Gene Symbol"]
    # d_score.reset_index(drop=False, inplace=True)
    drug_scores = pd.DataFrame()
    for index, row in perturbation_score.iterrows():
        target = row["Gene Symbol"]
        drugs = drug_db.loc[drug_db["target"] == target, :].copy()  # Use `.copy()` to prevent `SettingWithCopyWarning`
        drugs["score"] = row["score"]
        drug_scores = pd.concat([drug_scores, drugs], ignore_index=True)

    drug_scores.drop_duplicates(inplace=True)
    drug_scores = drug_scores[drug_scores["moa"].str.lower().str.contains("inhibitor")]
    return drug_scores


def main(argv):
    parser = argparse.ArgumentParser(
        prog="knock_out_simulation.py",
        description="This script is responsible for mapping drug targets in metabolic models, performing knock out simulations, and comparing simulation results with disease genes. It also identified drug targets and repurposable drugs.",
        epilog="For additional help, please post questions/issues in the MADRID GitHub repo at "
        "https://github.com/HelikarLab/COMO",
    )
    parser.add_argument(
        "-m",
        "--context-model",
        type=str,
        required=True,
        dest="model",
        help="The context-specific model file, (must be .mat, .xml, or .json",
    )
    parser.add_argument(
        "-c",
        "--context-name",
        type=str,
        required=True,
        dest="context",
        help="Name of context, tissue, cell-type, etc",
    )
    parser.add_argument(
        "-d",
        "--disease-name",
        type=str,
        required=True,
        dest="disease",
        help="Name of disease",
    )
    parser.add_argument(
        "-up",
        "--disease-up",
        type=str,
        required=True,
        dest="disease_up",
        help="The name of the disease up-regulated file",
    )
    parser.add_argument(
        "-dn",
        "--disease-down",
        type=str,
        required=True,
        dest="disease_down",
        help="The name of the disease down-regulated file",
    )
    parser.add_argument(
        "-r",
        "--raw-drug-file",
        type=str,
        required=True,
        dest="raw_drug_file",
        help="The name of the raw drug file",
    )
    parser.add_argument(
        "-f",
        "--reference-flux-file",
        type=str if ("--reference-flux-file" in argv or "-f" in argv) else type(None),
        required=False,
        default=None,
        dest="ref_flux_file",
        help="The name of the reference flux file",
    )
    parser.add_argument(
        "-a",
        "--test-all",
        action="store_true",
        required=False,
        default=False,
        dest="test_all",
        help="Test all genes, even ones predicted to have little no effect.",
    )
    parser.add_argument(
        "-p",
        "--parsimonious",
        action="store_true",
        required=False,
        default=False,
        dest="pars_flag",
        help="Use parsimonious FBA for optimal reference solution (only if not providing flux file)",
    )
    parser.add_argument(
        "-s",
        "--solver",
        type=str,
        required=False,
        default="gurobi",
        dest="solver",
        help="The solver to use for FBA. Options are: gurobi or glpk",
    )

    args = parser.parse_args()
    tissue_spec_model_file = Path(args.model)
    context = args.context
    disease = args.disease
    disease_up_file = Path(args.disease_up)
    disease_down_file = Path(args.disease_down)
    raw_drug_filename = args.raw_drug_file
    ref_flux_file = args.ref_flux_file
    test_all = args.test_all
    pars_flag = args.pars_flag
    solver = args.solver

    output_dir = Path(configs.data_dir, "results", context, disease)
    inhibitors_filepath = Path(output_dir, f"{context}_{disease}_inhibitors.tsv")
    biodbnet = BioDBNet(cache=False)
    thread_pool = ThreadPoolExecutor(max_workers=1)

    print(f"Output directory: '{output_dir.as_posix()}'")
    print(f"Tissue Specific Model file is at: {tissue_spec_model_file.as_posix()}")
    print(f"Tissue specific inhibitors is at: {inhibitors_filepath.as_posix()}")

    if not tissue_spec_model_file.exists():
        raise FileNotFoundError(f"Model file not found at {tissue_spec_model_file.as_posix()}")
    elif tissue_spec_model_file.suffix == ".mat":
        future = thread_pool.submit(cobra.io.load_matlab_model, infile_path=tissue_spec_model_file.as_posix())  # type: ignore
    elif tissue_spec_model_file.suffix in (".xml", ".sbml"):
        future = thread_pool.submit(cobra.io.read_sbml_model, filename=tissue_spec_model_file.as_posix())  # type: ignore
    elif tissue_spec_model_file.suffix == ".json":
        future = thread_pool.submit(cobra.io.load_json_model, filename=tissue_spec_model_file.as_posix())  # type: ignore
    else:
        raise NameError("Reference model  must be in 'mat', 'xml', 'sbml', or 'json' format.")

    raw_drug_filepath = Path(configs.data_dir, raw_drug_filename)
    reformatted_drug_filepath = raw_drug_filepath.with_stem(f"{raw_drug_filepath.stem}_processed")
    drug_info_df: pd.DataFrame
    if reformatted_drug_filepath.exists():
        print(f"Found preprocessed Repurposing Hub tsv file at: {reformatted_drug_filepath}")
        drug_info_df = pd.read_csv(reformatted_drug_filepath, sep="\t")
    else:
        print("Preprocessing raw Repurposing Hub DB file...")
        drug_info_df = repurposing_hub_preproc(drug_info_filepath=raw_drug_filepath, biodbnet=biodbnet)
        drug_info_df.to_csv(reformatted_drug_filepath, index=False, sep="\t")
        print(f"Preprocessed Repurposing Hub tsv file written to: {reformatted_drug_filepath.as_posix()}")

    cobra_model: cobra.Model = future.result()
    cobra_model.solver = solver
    thread_pool.shutdown()

    knockout_results = knock_out_simulation(
        model=cobra_model,
        inhibitors_filepath=inhibitors_filepath,
        drug_db=drug_info_df,
        reference_flux_filepath=ref_flux_file,
        test_all=test_all,
        pars_flag=pars_flag,
    )

    knockout_results.flux_solution_diffs.to_csv(output_dir / "flux_diffs_KO.csv")
    knockout_results.flux_solution_ratios.to_csv(output_dir / "flux_ratios_KO.csv")

    gene_pairs_down = create_gene_pairs(
        configs.data_dir,
        knockout_results.model,
        knockout_results.gene_ind2genes,
        knockout_results.flux_solution,
        knockout_results.flux_solution_ratios,
        knockout_results.flux_solution_diffs,
        knockout_results.genes_with_metabolic_effects,
        disease_genes_filename=disease_down_file,
    )
    gene_pairs_down.to_csv(os.path.join(output_dir, f"{context}_Gene_Pairs_Inhi_Fratio_DOWN.txt"), index=False)

    gene_pairs_up = create_gene_pairs(
        configs.data_dir,
        knockout_results.model,
        knockout_results.gene_ind2genes,
        knockout_results.flux_solution,
        knockout_results.flux_solution_ratios,
        knockout_results.flux_solution_diffs,
        knockout_results.genes_with_metabolic_effects,
        disease_genes_filename=disease_up_file,
    )
    gene_pairs_up.to_csv(os.path.join(output_dir, f"{context}_Gene_Pairs_Inhi_Fratio_UP.txt"), index=False)

    d_score_down = score_gene_pairs(
        gene_pairs_down,
        os.path.join(output_dir, f"{context}_d_score_DOWN.csv"),
        input_reg="down",
    )
    d_score_up = score_gene_pairs(
        gene_pairs_up,
        os.path.join(output_dir, f"{context}_d_score_UP.csv"),
        input_reg="up",
    )
    perturbation_score: pd.DataFrame = (d_score_up + d_score_down).sort_values(by="score", ascending=False)
    perturbation_score.to_csv(os.path.join(output_dir, f"{context}_d_score.csv"))
    perturbation_score.reset_index(drop=False, inplace=True)

    drug_score = drug_repurposing(drug_db=drug_info_df, perturbation_score=perturbation_score, biodbnet=biodbnet)
    drug_score_file = os.path.join(output_dir, f"{context}_drug_score.csv")
    drug_score.to_csv(drug_score_file, index=False)
    print(f"Gene D score mapped to repurposing drugs saved to {drug_score_file}")

    print(f"\nFinished {disease}!")


if __name__ == "__main__":
    main(sys.argv[1:])
