#!/usr/bin/python3
import argparse
import cobra
import cobra.flux_analysis
import copy
import numpy as np
import os
import pandas as pd
import re
import sys
import cobra
import copy
from cobra.flux_analysis import moma, pfba
from project import configs
# from instruments import fetch_entrez_gene_id

import async_bioservices
from async_bioservices import database_convert
from async_bioservices.input_database import InputDatabase
from async_bioservices.output_database import OutputDatabase


def knock_out_simulation(model, inhibitors_filepath, drug_db, ref_flux_file, test_all, pars_flag):
    if ref_flux_file is not None:
        try:
            ref_flux_df = pd.read_csv(ref_flux_file)
            ref_flux_df.set_index("rxn", inplace=True)
            ref_flux = ref_flux_df["flux"].squeeze()
        except FileNotFoundError:
            print(f"Reference flux file not found at {ref_flux_file}")
            sys.exit()
        except KeyError:
            print("Reference flux file must be a CSV file with the columns 'rxn' and 'flux' and row number equal to "
                  "the given context-specific model!")
            sys.exit()

        ref_sol = cobra.core.solution.Solution(model.objective, "OPTIMAL", ref_flux)
    else:
        if pars_flag:
            ref_sol = cobra.flux_analysis.pfba(model)
        else:
            ref_sol = model.optimize()

    if os.path.isfile(inhibitors_filepath):
        print(f"Inhibitors file found at:\n{inhibitors_filepath}")
        DT_genes = pd.read_csv(os.path.join(configs.datadir, inhibitors_filepath), header=None, sep="\t")
        DT_genes.rename(columns={0: "Gene ID"}, inplace=True)
        DT_genes["Gene ID"] = DT_genes["Gene ID"].astype(str)
    else:
        # if inhibitors file does not exist, create a new one
        # only keep inhibitor
        drug_db = drug_db[drug_db["MOA"].str.lower().str.contains("inhibitor") == True]
        DT_genes = pd.DataFrame(columns=["Gene ID"])
        DT_genes["Gene ID"] = drug_db["ENTREZ_GENE_ID"].astype(str)
        DT_genes.replace("-", np.nan, inplace=True)
        DT_genes.dropna(axis=0, inplace=True)
        DT_genes.to_csv(inhibitors_filepath, header=False, sep="\t")
        print(f"Inhibitors file written to:\n{inhibitors_filepath}")

    gene_ind2genes = [x.id for x in model.genes]
    gene_ind2genes = set(gene_ind2genes)
    print(f"{len(gene_ind2genes)} total genes in model.")
    DT_model = list(set(DT_genes["Gene ID"].tolist()).intersection(gene_ind2genes))
    print(f"{len(DT_model)} drug target genes in model")

    model_opt = cobra.flux_analysis.moma(model, solution=ref_sol, linear=False).to_frame()
    model_opt[abs(model_opt) < 1e-8] = 0.0

    has_effects_gene = []
    for id in DT_model:
        gene = model.genes.get_by_id(id)
        for rxn in gene.reactions:
            gene_reaction_rule = rxn.gene_reaction_rule
            gene_ids = re.findall(r"\d+", gene_reaction_rule)
            for gene_id in gene_ids:
                if gene_id == id:
                    boolval = "False"
                else:
                    # boolval = "{}".format(model.genes.get_by_id(gene_id)._functional)
                    boolval = "{}".format(model.genes.get_by_id(gene_id).functional)
                gene_reaction_rule = gene_reaction_rule.replace(
                    "{}".format(gene_id), boolval, 1
                )
            if not eval(gene_reaction_rule) or test_all:
                has_effects_gene.append(id)
                break
    print(f"{len(has_effects_gene)} drug target genes with metabolic effects in model")
    flux_solution = pd.DataFrame()
    for id in has_effects_gene:
        print(f"Peforming knock-out simulation for {id}")
        model_cp = copy.deepcopy(model)  # using model_opt instead bc it makes more sense?
        gene = model_cp.genes.get_by_id(id)
        gene.knock_out()
        opt_model = cobra.flux_analysis.moma(model_cp, solution=ref_sol, linear=False).to_frame()
        flux_solution[id] = opt_model["fluxes"]
        del model_cp

    # flux_solution
    flux_solution[abs(flux_solution) < 1e-8] = 0.0

    flux_solution_ratios = flux_solution.div(model_opt["fluxes"], axis=0)  # ko / original : inf means
    # flux_solution_ratios
    flux_solution_diffs = flux_solution.sub(model_opt["fluxes"], axis=0)  # ko - original
    # flux_solution_diffs

    # has_effects_gene
    return (
        model,
        gene_ind2genes,
        has_effects_gene,
        flux_solution,
        flux_solution_ratios,
        flux_solution_diffs,
    )


def create_gene_pairs(
        datadir,
        model,
        gene_ind2genes,
        flux_solution,
        flux_solution_ratios,
        flux_solution_diffs,
        has_effects_gene,
        disease_down,
):
    disease_down = pd.read_csv(os.path.join(datadir, disease_down))
    DAG_dis_genes = pd.DataFrame()  # data analysis genes
    DAG_dis_genes["Gene ID"] = disease_down.iloc[:, 0].astype(str)
    # DAG_dis_genes
    DAG_dis_met_genes = set(DAG_dis_genes["Gene ID"].tolist()).intersection(
        gene_ind2genes
    )
    # DAG_dis_met_genes

    DAG_dis_met_rxn_ind = []
    gene_i = []
    for id in DAG_dis_met_genes:
        gene = model.genes.get_by_id(id)
        for rxn in gene.reactions:
            DAG_dis_met_rxn_ind.append(rxn.id)
            gene_i.append(id)

    # DAG_dis_met_rxn_ind
    gene_df = pd.DataFrame(gene_i, columns=["Gene IDs"], index=DAG_dis_met_rxn_ind)
    # gene_df

    dag_rxn_flux_ratio = flux_solution_ratios.loc[DAG_dis_met_rxn_ind]
    dag_rxn_flux_diffs = flux_solution_diffs.loc[DAG_dis_met_rxn_ind]
    dag_rxn_flux_value = flux_solution.loc[DAG_dis_met_rxn_ind]
    # dag_rxn_flux_ratio

    gene_mat_out = []
    # gene_i = DAG_dis_met_genes
    # Rind_i = DAG_dis_met_rxn_ind
    for id in has_effects_gene:
        pegene = pd.DataFrame()
        pegene["Gene IDs"] = gene_df["Gene IDs"].copy()
        pegene["rxn_fluxRatio"] = dag_rxn_flux_ratio[id].copy()
        rxn_fluxDiffs = dag_rxn_flux_diffs[id].copy()
        rxn_fluxValue = dag_rxn_flux_value[id].copy()
        pegene["Gene"] = id
        pegene = pegene.loc[
            (~pegene["rxn_fluxRatio"].isna())
            & (abs(rxn_fluxDiffs) + abs(rxn_fluxValue) > 1e-8)
            ]
        # pegene.dropna(axis=0,subset=['rxn_fluxRatio'],inplace=True)
        pegene.index.name = "reaction"
        pegene.reset_index(drop=False, inplace=True)
        gene_mat_out.append(pegene)

    gene_pairs = pd.concat(gene_mat_out, ignore_index=True)
    return gene_pairs


def score_gene_pairs(gene_pairs, filename, input_reg):
    p_model_genes = gene_pairs.Gene.unique()
    d_score = pd.DataFrame([], columns=["score"])
    for p_gene in p_model_genes:
        data_p = gene_pairs.loc[gene_pairs["Gene"] == p_gene].copy()
        total_aff = data_p["Gene IDs"].unique().size
        n_aff_down = (data_p.loc[abs(data_p["rxn_fluxRatio"]) < 0.9, "Gene IDs"].unique().size)
        n_aff_up = (data_p.loc[abs(data_p["rxn_fluxRatio"]) > 1.1, "Gene IDs"].unique().size)
        if input_reg == "up":
            d_s = (n_aff_down - n_aff_up) / total_aff
        else:
            d_s = (n_aff_up - n_aff_down) / total_aff

        d_score.at[p_gene, "score"] = d_s

    d_score.index.name = "Gene"
    d_score.to_csv(os.path.join(configs.datadir, filename))
    return d_score


def score_gene_pairs_diff(gene_pairs, file_full_path):
    p_model_genes = gene_pairs.Gene.unique()
    d_score = pd.DataFrame([], columns=["score"])
    for p_gene in p_model_genes:
        data_p = gene_pairs.loc[gene_pairs["Gene"] == p_gene].copy()
        total_aff = data_p["Gene IDs"].unique().size
        n_aff_down = (
            data_p.loc[data_p["rxn_fluxRatio"] < -1e-8, "Gene IDs"].unique().size
        )
        n_aff_up = data_p.loc[data_p["rxn_fluxRatio"] > 1e-8, "Gene IDs"].unique().size
        d_s = (n_aff_down - n_aff_up) / total_aff
        d_score.at[p_gene, "score"] = d_s

    d_score.index.name = "Gene"
    d_score.to_csv(file_full_path)
    return d_score


def load_Inhi_Fratio(filepath):
    temp2 = pd.read_csv(filepath)
    temp2.rename(
        columns={
            "gene_mat_out1": "Gene",
            "gene_mat_out2": "Gene IDs",
            "gene_mat_out3": "rxn_fluxRatio",
        },
        inplace=True,
    )
    temp2.Gene = temp2.Gene.astype(str)
    temp2["Gene IDs"] = temp2["Gene IDs"].astype(str)
    return temp2


def repurposing_hub_preproc(drug_file):
    drug_db = pd.read_csv(drug_file, sep="\t")
    drug_db_new = pd.DataFrame()
    for index, row in drug_db.iterrows():
        if pd.isnull(row["target"]):
            continue
        for target in row["target"].split("|"):
            drug_db_new = pd.concat(
                [drug_db_new, pd.DataFrame([{"Name": row["pert_iname"],
                                             "MOA": row["moa"],
                                             "Target": target.strip(),
                                             "Phase": row["clinical_phase"]
                                             }])],
                ignore_index=True
            )

    entrez_ids = database_convert.fetch_gene_info(
        input_values=drug_db_new["Target"].tolist(),
        input_db=InputDatabase.GENE_SYMBOL,
        output_db=[
            OutputDatabase.GENE_ID,
            OutputDatabase.ENSEMBL_GENE_ID
        ]
    )

    # entrez_ids = fetch_entrez_gene_id(drug_db_new["Target"].tolist(), input_db="Gene Symbol")
    entrez_ids.reset_index(drop=False, inplace=True)
    drug_db_new["ENTREZ_GENE_ID"] = entrez_ids["Gene ID"]
    drug_db_new = drug_db_new[["Name", "MOA", "Target", "ENTREZ_GENE_ID", "Phase"]]
    return drug_db_new


def drug_repurposing(drug_db, d_score):
    d_score["Gene"] = d_score["Gene"].astype(str)

    d_score_gene_sym = database_convert.fetch_gene_info(
        input_values=d_score["Gene"].tolist(),
        input_db=InputDatabase.GENE_ID,
        output_db=[OutputDatabase.GENE_SYMBOL]
    )

    # d_score_gene_sym = fetch_entrez_gene_id(d_score["Gene"].tolist(), input_db="Gene ID", output_db=["Gene Symbol"])
    d_score.set_index("Gene", inplace=True)
    d_score["Gene Symbol"] = d_score_gene_sym["Gene Symbol"]
    d_score.reset_index(drop=False, inplace=True)
    d_score_new = pd.DataFrame()
    for index, row in d_score.iterrows():
        target = row["Gene Symbol"]
        drugs = drug_db.loc[drug_db["Target"] == target, :]
        drugs.is_copy = False
        drugs["d_score"] = row["score"]
        d_score_new = pd.concat([d_score_new, drugs], ignore_index=True)

    d_score_new.drop_duplicates(inplace=True)
    d_score_trim = d_score_new[
        d_score_new["MOA"].str.lower().str.contains("inhibitor") == True
        ]

    return d_score_trim


def main(argv):
    drug_raw_file = "Repurposing_Hub_export.txt"

    parser = argparse.ArgumentParser(
        prog="knock_out_simulation.py",
        description="This script is responsible for mapping drug targets in metabolic models, performing knock out simulations, and comparing simulation results with disease genes. It also identified drug targets and repurposable drugs.",
        epilog="For additional help, please post questions/issues in the MADRID GitHub repo at "
               "https://github.com/HelikarLab/MADRID or email babessell@gmail.com",
    )
    parser.add_argument(
        "-m",
        "--context-model",
        type=str,
        required=True,
        dest="model",
        help="The context-specific model file, (must be .mat, .xml, or .json"
    )
    parser.add_argument(
        "-c",
        "--context-name",
        type=str,
        required=True,
        dest="context",
        help="Name of context, tissue, cell-type, etc"
    )
    parser.add_argument(
        "-d",
        "--disease-name",
        type=str,
        required=True,
        dest="disease",
        help="Name of disease"
    )
    parser.add_argument(
        "-up",
        "--disease-up",
        type=str,
        required=True,
        dest="disease_up",
        help="The name of the disease up-regulated file"
    )
    parser.add_argument(
        "-dn",
        "--disease-down",
        type=str,
        required=True,
        dest="disease_down",
        help="The name of the disease down-regulated file"
    )
    parser.add_argument(
        "-r",
        "--raw-drug-file",
        type=str,
        required=True,
        dest="raw_drug_file",
        help="The name of the raw drug file"
    )
    parser.add_argument(
        "-f",
        "--reference-flux-file",
        type=str if ("--reference-flux-file" in argv or "-f" in argv) else type(None),
        required=False,
        default=None,
        dest="ref_flux_file",
        help="The name of the reference flux file"
    )
    parser.add_argument(
        "-a",
        "--test-all",
        action="store_true",
        required=False,
        default=False,
        dest="test_all",
        help="Test all genes, even ones predicted to have little no effect."
    )
    parser.add_argument(
        "-p",
        "--parsimonious",
        action="store_true",
        required=False,
        default=False,
        dest="pars_flag",
        help="Use parsimonious FBA for optimal reference solution (only if not providing flux file)"
    )

    args = parser.parse_args()
    tissue_spec_model_file = args.model
    context = args.context
    disease = args.disease
    disease_up_file = args.disease_up
    disease_down_file = args.disease_down
    drug_raw_file = args.raw_drug_file
    ref_flux_file = args.ref_flux_file
    test_all = args.test_all
    pars_flag = args.pars_flag

    output_dir = os.path.join(configs.datadir, "results", context, disease)
    inhibitors_file = os.path.join(output_dir, f"{context}_{disease}_inhibitors.tsv")

    print(f"Output directory: '{output_dir}'")
    print(f"Tissue Specific Model file is at: {tissue_spec_model_file}")
    print(f"Tissue specific inhibitors is at: {inhibitors_file}")

    if tissue_spec_model_file[-4:] == ".mat":
        cobra_model = cobra.io.load_matlab_model(tissue_spec_model_file)
    elif tissue_spec_model_file[-4:] == ".xml":
        cobra_model = cobra.io.load_sbml_model(tissue_spec_model_file)
    elif tissue_spec_model_file[-5:] == ".json":
        cobra_model = cobra.io.load_json_model(tissue_spec_model_file)
    else:
        raise NameError("reference model format must be .xml, .mat, or .json")

    cobra_model.solver = "gurobi"

    # preprocess repurposing hub data
    drug_tsv_file = "Repurposing_Hub_Preproc.tsv"
    drug_raw_file = os.path.join(configs.datadir, drug_raw_file)
    drug_file = os.path.join(configs.datadir, drug_tsv_file)
    if not os.path.isfile(drug_file):
        print("Preprocessing raw Repurposing Hub DB file...")
        drug_db = repurposing_hub_preproc(drug_raw_file)
        drug_db.to_csv(drug_file, index=False, sep="\t")
        print(f"Preprocessed Repurposing Hub tsv file written to:\n{drug_file}")
    else:
        print(f"Found preprocessed Repurposing Hub tsv file at:\n{drug_file}")
        drug_db = pd.read_csv(drug_file, sep="\t")

    # Knock Out Simulation
    model, gene_ind2genes, has_effects_gene, fluxsolution, flux_solution_ratios, flux_solution_diffs = \
        knock_out_simulation(
            model=cobra_model,
            inhibitors_filepath=inhibitors_file,
            drug_db=drug_db,
            ref_flux_file=ref_flux_file,
            test_all=test_all,
            pars_flag=pars_flag
        )

    flux_solution_diffs.to_csv(os.path.join(output_dir, "flux_diffs_KO.csv"))
    flux_solution_ratios.to_csv(os.path.join(output_dir, "flux_ratios_KO.csv"))

    gene_pairs_down = create_gene_pairs(
        configs.datadir,
        model,
        gene_ind2genes,
        fluxsolution,
        flux_solution_ratios,
        flux_solution_diffs,
        has_effects_gene,
        disease_down=disease_down_file,
    )

    gene_pairs_down.to_csv(os.path.join(output_dir, "Gene_Pairs_Inhi_Fratio_DOWN.txt"), index=False)

    gene_pairs_up = create_gene_pairs(
        configs.datadir,
        model,
        gene_ind2genes,
        fluxsolution,
        flux_solution_ratios,
        flux_solution_diffs,
        has_effects_gene,
        disease_down=disease_up_file,
    )
    gene_pairs_up.to_csv(os.path.join(output_dir, "Gene_Pairs_Inhi_Fratio_UP.txt"), index=False)
    d_score_down = score_gene_pairs(gene_pairs_down, os.path.join(output_dir, "d_score_DOWN.csv"), input_reg="down")
    d_score_up = score_gene_pairs(gene_pairs_up, os.path.join(output_dir, "d_score_UP.csv"), input_reg="up")
    pertubation_effect_score = (d_score_up + d_score_down).sort_values(by="score", ascending=False)
    pertubation_effect_score.to_csv(os.path.join(output_dir, "d_score.csv"))
    pertubation_effect_score.reset_index(drop=False, inplace=True)

    # last step: output drugs based on d score
    drug_score = drug_repurposing(drug_db, pertubation_effect_score)
    drug_score_file = os.path.join(output_dir, "drug_score.csv")
    drug_score.to_csv(drug_score_file, index=False)
    print("Gene D score mapped to repurposing drugs saved to\n{}".format(drug_score_file))


if __name__ == "__main__":
    main(sys.argv[1:])
