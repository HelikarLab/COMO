#!/usr/bin/python3
import argparse
import os
import re
import sys
import getopt
import time
import pandas as pd
import numpy as np
import cobra
import copy
from cobra.flux_analysis import (
    single_gene_deletion,
    single_reaction_deletion,
    double_gene_deletion,
    double_reaction_deletion,
    moma,
)
from project import configs

# from transcriptomic_gen import *
# from proteomics_gen import *
from instruments import fetch_entrez_gene_id


def knock_out_simulation(datadir, model_file, inhibitors, drugDB):

    if model_file[-4:] == ".xml":
        model = cobra.io.read_sbml_model(model_file)
    elif model_file[-4:] == ".mat":
        model = cobra.io.load_matlab_model(model_file)
    elif model_file[-5:] == ".json":
        model = cobra.io.load_json_model(model_file)
    else:
        print("Unsupported File Format of Model: {}".format(model_file))
        return None

    inhibitorsFullpath = os.path.join(datadir, inhibitors)
    if os.path.isfile(inhibitorsFullpath):
        DT_genes = pd.read_csv(os.path.join(datadir, inhibitors), header=None)
        DT_genes.rename(columns={0: "Gene ID"}, inplace=True)
        DT_genes["Gene ID"] = DT_genes["Gene ID"].astype(str)
    else:
        # if inhibitors file does not exist, create a new one
        # only keep inhibitor
        drugDB = drugDB[drugDB["MOA"].str.lower().str.contains("inhibitor") == True]
        DT_genes = pd.DataFrame(columns=["Gene ID"])
        DT_genes["Gene ID"] = drugDB["ENTREZ_GENE_ID"].astype(str)
        DT_genes.replace("-", np.nan, inplace=True)
        DT_genes.dropna(axis=0, inplace=True)
        DT_genes.to_csv(inhibitorsFullpath, header=False)

    geneInd2genes = [x.id for x in model.genes]
    print(len(geneInd2genes))
    geneInd2genes = set(geneInd2genes)
    print(len(geneInd2genes))
    DT_model = set(DT_genes["Gene ID"].tolist()).intersection(geneInd2genes)
    print(len(DT_model))
    DT_model = list(DT_model)

    model_opt = moma(model).to_frame()
    # model_opt = model.optimize().to_frame()
    model_opt[abs(model_opt) < 1e-8] = 0.0

    HasEffects_Gene = []
    for id in DT_model:
        gene = model.genes.get_by_id(id)
        for rxn in gene.reactions:
            gene_reaction_rule = rxn.gene_reaction_rule
            gene_ids = re.findall(r"\d+", gene_reaction_rule)
            for gene_id in gene_ids:
                if gene_id == id:
                    boolval = "False"
                else:
                    boolval = "{}".format(model.genes.get_by_id(gene_id)._functional)
                gene_reaction_rule = gene_reaction_rule.replace(
                    "{}".format(gene_id), boolval, 1
                )
            if not eval(gene_reaction_rule):
                HasEffects_Gene.append(id)
                break

    fluxsolution = pd.DataFrame()
    for id in HasEffects_Gene:
        model_cp = copy.deepcopy(model)
        gene = model_cp.genes.get_by_id(id)
        gene.knock_out()
        opt_model = moma(model_cp).to_frame()
        # opt_model = model_cp.optimize().to_frame() #FBA
        fluxsolution[id] = opt_model["fluxes"]
        del model_cp

    # fluxsolution
    fluxsolution[abs(fluxsolution) < 1e-8] = 0.0

    fluxSolutionRatios = fluxsolution.div(model_opt["fluxes"], axis=0)
    # fluxSolutionRatios
    fluxSolutionDiffs = fluxsolution.sub(model_opt["fluxes"], axis=0)
    # fluxSolutionDiffs

    # HasEffects_Gene
    return (
        model,
        geneInd2genes,
        HasEffects_Gene,
        fluxsolution,
        fluxSolutionRatios,
        fluxSolutionDiffs,
    )


def create_gene_pairs(
    datadir,
    model,
    geneInd2genes,
    fluxsolution,
    fluxSolutionRatios,
    fluxSolutionDiffs,
    HasEffects_Gene,
    Disease_Down,
):
    Disease_down = pd.read_csv(os.path.join(datadir, Disease_Down))
    DAG_dis_genes = pd.DataFrame()
    DAG_dis_genes["Gene ID"] = Disease_down.iloc[:, 0].astype(str)
    # DAG_dis_genes

    DAG_dis_met_genes = set(DAG_dis_genes["Gene ID"].tolist()).intersection(
        geneInd2genes
    )
    # DAG_dis_met_genes

    DAG_dis_met_rxnInd = []
    Gene_i = []
    for id in DAG_dis_met_genes:
        gene = model.genes.get_by_id(id)
        for rxn in gene.reactions:
            DAG_dis_met_rxnInd.append(rxn.id)
            Gene_i.append(id)

    # DAG_dis_met_rxnInd
    Gene_df = pd.DataFrame(Gene_i, columns=["Gene IDs"], index=DAG_dis_met_rxnInd)
    # Gene_df

    DAG_rxn_fluxRatio = fluxSolutionRatios.loc[DAG_dis_met_rxnInd]
    DAG_rxn_fluxDiffs = fluxSolutionDiffs.loc[DAG_dis_met_rxnInd]
    DAG_rxn_fluxValue = fluxsolution.loc[DAG_dis_met_rxnInd]
    # DAG_rxn_fluxRatio

    gene_mat_out = []
    # Gene_i = DAG_dis_met_genes
    # Rind_i = DAG_dis_met_rxnInd
    for id in HasEffects_Gene:
        pegene = pd.DataFrame()
        pegene["Gene IDs"] = Gene_df["Gene IDs"].copy()
        pegene["rxn_fluxRatio"] = DAG_rxn_fluxRatio[id].copy()
        rxn_fluxDiffs = DAG_rxn_fluxDiffs[id].copy()
        rxn_fluxValue = DAG_rxn_fluxValue[id].copy()
        pegene["Gene"] = id
        pegene = pegene.loc[
            (~pegene["rxn_fluxRatio"].isna())
            & (abs(rxn_fluxDiffs) + abs(rxn_fluxValue) > 1e-6)
        ]
        # pegene.dropna(axis=0,subset=['rxn_fluxRatio'],inplace=True)
        pegene.index.name = "reaction"
        pegene.reset_index(drop=False, inplace=True)
        gene_mat_out.append(pegene)

    Gene_Pairs = pd.concat(gene_mat_out, ignore_index=True)
    return Gene_Pairs


def score_gene_pairs(Gene_Pairs, filename):
    p_model_genes = Gene_Pairs.Gene.unique()
    d_score = pd.DataFrame([], columns=["score"])
    for p_gene in p_model_genes:
        data_p = Gene_Pairs.loc[Gene_Pairs["Gene"] == p_gene].copy()
        # print(data_p)
        total_aff = data_p["Gene IDs"].unique().size
        # print(total_aff)
        n_aff_down = (
            data_p.loc[abs(data_p["rxn_fluxRatio"]) < 0.99, "Gene IDs"].unique().size
        )
        # print(n_aff_down)
        n_aff_up = (
            data_p.loc[abs(data_p["rxn_fluxRatio"]) > 1.0, "Gene IDs"].unique().size
        )
        # print(n_aff_up)
        d_s = (n_aff_down - n_aff_up) / total_aff
        # print(d_s)
        d_score.at[p_gene, "score"] = d_s

    d_score.index.name = "Gene"
    d_score.to_csv(os.path.join(configs.datadir, filename))
    return d_score


def score_gene_pairs_diff(Gene_Pairs, fileFullPath):
    p_model_genes = Gene_Pairs.Gene.unique()
    d_score = pd.DataFrame([], columns=["score"])
    for p_gene in p_model_genes:
        data_p = Gene_Pairs.loc[Gene_Pairs["Gene"] == p_gene].copy()
        # print(data_p)
        total_aff = data_p["Gene IDs"].unique().size
        # print(total_aff)
        n_aff_down = (
            data_p.loc[data_p["rxn_fluxRatio"] < -1e-8, "Gene IDs"].unique().size
        )
        # print(n_aff_down)
        n_aff_up = data_p.loc[data_p["rxn_fluxRatio"] > 1e-8, "Gene IDs"].unique().size
        # print(n_aff_up)
        d_s = (n_aff_down - n_aff_up) / total_aff
        # print(d_s)
        d_score.at[p_gene, "score"] = d_s

    d_score.index.name = "Gene"
    d_score.to_csv(fileFullPath)
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


def repurposing_hub_preproc(drugFile):
    drugDB = pd.read_csv(drugFile, sep="\t")
    drugDB_new = pd.DataFrame()
    for index, row in drugDB.iterrows():
        if pd.isnull(row["Target"]):
            continue
        for target in row["Target"].split(","):
            drugDB_new = drugDB_new.append(
                {
                    "Name": row["Name"],
                    "MOA": row["MOA"],
                    "Target": target.strip(),
                    "Phase": row["Phase"],
                },
                ignore_index=True,
            )
    entrez_IDs = fetch_entrez_gene_id(
        drugDB_new["Target"].tolist(), input_db="Gene Symbol"
    )
    entrez_IDs.reset_index(drop=False, inplace=True)
    drugDB_new["ENTREZ_GENE_ID"] = entrez_IDs["Gene ID"]
    drugDB_new = drugDB_new[["Name", "MOA", "Target", "ENTREZ_GENE_ID", "Phase"]]
    return drugDB_new


def drug_repurposing(drugDB, d_score):
    d_score["Gene"] = d_score["Gene"].astype(str)
    d_score_geneSym = fetch_entrez_gene_id(
        d_score["Gene"].tolist(), input_db="Gene ID", output_db=["Gene Symbol"]
    )
    d_score.set_index("Gene", inplace=True)
    d_score["Gene Symbol"] = d_score_geneSym["Gene Symbol"]
    d_score.reset_index(drop=False, inplace=True)
    d_score_new = pd.DataFrame()
    for index, row in d_score.iterrows():
        target = row["Gene Symbol"]
        # print(target)
        drugs = drugDB.loc[drugDB["Target"] == target, :]
        # print(drugs)
        drugs["d score"] = row["score"]
        d_score_new = d_score_new.append(drugs, ignore_index=True)

    d_score_new.drop_duplicates(inplace=True)
    d_score_trim = d_score_new[
        d_score_new["MOA"].str.lower().str.contains("inhibitor") == True
    ]
    return d_score_trim


def main(argv):
    drug_raw_file = "Repurposing_Hub_export.txt"

    # TODO: Fix this description
    parser = argparse.ArgumentParser(
        prog="knock_out_simulation.py",
        description="Description goes here",
        epilog="For additional help, please post questions/issues in the MADRID GitHub repo at "
        "https://github.com/HelikarLab/MADRID or email babessell@gmail.com",
    )
    parser.add_argument(
        "-t",
        "--tissue-file",
        type=str,
        required=True,
        dest="tissue_file",
        help="The name of the input tissue file",
    )
    parser.add_argument(
        "-i",
        "--inhibitors-file",
        type=str,
        required=True,
        dest="inhibitors_file",
        help="The name of the input inhibitors file",
    )
    parser.add_argument(
        "-u",
        "--disease-up",
        type=str,
        required=True,
        dest="disease_up",
        help="The name of the diease up-regulated file",
    )
    parser.add_argument(
        "-d",
        "--disease-down",
        type=str,
        required=True,
        dest="disease_down",
        help="The name of the disease down-regulated file",
    )
    parser.add_argument(
        "-e",
        "--data-directory",
        type=str,
        required=True,
        dest="data_directory",
        help="The data directory",
    )
    parser.add_argument(
        "-r",
        "--raw-drug-file",
        type=str,
        required=True,
        dest="raw_drug_file",
        help="The name of the raw drug file",
    )
    args = parser.parse_args()
    tissue_spec_model_file = args.tissue_file
    inhibitors_file = args.inhibitors_file
    disease_up_file = args.disease_up
    disease_down_file = args.disease_down
    data_directory = args.data_directory
    drug_raw_file = args.raw_drug_file

    # TODO: Remove this after verifying argparse works
    # try:
    #     opts, args = getopt.getopt(
    #         argv,
    #         "ht:i:u:d:f:r:",
    #         ["tfile=", "ifile=", "upfile=", "downfile=", "output_dir=", "drugfile="],
    #     )
    # except getopt.GetoptError:
    #     print(
    #         "python3 knock_out_simulation.py -t <tissue_model_file> -i <inhibitor_file> -u <up_reg_file> -d <down_reg_file> -f <output_dir> -r <repurpose drug raw file>"
    #     )
    #     sys.exit(2)
    # for opt, arg in opts:
    #     if opt in ("-t", "--tfile"):
    #         tissue_spec_model_file = arg
    #     elif opt in ("-i", "--ifile"):
    #         inhibitors_file = arg
    #     elif opt in ("-u", "--upfile"):
    #         Disease_Up_file = arg
    #     elif opt in ("-d", "--downfile"):
    #         Disease_Dn_file = arg
    #     elif opt in ("-f", "--output_dir"):
    #         datadir = arg
    #     elif opt in ("-r", "--drugfile"):
    #         drug_raw_file = arg
    # datadir = os.path.join(configs.datadir, folder)

    print(f"Output directory: '{data_directory}'")
    print(f"Tissue Specific Model file is at: {tissue_spec_model_file}")
    print(f"Tissue specific inhibitors is at: {inhibitors_file}")

    # preprocess repurposing hub data
    drug_csv_file = "Repurposing_Hub_Preproc.csv"
    drugRawFile = os.path.join(configs.datadir, drug_raw_file)
    drugFile = os.path.join(data_directory, drug_csv_file)
    if not os.path.isfile(drugFile):
        drugDB = repurposing_hub_preproc(drugRawFile)
        drugDB.to_csv(drugFile, index=False)
    else:
        drugDB = pd.read_csv(drugFile)

    # Knock Out Simulation
    (
        model,
        geneInd2genes,
        HasEffects_Gene,
        fluxsolution,
        fluxSolutionRatios,
        fluxSolutionDiffs,
    ) = knock_out_simulation(
        datadir=configs.datadir,
        model_file=tissue_spec_model_file,
        inhibitors=inhibitors_file,
        drugDB=drugDB,
    )
    Gene_Pairs_down = create_gene_pairs(
        configs.datadir,
        model,
        geneInd2genes,
        fluxsolution,
        fluxSolutionRatios,
        fluxSolutionDiffs,
        HasEffects_Gene,
        Disease_Down=disease_down_file,
    )
    Gene_Pairs_down.to_csv(
        os.path.join(data_directory, "Gene_Pairs_Inhi_Fratio_DOWN.txt"), index=False
    )
    Gene_Pairs_up = create_gene_pairs(
        configs.datadir,
        model,
        geneInd2genes,
        fluxsolution,
        fluxSolutionRatios,
        fluxSolutionDiffs,
        HasEffects_Gene,
        Disease_Down=disease_up_file,
    )
    Gene_Pairs_up.to_csv(
        os.path.join(data_directory, "Gene_Pairs_Inhi_Fratio_UP.txt"), index=False
    )
    # print(geneInd2genes)
    # print(fluxSolutionRatios)
    # print(HasEffects_Gene)
    # Gene_Pairs_down = load_Inhi_Fratio(os.path.join(datadir,'Gene_Pairs_Inhi_Fratio_DOWN.txt'))
    # Gene_Pairs_up = load_Inhi_Fratio(os.path.join(datadir,'Gene_Pairs_Inhi_Fratio_UP.txt'))
    d_score_down = score_gene_pairs(
        Gene_Pairs_down, os.path.join(data_directory, "d_score_DOWN.csv")
    )
    # d_score_down = score_gene_pairs_diff(Gene_Pairs_down, 'd_score_DOWN.csv')
    d_score_up = score_gene_pairs(
        Gene_Pairs_up, os.path.join(data_directory, "d_score_UP.csv")
    )
    # d_score_up = score_gene_pairs_diff(Gene_Pairs_up, 'd_score_UP.csv')
    PES = (d_score_up - d_score_down).sort_values(by="score", ascending=False)
    PES.to_csv(os.path.join(data_directory, "d_score.csv"))
    print(d_score_down)
    print(d_score_up)
    print(PES)
    PES.reset_index(drop=False, inplace=True)
    # PES = pd.read_csv(os.path.join(datadir,'d_score.csv'))

    # last step: output drugs based on d score
    drug_score = drug_repurposing(drugDB, PES)
    drugScoreFile = os.path.join(data_directory, "drug_score.csv")
    drug_score.to_csv(drugScoreFile, index=False)
    print("Gene D score mapped to repurposing drugs saved to\n{}".format(drugScoreFile))


if __name__ == "__main__":
    main(sys.argv)
