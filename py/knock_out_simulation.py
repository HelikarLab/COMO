#!/usr/bin/python3
import os
import sys
import time
import pandas as pd
import numpy as np
import cobra
from cobra.flux_analysis import (single_gene_deletion, single_reaction_deletion,
                                 double_gene_deletion, double_reaction_deletion)

from transcriptomic_gen import *
from proteomics_gen import *

def knock_out_simulation(datadir, model_file, inhibitors):

    if model_file[-4:] == '.xml':
        model = cobra.io.read_sbml_model(os.path.join(datadir, model_file))
    elif model_file[-4:] == '.mat':
        model = cobra.io.load_matlab_model(os.path.join(datadir, model_file))
    else:
        print("Unsupported File Format of Model: {}".format(model_file))
        return None

    DT_genes = pd.read_csv(os.path.join(datadir,inhibitors), header=None)
    DT_genes.rename(columns={0:'Gene ID'},inplace=True)
    DT_genes['Gene ID'] = DT_genes['Gene ID'].astype(str)

    geneInd2genes = [x.id for x in model.genes]
    print(len(geneInd2genes))
    geneInd2genes = set(geneInd2genes)
    print(len(geneInd2genes))
    DT_model = set(DT_genes['Gene ID'].tolist()).intersection(geneInd2genes)
    print(len(DT_model))
    DT_model = list(DT_model)

    model_opt = model.optimize().to_frame()

    fluxsolution = pd.DataFrame()
    for id in DT_model:
        gene = model.genes.get_by_id(id)
        gene.knock_out()
        opt_model = model.optimize().to_frame()
        fluxsolution[id]=opt_model['fluxes']

    # fluxsolution

    fluxSolutionRatios = fluxsolution.div(model_opt['fluxes'],axis=0)
    # fluxSolutionRatios
    fluxSolutionDiffs = fluxsolution.sub(model_opt['fluxes'],axis=0)
    # fluxSolutionDiffs

    HasEffects_Gene = []
    for id in DT_model:
        gene = model.genes.get_by_id(id)
        for rxn in gene.reactions:
            if fluxSolutionDiffs.at[rxn.id,id] != 0:
                HasEffects_Gene.append(id)
                break
    # HasEffects_Gene
    return model, geneInd2genes, fluxSolutionRatios, HasEffects_Gene


def create_gene_pairs(datadir, model, geneInd2genes, fluxSolutionRatios, HasEffects_Gene, RA_Down):
    RA_down = pd.read_csv(os.path.join(datadir,RA_Down))
    DAG_dis_genes = pd.DataFrame()
    DAG_dis_genes['Gene ID'] = RA_down.iloc[:,0].astype(str)
    # DAG_dis_genes

    DAG_dis_met_genes = set(DAG_dis_genes['Gene ID'].tolist()).intersection(geneInd2genes)
    # DAG_dis_met_genes

    DAG_dis_met_rxnInd = []
    Gene_i = []
    for id in DAG_dis_met_genes:
        gene = model.genes.get_by_id(id)
        for rxn in gene.reactions:
            DAG_dis_met_rxnInd.append(rxn.id)
            Gene_i.append(id)

    # DAG_dis_met_rxnInd
    Gene_df = pd.DataFrame(Gene_i,columns=['Gene IDs'],index=DAG_dis_met_rxnInd)
    # Gene_df

    DAG_rxn_fluxRatio = fluxSolutionRatios.loc[DAG_dis_met_rxnInd]
    # DAG_rxn_fluxRatio

    gene_mat_out = []
    # Gene_i = DAG_dis_met_genes
    # Rind_i = DAG_dis_met_rxnInd
    for id in HasEffects_Gene:
        pegene = pd.DataFrame()
        pegene['Gene IDs'] = Gene_df['Gene IDs'].copy()
        pegene['rxn_fluxRatio'] = DAG_rxn_fluxRatio[id].copy()
        pegene['Gene'] = id
        pegene.dropna(axis=0,subset=['rxn_fluxRatio'],inplace=True)
        pegene.index.name='reaction'
        pegene.reset_index(drop=False, inplace=True)
        gene_mat_out.append(pegene)

    Gene_Pairs = pd.concat(gene_mat_out, ignore_index=True)
    return Gene_Pairs

def main(argv):
    print(projectdir)
    datadir = os.path.join(projectdir,'data')
    print(datadir)
    model, geneInd2genes, fluxSolutionRatios, HasEffects_Gene = knock_out_simulation(datadir=datadir,
                                      model_file='Th1_Cell_SpecificModel4manuscript.xml',
                                      inhibitors='Th1_inhibitors_Entrez.txt')
    Gene_Pairs_down = create_gene_pairs(datadir,
                                   model,
                                   geneInd2genes,
                                   fluxSolutionRatios,
                                   HasEffects_Gene,
                                   RA_Down='RA_DOWN.txt')
    Gene_Pairs_down.to_csv(os.path.join(datadir,'Gene_Pairs_Inhi_Fratio_DOWN.txt'),index=False)
    Gene_Pairs_up = create_gene_pairs(datadir,
                                   model,
                                   geneInd2genes,
                                   fluxSolutionRatios,
                                   HasEffects_Gene,
                                   RA_Down='RA_DOWN.txt')
    Gene_Pairs_up.to_csv(os.path.join(datadir,'Gene_Pairs_Inhi_Fratio_DOWN.txt'),index=False)
    print(geneInd2genes)
    print(fluxSolutionRatios)
    print(HasEffects_Gene)

if __name__ == "__main__":
   main(sys.argv[1:])
