# This is the Python version of ftINIT originally written in MATLAB.
# The tool was developed with the intention to develop genome-scale metabolic models from scRNA-seq  data.
# Source code:  https://github.com/SysBioChalmers/RAVEN/blob/cd82b5c6ba26a88bed2a9ca179c682a6434fd541/INIT/ftINIT.m

# ftINIT
#   Main function for generates a model using the ftINIT algorithm, based
#   on proteomics and/or transcriptomics and/or metabolomics and/or metabolic
#   tasks. The algorithm is not designed for running with metabolomics only.
#   The function prepINITModel needs to be run first for the template
#   model (such as the generic Human-GEM), but only need to be run once. This
#   function precalculates things independent of the omics to speed up the model
#   generation process, and outputs the prepData, which is input to this function.
#
#   prepData            The prepdata for the model.
#   tissue              tissue to score for. Should exist in either
#                       hpaData.tissues or transcrData.tissues
#   celltype            cell type to score for. Should exist in either
#                       hpaData.celltypes or transcrData.celltypes for this
#                       tissue (optional, default is to use the max values
#                       among all the cell types for the tissue.
#   hpaData             HPA data structure from parseHPA (optional if transcrData
#                       is supplied, default [])
#   transcrData         gene expression data structure (optional if hpaData is
#                       supplied, default []). Used to be called arrayData.
#       genes           cell array with the unique gene names
#       tissues         cell array with the tissue names. The list may not
#                       be unique, as there can be multiple cell types per
#                       tissue.
#       celltypes       cell array with the cell type names for each tissue
#       levels          GENESxTISSUES array with the expression level for
#                       each gene in each tissue/celltype. NaN should be
#                       used when no measurement was performed
#       threshold       a single value or a vector of gene expression
#                       thresholds, above which genes are considered to be
#                       "expressed". default = 1(optional, by default, the mean expression
#                       levels of each gene across all tissues in transcrData
#                       will be used as the threshold values)
#       singleCells     binary value selecting whether to use the
#                       single-cell algorithm to identify expressed genes.
#                       If used, specify cell subpopulations in CELLTYPES
#                       (optional, default [])
#       plotResults     true if single cell probability distributions
#                       should be plotted (optional, default = false)
#   metabolomicsData    cell array with metabolite names that the model
#                       should produce (optional, default [])
#   INITSteps           Specifies the steps in the algorithm. For more info,
#                       see INITStepDesc and getINITSteps.
#                       (optional, default getINITSteps(), which is the standard ftINIT).
#   removeGenes         if true, low-abundance genes will be removed from
#                       grRules, unless they are the only gene associated
#                       with a reaction, or a subunit of an enzyme complex
#                       (see "removeLowScoreGenes" function for details).
#                       If false, grRules will not be modified; however,
#                       genes that were associated only with removed
#                       reactions will not be present in the final model.
#                       (optional, default true).
#   useScoresForTasks   true if the calculated reaction scored should be
#                       used as weights when fitting to tasks (optional, default
#                       true)
#   paramsFT            parameter structure as used by getMILPParams. This
#                       is for the fitTasks step. For the INIT algorithm,
#                       see params (optional, default [])
#   verbose             if true, the MILP progression will be shown.
#                       (optional, default false)
#
#   model                   the resulting model structure
#   metProduction           array that indicates which of the
#                           metabolites in metabolomicsData that could be
#                           produced. Note that this is before the
#                           gap-filling process to enable defined tasks. To
#                           see which metabolites that can be produced in
#                           the final model, use canProduce.
#                           -2: metabolite name not found in model
#                           -1: metabolite found, but it could not be produced
#                           1: metabolite could be produced
#   addedRxnsForTasks       cell array of the reactions which were added in
#                           order to perform the tasks
#   deletedRxnsInINIT       cell array of reactions deleted because they
#                           could not carry flux (INIT requires a
#                           functional input model)
#   fullMipRes              The solver results from the last MILP step run
#
#   This is the main function for automatic reconstruction of models based
#   on the ftINIT algorithm ().
#
#   NOTE: Exchange metabolites should normally not be removed from the model
#   when using this approach, since checkTasks/fitTasks rely on putting specific
#   constraints for each task. The INIT algorithm will remove exchange metabolites
#   if any are present. Use importModel(file,false) to import a model with
#   exchange metabolites remaining.
#
# Usage: [model, metProduction, addedRxnsForTasks, deletedRxnsInINIT, ...
#               fullMipRes] = ...
#               ftINIT(prepData, tissue, celltype, hpaData, transcrData, ...
#               metabolomicsData, INITSteps, removeGenes, useScoresForTasks, ...
#               paramsFT)

###----------------------------------------------------CODE----------------------------------------------------###

from typing import List, Tuple, Optional
import numpy as np
import scipy.sparse as sp
from cobra import Model
from cobra.util import solver

from reconstruction_methods.ftinit.ftinit_fill_gaps_for_all_tasks import ftinit_fill_gaps_for_all_tasks
from reconstruction_methods.ftinit.ftinit_internal_alg import ftinit_internal_alg
from reconstruction_methods.ftinit.getinitsteps import get_initstep
from reconstruction_methods.ftinit.group_rxn_scores import group_rxn_scores
from reconstruction_methods.ftinit.remove_low_score_genes import remove_low_score_genes
from reconstruction_methods.ftinit.reverse_rxns import reverse_rxns
from reconstruction_methods.ftinit.score_complex_model import score_complex_model

def run_ftinit(prep_data, tissue: str, celltype: Optional[str]=None, hpa_data=None, transcr_data=None, metabolomics_data=None,
               INIT_steps=None, remove_genes: bool=True, use_score_for_tasks: bool=True, params_ft: Optional[dict]=None, verbose: bool=False)-> Tuple[Model, np.ndarray,List[str],List[str],dict]:

    if INIT_steps is None:
        INIT_steps = get_initstep([],'1+1')

    if metabolomics_data:
        if len(set(m.upper() for m in metabolomics_data)) != len(metabolomics_data):
            raise ValueError("Metabolomics contains the same metabolite multiple times")

        met_data = sp.lil_matrix((len(metabolomics_data), len(prepData['min_model'].reactions)))

        ref_model = prepData['ref_model']
        for i, met_name in enumerate(metabolomics_data):
            met_sel = [met.upper() == met_name.upper() for met in ref_model.metabolites.list_attr("name")]
            S = ref_model.solver.matrix
            prod_rxns_sel = np.any(S[met_sel,:] > 0, axis=0) | \
                            (np.any(S[met_sel,:] < 0, axis=0) & ref_model.reversible)
            # convert rxns from ref_model to min_model
            group_ids = prep_data['group_ids']
            min_rxns = prep_data['min_model'].reactions.list_attr("id")
            ref_rxns = ref_model.reactions.list_attr("id")
            ia_ib = [{i, ref_rxns.index(r.id)) for i, r in enumerate(prep_data['min_model'].reactions) if r.id in ref_rxns]
            grp_ids_merged = [np.nan] * len(min_rxns)
            for i,j in ia_ib:
                grp_ids_merged[i] = group_ids[j]

                group_ids_pos = list(set(group_ids[i] for i, val in enumerate(prod_rxns_sel) if val))
                group_ids_pos = [g for g in group_ids_pos if g != 0]

                direct_match = [r.id in ref_model.reactions for r in prep_data['min_model'].reactions]
                met_data[i,:] = [(g in group_ids_pos) or dm for g, dm in zip(grp_ids_merged,direct_match)]
            met_data = sp.csr_matrix(met_data)
        else:
            met_data = None

        orig_rxn_scores, _ = score_complex_model(prep_data['ref_model'], hpa_data,transcr_data,tissue,celltype)
        orig_rxn_scores = np.where((orig_rxn_scores > 0.1) & (orig_rxn_scores <= 0), -0.1, orig_rxn_scores)
        orig_rxn_scores = np.where((orig_rxn_scores < 0.1) & (orig_rxn_scores > 0), 0.1, orig_rxn_scores)

        rxn_turned_on = np.zeros(len(prep_data['min_model'].reactions), dtype=bool)
        fluxes = np.full(len(prep_data['min_model'].reactions),0.1)
        rxns_to_ignore_last_step = [1] * 8

        full_mip_res = {}
        for i, step in enumerate(INIT_steps):
            print(f"ftINIT: Running step {i+1} of {len(INIT_steps)}")

            if any(a-b < 0 for a,b in zip(rxns_to_ignore_last_step, step['rxns_to_ignore'])):
                raise ValueError("Rxns_to_ignore_mask must cover prior ignored rxns.")

            rxns_to_ignore_last_step = step['rxns_to_ignore_mask']

            mm = prep_data['min_model'].copy()

            if step.get('mets_to_ignore') and step['mets_to_ignore'].get('simple_mets'):
                mets_to_erm = [m.name in step['mets_to_ignore']['simple_mets']['mets'] for m in mm.metabolites]
                somps_to_keep = [mm.compartments.index(c)]

