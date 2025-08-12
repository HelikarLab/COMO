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
from reconstruction_methods.ftinit.remove_reactions import remove_reactions
from reconstruction_methods.ftinit.get_exchange_rxns import get_exchange_rxns
from reconstruction_methods.ftinit.remove_mets import remove_mets

def run_ftinit(prep_data, tissue: str, celltype: Optional[str]=None, hpa_data=None, transcr_data=None, metabolomics_data=None,
               INIT_steps=None, remove_genes: bool=True, use_score_for_tasks: bool=True, params_ft: Optional[dict]=None, verbose: bool=False)-> Tuple[Model, np.ndarray,List[str],List[str],dict]:
    prep_data = prep_data.copy()
    use_score_for_tasks = use_score_for_tasks
    if INIT_steps is None:
        INIT_steps = get_initstep([],'1+1')

    if metabolomics_data:
        if len(set(m.upper() for m in metabolomics_data)) != len(metabolomics_data):
            raise ValueError("Metabolomics contains the same metabolite multiple times")

        met_data = sp.lil_matrix((len(metabolomics_data), len(prep_data['min_model'].reactions)))

        ref_model = prep_data['ref_model']
        for i, met_name in enumerate(metabolomics_data):
            met_sel = [met.upper() == met_name.upper() for met in ref_model.metabolites.list_attr("name")]
            S = ref_model.solver.matrix
            prod_rxns_sel = np.any(S[met_sel,:] > 0, axis=0) | \
                            (np.any(S[met_sel,:] < 0, axis=0) & ref_model.reversible)
            # convert rxns from ref_model to min_model
            group_ids = prep_data['group_ids']
            min_rxns = prep_data['min_model'].reactions.list_attr("id")
            ref_rxns = ref_model.reactions.list_attr("id")
            ia_ib = [(i, ref_rxns.index(r.id)) for i, r in enumerate(prep_data['min_model'].reactions) if r.id in ref_rxns]
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

        rxns_turned_on = np.zeros(len(prep_data['min_model'].reactions), dtype=bool)
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
                mets_to_rem = [m.name in step['mets_to_ignore']['simple_mets']['mets'] for m in mm.metabolites]
                comps_to_keep = [mm.compartments.index(c) for c in step['mets_to_ignore']['simple_mets']['comps_to_keep']]
                met_comp_mask = [m.compartment not in comps_to_keep for m in mm.metabolites]
                for i, remove in enumerate(np.logical_and(mets_to_rem, met_comp_mask)):
                    if remove:
                        mm.solver.matrix[i,:] = 0
            rxns_to_ignore = get_rxns_from_pattern(step['rxns_to_ignore_mask'], prep_data)
            rxn_scores = group_rxn_scores(mm, orig_rxn_scores, prep_data['ref_model'].reactions.list_attr("id"), prep_data['group_ids'], rxns_to_ignore)

            essential_rxns = prep_data['essential_rxns'][:]
            to_rev = np.zeros(len(mm.reactions),dtype=bool)

            if step['how_to_use_results'] == 'exclude':
                rxn_scores[rxns_turned_on] = 0
            elif step['how_to_use_prev_results'] == 'essential':
                rev = [rxn.reversibility for rxn in mm.reactions]
                to_rev = np.logical_and(rxns_turned_on, rev & (fluxes < 0))
                mm = reverse_rxns(mm, [r.id for i, r in enumerate(mm.reations) if to_rev[i]])
                for i,r in enumerate(mm.reactions):
                    if rxns_turned_on[i]:
                        r.lower_bound = 0
                        r.reversibility = False
                essential_rxns += [r.id for i, r in enumerate(mm.reactions) if rxns_turned_on[i]]

            success = False
            first = True
            mip_gap = 1

            for params in step['milp_params']:
                params.setdefault('milp_gap', 0.0004)
                params.setdefault('time_limit', 5000)
                if not first:
                    params['milp_gap'] = min(max(params['milp_gap'], step['abs_mip_gap'][0]/abs(last_obj_val)),1)
                    params['seed'] = 1234
                    if mip_gap <= params['milp_gap']:
                        success = True
                        break
                first = False
                start_vals = full_mip_res.get('full', None)
                try:
                    deleted_rxns_in_INIT1, met_production, full_mip_res, rxns_turned_on1,fluxes1 = ftinit_internal_alg(mm,rxn_scores, met_data, essential_rxns, 5, step['allow_met_secr'],step['pos_rev_off'], params,start_vals,verbose)
                    fluxes1[to_rev] = -fluxes1[to_rev]
                    mip_gap = full_mip_res['obj']
                except Exception:
                    mip_gap = float('inf')
                    last_obj_val = float('inf')

                success = mip_gap <= params['milp_gap']

            if not success:
                raise RuntimeError(f"Failed to find good enough solution. MIPGap: {mip_gap}")

            rxns_turned_on = np.logical_or(rxns_turned_on, rxns_turned_on1)
            fluxes = np.where(np.abs(fluxes1) > 1e-6, fluxes1, fluxes)

        essential = [r.id in prep_data['essential_rxns'] for r in prep_data['min_model'].reactions]
        rxns_to_ign = rxn_scores == 0
        deleted_sel = ~(rxns_turned_on | rxns_to_ign | essential)
        deleted_rxns_in_INIT = [r.id for i, r in enumerate(prep_data['min_model'].reactions) if deleted_sel[i]]

        group_ids_removed =[prep_data['group_ids'][prep_data['ref_model'].reactions.index(r)]
                            for r in deleted_rxns_in_INIT if prep_data['group_ids'][prep_data['ref_model'].reactions.index(r)] != 0]
        rxns_to_rem = set([r.id for i, r in enumerate(prep_data['ref_model'].reactions) if prep_data['group_ids'][i] in group_ids_removed] + deleted_rxns_in_INIT)

        init_model = remove_reactions(prep_data['ref_model'], rxns_to_rem, False, True,True)
        unused_mets = [m.id for m in init_model.metabolites if all(v == 0 for v in init_model.metabolite_coefficients(m).values())]
        init_model = remove_mets(init_model, set(unused_mets) - set(prep_data['essential_mets_for_tasks']))

        init_model.id = 'init_model'

        if prep_data.get('task_struct'):
            exch_rxns = get_exchange_rxns(prep_data['ref_model'])
            ref_model_no_exc = remove_reactions(prep_data['ref_model_with_BM'], exch_rxns, remove_genes=False, remove_unused=True)
            exch_rxns = get_exchange_rxns(init_model)
            init_model_no_exc = remove_reactions(close_model(init_model), exch_rxns, remove_genes=False, remove_unused=True)

            if use_score_for_tasks:
                ref_rxns = prep_data['ref_model'].reactions.list_attr("id")
                ref_rxns_no_exc = ref_model_no_exc.reactions.list_attr("id")
                rxn_scores_2nd = np.full(len(ref_rxns_no_exc), np.nan)
                for i, r in enumerate(ref_rxns_no_exc):
                    if r in ref_rxns:
                        rxn_scores_2nd[i] = orig_rxn_scores[ref_rxns.index(r)]
                out_model, added_rxn_mat = ftinit_fill_gaps_for_all_tasks(init_model_no_exc, ref_model_no_exc, None, True, np.minimum(rxn_scores_2nd, -0.1), prep_data['task_struct'], params_ft, verbose)
            else:
                out_model, added_rxn_mat = ftinit_fill_gaps_for_all_tasks(init_model_no_exc, ref_model_no_exc, None, True,None,prep_data['task_strct'], params_ft,verbose)
            added_rxns_for_tasks = [r.id for i, r in enumerate(ref_model_no_exc.reactions) if any(added_rxn_mat[i,:])]
        else:
            out_model = init_model
            added_rxns_for_tasks = []

        exch_rxns = get_exchange_rxns(prep_data['ref_model'])
        deleted_rxns_in_INIT = list(set(prep_data['ref_model'].reactions.list_attr("id")) - set(out_model.reactions.list_attr("id")) - set(added_rxns_for_tasks) - set(exch_rxns))

        if remove_genes:
            _, gene_scores = score_complex_model(out_model, hpa_data, transcr_data, tissue, celltype)
            out_model = remove_low_score_genes(out_model, gene_scores)
        return out_model, met_production, added_rxns_for_tasks, deleted_rxns_in_INIT, full_mip_res

    def get_rxns_from_pattern(pattern: List[int], prep_data) -> np.ndarray:
        flags = [
            'to_ignore_exch', 'to_ignore_import_rxns', 'to_ignore_simple_transp',
            'to_ignore_adv_transp', 'to_ignore_spont', 'to_ignore_s',
            'to_ignore_custom_rxns', 'to_ignore_all_without_grps'
        ]
        return np.any([pattern[i] and prep_data[flags[i]] for i in range(8)], axis=0)