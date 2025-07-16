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

import inspect
import numpy as np
from jupyter_server.auth import passwd
from networkx.classes import is_empty
from numpy.ma.extras import unique
from scipy.sparse import csr_matrix
from xlrd.formula import num2strg

from reconstruction_methods.ftinit.getinitsteps import get_initstep
from reconstruction_methods.ftinit.score_complex_model import score_complex_model


def run_ftinit(prepData, tissue, celltype, hpaData, transcrData, metabolomicsData, INITSteps, removeGenes, useScoreForTasks, paramsFT, verbose):
    sig = inspect.signature(run_ftinit)
    num_args = len(sig.parameters)

    if num_args < 5:
        transcrData = []
    if num_args < 6:
        metabolomicsData = []
    if num_args < 7 or is_empty(INITSteps) :
        INITSteps = get_initstep([],'1+1')
    if num_args < 8 or is_empty(removeGenes):
        removeGenes = True
    if num_args < 9 or is_empty(useScoreForTasks):
        useScoreForTasks = True
    if num_args < 10:
        paramsFT = []
    if num_args < 11:
        verbose = False

    # Handle detected mets:
    # Prevviously, this was handled by giving a bonus for secreting those metabolites, but that doesn't work since
    # the metabolite secretion and uptake can be lost when we merge linearly dependent reaction.
    # Instead, we need to figure out which reactions either produce or take up the mets.
    # We then give a bonus if any of them carry flux.
    # To simplify things, we focus on reactions that produce the metabolite (since there must be once such reaction).
    # It is still a bit complicated though. In this step, we focus on identifying producer reactions. We further reason
    # that the direction doesn't matter - we can force one of these reactions in any direction - if it becomes a
    # consumer,it will automatically force another producer on as well (otherwise we'll have a net consumption).

    if metabolomicsData:
        # Check for duplicates (case-insensitive)
        if len(set(m.upper() for m in metabolomicsData)) != len(metabolomicsData):
            raise ValueError("Metabolomics contains the same metabolite multiple times")

        num_mets = len(metabolomicsData)
        num_rxns = len(prepData["minModel"]["rxns"])
        metData = np.zeros((num_mets, num_rxns), dtype=bool)

        for i, met_name in enumerate(metabolomicsData):
            # Match metabolites in a case-insensitive way
            metSel = np.array([name.upper()== met_name.upper() for name in prepData["refModel"]["metNames"]])

            S = prepData["refModel"]["S"]
            rev = np.array(prepData["refModel"]["rev"], dtype=bool)

            prodRxnsSel = np.any(S[metSel, :] > 0, axis=0) | (np.any(S[metSel, :] < 0, axis=0) & rev)

            # Map production reactions from refModel to minModel
            ref_rxns = prepData["refModel"]["rxns"]
            min_rxns = prepData["minModel"]["rxns"]
            group_ids = prepData["groupIds"]

            rxn_map = {rxn: idx for idx, rxn in enumerate(min_rxns)}
            ia = [rxn_map[rxn] for rxn in ref_rxns if rxn in rxn_map]
            ib = [ref_rxns.index(rxn) for rxn in ref_rxns if rxn in rxn_map]

            grpIdsMerged = np.full(len(min_rxns), np.nan)
            for m_idx, r_idx in zip(ia, ib):
                grpIdsMerged[m_idx] = group_ids

            groupIdsPos = np.unique(group_ids[prodRxnsSel])
            groupIdsPos = groupIdsPos[groupIdsPos != 0]

            posRxns = np.array(ref_rxns)[prodRxnsSel]
            directMatch = np.isin(min_rxns, posRxns)

            metData[i, :] = np.isin(grpIdsMerged, groupIdsPos) | directMatch

        metData = csr_matrix(metData)
    else:
        metData = None

    # Get rxn scores and adapt them to the minimized model
    """
    origRxnScores = score_complex_model(prepData.refModel, hpaData,transcrData,tissue,celltype)
    origRxnScores = 
    origRxnScores =
    
    rxnsTurnedOn = 
    fluxes = 
    
    rxnsToIgnoreLastStep = [[1],[1],[1],[1],[1],[1],[1],[1]]
    """

    # We assume that all essential rxns are irrev - this is taken care of in prepINITMode. We then use an initial flux "
    # from last run" of 0.1 for all reactions. This is used for knowing what flux should be forced through an essential rxn.

    #! have to figure out what 'ones' does in matlab
    fluxes = ones(len(prepData.minModel.rxns)) * 0.1

    for initStep in INITSteps:
        print(f'ftINIT: Running step {num2strg(initStep)}')
        stp = INITSteps(initStep)

        #! [rewrite] if any ((rxnsToIgnoreLastStep - stp.RxnsToIgnoreMask) < 0
            print('RxnsToIgnoreMask may not cover rxns not covered in previous steps, but the other way is fine.')
        rxnsToIgnoreLastStep = stp.RxnsToIgnoreMask

        mm = prepData.minModel

        if is_empty(stp.MetsToIgnore):
            pass
    # Set up the reaction scores and essential rxns
    """
    Code comes here
    """
    # Handle the results from previous steps ('ignore', 'exclude', 'essential')



