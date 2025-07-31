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
from statistics import mean
from certifi.core import where
from jupyter_server.auth import passwd
from networkx.classes import is_empty
from numpy.ma.extras import unique
from rpy2.robjects.lib.dplyr import setdiff
from scipy.sparse import csr_matrix
from xlrd.formula import num2strg

from reconstruction_methods.ftinit.ftinit_fill_gaps_for_all_tasks import ftinit_fill_gaps_for_all_tasks
from reconstruction_methods.ftinit.ftinit_internal_alg import ftinit_internal_alg
from reconstruction_methods.ftinit.getinitsteps import get_initstep
from reconstruction_methods.ftinit.group_rxn_scores import group_rxn_scores
from reconstruction_methods.ftinit.remove_low_score_genes import remove_low_score_genes
from reconstruction_methods.ftinit.reverse_rxns import reverse_rxns
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
    # Previously, this was handled by giving a bonus for secreting those metabolites, but that doesn't work since
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
    origRxnScores = score_complex_model(prepData["refModel"], hpaData, transcrData,tissue,celltype)
    origRxnScores = np.where((origRxnScores > 0.1) & (origRxnScores <= 0), -0.1, origRxnScores)
    origRxnScores = np.where((origRxnScores < 0.1) & (origRxnScores > 0), 0.1, origRxnScores)

    # Initialize boolean array and flux array
    rxnsTurnedOn = np.zeros(len(prepData["minModel"]["rxns"]), dtype=bool)
    fluxes = np.zeros(len(prepData["minModel"]["rxns"]))
    
    rxnsToIgnoreLastStep = [[1],[1],[1],[1],[1],[1],[1],[1]]

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

        if "MetsToIgnore" in stp and stp["MetsToIgnore"] is not None:
            simple_mets_info = stp["MetsToIgnore"].get("simpleMets", {})
            if simple_mets_info:
                mets_to_rem = np.isin(mm.metNames, simple_mets_info["mets"])
                comps_to_keep = [i for i, c in enumerate(mm.comps) if c in simple_mets_info["compsToKeep"]]
                comps_to_keep_set = set(comps_to_keep)
                # Filter our metabolites that belong to compartments we want to keep
                mets_to_rem = np.logical_and(
                    mets_to_rem,
                    ~np.isin(mm.metComps, comps_to_keep)
                )
                # Set corresponding rows in the stoichiometric matrix to zero
                mm.S[mets_to_rem, :] = 0

        # Set up the reaction scores and essential rxns
        rxnToIgnore = getRxnsFromPattern(stp["RxnsToIgnoreMask"], prepData)
        rxnScores = group_rxn_scores(prepData["minModel"], origRxnScores, prepData["refModel"].rxns, prepData["groupIds"], rxnsToIgnore)

        essentialRxns = prepData["essentialRxns"]
        toRev = np.zeros(len(mm.rxns), dtype=bool)

        # Handle results from previous steps
        if stp["HowToUsePrevResults"] == "exclude":
            rxnScores[rxnsTurnedOn] = 0
        elif stp["HowToUsePrevResults"] == "essential":
            rev = np.array(mm.rev) == 1
            toRev = rxnsTurnedOn & rev & (fluxes < 0)
            mm = reverse_rxns(mm, [mm.rxns[i] for i, rev_it in enumerate(toRev) if rev_it])

            # Then make reactions irreversible
            mm.rev[rxnsTurnedOn] = 0
            mm.lb[rxnsTurnedOn] = 0
            essentialRxns = list(set(prepData["essentialRxns"] + [mm.rxns[i] for i, on in enumerate(rxnsTurnedOn) if on]))

        # Prep to run MILP iterations
        minGap = 1
        first = True
        success = False
        fullMipRes = {}
        for rn, params in enumerate(stp["MILPParams"]):
            if "MIPGap" not in params:
                params["MIPGap"] = 0.0004
            if "TimeLimit" not in params:
                params["TimeLimit"] = 5000

            if not first:
                # There is sometimes a problem with that the objective function becomes close to zero, which leads to that
                # a small percentage of that (which is the MIPGap sent in) is very small and the MILP hence takes a lot of time to finish.
                # We also therefore use an absolute MIP gap, converted to a percentage using the last value of the objective function.
                last_abs_obj = abs(lastObjVal) if abs(lastObjVal) > 1e-12 else 1e-12 # Avoid divide-by-zero
                params["MIPGap"] = min(max(params["MIPGap"], stp["AbsMIPGaps"][rn] / last_abs_obj), 1)
                params["seed"] = 1234 # use another seed, may work better
                if mipGap <= params["MIPGap"]:
                    success = True
                    break # we are done - this will not happen the first time
                else:
                    print(f"MipGap too high, trying with a different run. MipGap = {mipGap} New MIPGap Limit = {params[MIPGap]}")
            first = False

            try:
                # The prodweight for metabolomics is currently set to 5 - 0.5 was default in the old version, which I deemed very small?
                # There could be a need to specify this somewhere in the call at some point.
                # This value has not been evaluated, but is assumed in the test cases - if changed, update the test case
                startVals = fullMipRes.get("full", None)
                result = ftinit_internal_alg(mm, rxnScores, metData, essentialRxns, 5, stp["AllowMetSecr"], stp["PosRevOff"], params, startVals, fluzes, verbose)

                deletedRxnsInINTI1, metProduction, fullMipRes, rxnsTurnedOn1, fluxes1 = result

                # Reverse flux direction for reactions we reversed
                fluxes1[toRev] = -fluxes1[toRev]

                mipGap = fullMipRes["mipGap"]
                lastObjVal = fullMipRes["obj"]

            except Exception as e:
                mipGap = float("inf")
                lastObjVal = float("inf") # we need to set something here, Ing leads to that this doesn't come into play

            success = mipGap <= params["MIPGap"]

        if not success:
            raise RuntimeError(f"Failed to find good enough solution within the time frame. MIPGap: {mipGap}")

        # Save rxnsTurnedOn and their fluxes for the next step
        rxnsTurnedOn = rxnsTurnedOn | np.array(rxnsTurnedOn1, dtype=bool)

        # The fluxes are a bit tricky - what if they change direction between the steps?
        #  The fluxes are used to determine the direction in which reactions are forced on (to simplify the problem it is
        # good if they are unidirectional).
        # We use the following strategy:
        # 1. Use the fluxes from the most recent step.
        # 2. If any flux is very low there. (i.e. basically zero), use the flux from the previous steps
        # This could in theory cause problems, but seems to work well practically
        fluxesOld = fluxes
        fluxes = np.array(fluxes1)
        # make sure that all reactions that are on actually has a flux - otherwise
        # things could go bad, since the flux will be set to essential in a random direction
        # This sometimes happens for rxns with negative score - let's just accept that.
        low_flux = np.abs(fluxes) < 1e-7
        fluxes[low_flux] = fluxesOld[low_flux]

        # get the essential reactions
        essential = np.isin(prepData["minModel"].rxns, prepData["essentialRxns"])
        # So, we only add reactions where the linearly merged scores are zero for all linearly dependent reactions
        # (this cannot happen by chance, taken care of in the function groupRxnScores)
        rxnsToIgn = rxnScores == 0
        deletedRxnsInINITSel = ~np.isin(rxnsTurnedOn, rxnsToIgn, essential)
        deletedRxnsInINIT = prepData["minModel"].rxns[deletedRxnsInINITSel]

        # Here we need to figure out which original reactions (before the linear merge) that were removed. These are all
        # reactions with the same group ids as the removed reactions.

        """
        rewrite the code below
        """
        groupIdsRemoved = prepData["groupIds"] in (prepData["refModel"].rxns, deletedRxnsInINIT)
        groupIdsRemoved = where(groupIdsRemoved != 0)
        rxnsToRem = union(prepData["refModel"].rxns[groupIdsRemoved], deletedRxnsInINIT)

        initModel = removeReactions(prepData["refModel"], rxnsToRem,False,True)

        # remove metabolites separately to avoid removing those needed for tasks
        unusedMets = initModel.mets(all(initModel.S == 0,2))
        initModel = removeMets(initModel, setdiff(unusedMets, prepDatap["essentialMetsForTasks"]))

        # if printReport:
        #     printScores(initModel, 'INIT model statistics', hpaData, transcrData, tissue, celltype)
        #     printScores(removeReactions(cModel,setdiff(cModel.rxnsm rxnsToRem), True, True), 'Reactions deleted by INIT', hpaData, transcrData, tissue, celltype)

        # The full model has exchange reactions in it. ftINITFillGapsForAllTasks calls ftINITFillGaps, which automatically removes
        # exchange metabolites (because it assumes that the reactions are constrained when appropriate). In this case the uptakes/outputs
        # are retrieved from the task sheet instead. To prevent exchange reactions being used to fill gaps, they are deleted
        # from the reference model here.
        initModel["id"] = 'INITModel'

        # If gaps in the model should be filled using a task list
        if prepData["tasksStruct"]:
            # Remove exchange reactions and reactions already included in the INIT model
            # We changed strategy and instead include all rxns except the exchane rxns in the ref model
            # But we do keep the exchange rxns that are essential.
            # Let's test to remove all, that should work

            # At this stage the model is fully connected and most of the gebes with good scores should have been included. The final gap-filling
            # should take the scores of the genes into account, so that "rather bad" reactions are preferred to "very bad" reactions.
            # However, reactions with positive scores will be included even if they are not connected in the current formation.
            # Therefore, such reactions will have to be assigned a small negative score instead.
            exchRxns = getExchangeRxns(prepData["refModel"])
            refModelNoExc = removeReactions(prepData.refModelWithBM,exchRxns,False,True)
            exchRxns = getExchangeRxns(initModel)
            initModelNoExc = removeReactions(closeModel(initModel),exchRxns,False,True)

            if useScoreForTasks:
                # map the rxn scores to the model without exchange rxns
                pass
            else:
                outModel, addedRxnMat = ftinit_fill_gaps_for_all_tasks(initModelNoExc,refModelNoExc,[],True,[],prepData["taskStruct"],paramsFT,verbose)
            addedRxnsForTasks = refModelNoExc.rxns[addedRxnMat,2]
        else:
            outModel = initModel
            addedRxnMat = []
            addedRxnsForTasks = []


        # The model can now perform all the tasks defined in the task list.
        model = outModel

        # At this stage the model will contain some exchange reactions but probably not all (and maybe zero). This can be inconvenient,
        # so all exchange reactions from the reference model are added, except for those which involve metabolites that are not in the model.

        # Start from the original model, and just remove the reactions that are no longer there (and keep exchange rxns).
        # from the problem is not complete, it doesn't have GRPs etc.
        # The logic below is a bit complicated. We identify the reactions that should be removed from the full model as reactions
        # that have been removed in the init model except the ones that were added back. In addition, we make reactions that
        # have been removed in the init model. except the ones that were added back. In addition, we make sure that no exchange rxns
        # are removed - they can be removed in the init model if they were linearly merged with other reactions that were decided to be
        # removed from the model. We want to keep all exchange rxns to make syre the tasks can be performed also without manipulating the b vector
        # in the model (which is what is done in the gap-filling).

        exchRxns = getExchangeRxns(prepData["refModel"])
        deleteRxnsInINIT = setdiff(prepData["refModel"].rxns, union(union(initModel.rxns, addedRxnsForTasks), exchRxns))
        outModel = removeReactions(prepData["refModel"], deleteRxnsInINIT, True) # we skip removing the genes for now, I'm not sure it is desireable

        # If requested, attempt to remove negative-score genes from the model, depending on their role (isozyme or complex subunit)
        # in each grRule. See the "removeLowScoreGenes" function for more details, and to adjust any default parameters therein.

        if (removeGenes):
            ~,geneScores = score_complex_model(outModel, hpaData, transcrData,tissue,celltype)
            outModel = remove_low_score_genes(outModel,geneScores)

        model = outModel

        # This is for printing a summary of a model

        def printScores(model, name, hpaData, transcrData, tissue, celltype):
            a,b = score_complex_model(model, hpaData, transcrData,tissue,celltype)
            rxnS = mean(a)
            geneS = mean(b,'omitnan')
            print(f"{name}:\n")
            print(f"\t {num2strg(sum(model.rxns))} reactions, {num2strg(sum(model.genes))} genes\n")
            print(f"\tMean reaction score: {num2strg(rxnS)}\n")
            print(f"\tMean gene score: {num2strg(geneS)}\n")
            print(f"\tReactions with positive scores: {num2strg(100*sum(a>0)/sum(a!=0))}%\n")
            return rxnS, geneS

        def getRxnsFromPattern(rxnsToIgnorePattern, prepData):
            rxns_to_ignore = np.zeros(len(prepData["toIgnoreExch"]), dtype=bool)
            if rxns_to_ignore[0]:
                rxns_to_ignore |= prepData["toIgnoreExch"]
            if rxns_to_ignore[1]:
                rxns_to_ignore |= prepData["toIgnoreImportRxns"]
            if rxns_to_ignore[2]:
                rxns_to_ignore |= prepData["toIgnoreSimpleTransp"]
            if rxns_to_ignore[3]:
                rxns_to_ignore |= prepData["toIgnoreAdvTransp"]
            if rxns_to_ignore[4]:
                rxns_to_ignore |= prepData["toIgnoreSpont"]
            if rxns_to_ignore[5]:
                rxns_to_ignore |= prepData["toIgnoreS"]
            if rxns_to_ignore[6]:
                rxns_to_ignore |= prepData["toIgnoreCustomRxns"]
            if rxns_to_ignore[7]:
                rxns_to_ignore |= prepData["toIgnoreAllWithoutGPRs"]
            return rxns_to_ignore






