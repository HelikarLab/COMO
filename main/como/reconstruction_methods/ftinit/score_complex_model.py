# This is the Python version of scoreComplexModel (used in ftinit) originally written in MATLAB.
# Source code: https://github.com/SysBioChalmers/RAVEN/blob/cd82b5c6ba26a88bed2a9ca179c682a6434fd541/INIT/scoreComplexModel.m

# score_complex_model
# Scores the reactions snd genes in a model containing complex gene rules based on expression data from
# HPA and/or gene arrays. It is highly recommended that the model grRules are "cleaned" and verified by the
# "cleanModelGeneRules" function prior to scoring the model with scoreComplexModel.

## INPUT
# model: a model structure
# hpaData: HPA data structure from parseHPA (optional if arrayData is supplied, default [])
# arrayData: gene expression data structure (optional if hpaData is supplied, default [])
#   genes: cell array with the unique gene names
#   tissues: cell array with the tissue names. The list may not be unique, as there can be multiple cell types per tissue
#   celltypes: cell array with the cell type names for each tissue
#   levels: GENESxTISSUES array with the expression level for each gene in each tissue/celltype. NaN should be used when
#   no measurement was performed
#   threshold: a single value or a vector of gene expression thresholds, above which genes are considered to be "expressed".
#       (optional, by default, the mean expression levels of each gene across all tissues in arrayData will be used as the threshold values)
# tissue: tissue to score for. should exist in either hpaData.tissues or arrayData.tissues
# celltype: cell type to score for. should exist in either hpaData.celltypes or arrayData.celltypes for this tissue
#       (optional, default is to used the max values among all the cell tyoes for the tissue. use [] if you want to supply more arguments)
# noGeneScore: score for reactions without genes (optional, default -2)
# isozymeScoring: determines how scores are calculated for reactions with multiple genes joined by "OR" expression(s)
#       ('min','max','median','average') (optional, default 'max')
# complexScoring: determines how scores are calculated for reactions with multiple genes joined by "AND" expression(s)
#       ('min','max','median','average') (optional, default 'min')
# multipleCellScoring: determines how scores are calculated when several cell types are used
#       ('max' or 'average') (optional, default 'max')
# hpaLevelScores: structure wuth numerical scores for the expression level categories from HPA. The structure shohuld have a "names" and a "scores" field
#       (optional, see code for default scores)

## OUTPUT
# rxnScores: scores for each of the reactions in model
# geneScores: scores for each of the genes in model. Gene which are not in the dataset(s) have -Inf as scores
# hpaScores: scores for each of the genes in model if only taking hpaData into account. Genes which are not in the dataset(s) have -Inf as scores
# arrayScores scores for each of the genes in model if only taking arrayData into account. Genes which are not in the dataset(s) have -Inf as scores

import inspect
import math
import re
from collections import defaultdict

import numpy as np
from networkx.classes import is_empty
from numpy.ma.extras import unique
from plotly.utils import node_generator
from scipy.sparse import csr_matrix
from sympy.codegen.cfunctions import isnan


def score_complex_model(
    Model, HpaData, ArrayData, Tissue, Celltype, NoGeneScore, IsozymeScoring, ComplexScoring, MultipleCellScoring, HpaLevelScores
):
    model = Model
    hpaData = HpaData
    arrayData = ArrayData
    tissue = Tissue
    celltype = Celltype
    noGeneScore = NoGeneScore
    isozymeScoring = IsozymeScoring
    complexScoring = ComplexScoring
    multipleCellScoring = MultipleCellScoring
    hpaLevelScores = HpaLevelScores
    sig = inspect.signature(score_complex_model)
    num_args = len(sig.parameters)

    if num_args < 5:
        celltype = []
    if num_args < 6 or is_empty(noGeneScore):
        noGeneScore = -2
    if num_args < 7 or is_empty(isozymeScoring):
        isozymeScoring = "max"
    if num_args < 8 or is_empty(complexScoring):
        complexScoring = "min"
    if num_args < 9 or is_empty(multipleCellScoring):
        multipleCellScoring = "max"
    if num_args < 10:
        hpaLevelScores.names = ["High", "Medium", "Low", "None", "Strong", "Moderate", "Weak", "Negative", "Not detected"]
        hpaLevelScores.scores = [20, 15, 10, -8, 20, 15, 10, -8, -8]
    if is_empty(hpaData) and is_empty(arrayData):
        pass
    if isozymeScoring.lower() not in ["min", "max", "median", "average"]:
        raise ValueError("Valid options for isozymeScoring are 'min', 'max', 'median', and 'average'.")
    if complexScoring.lower() not in ["min", "max", "median", "average"]:
        raise ValueError("Valid options for complexScoring are 'min', 'max', 'median', and 'average'.")
    if multipleCellScoring.lower() not in ["max", "average"]:
        raise ValueError("Valid options for multipleCellScoring are 'max' and 'average'.")

    # Throw an error if array data for only one tissue is supplied without specifying threshold values
    if arrayData and len(unique(arrayData)) < 2 and ("threshold" not in arrayData or is_empty(arrayData.threshold)):
        raise ValueError(
            " arrayData must contain measurements for at least two celltypes/tissues since the "
            "score is calculated based on the expression level compared to the overall average"
        )

    # Process arrayData.threshold if necessary
    if "threshold" in arrayData and len(arrayData["threshold"]) == 1:
        # if only a single gene threshold value is provided, then just duplicate this value for all genes.
        arrayData["threshold"] = arrayData["threshold"] * len(arrayData["genes"])

    # This is so the code can ignore which combination of input data is used
    if not arrayData:
        arrayData = {"genes": [], "tissues": [], "celltypes": [], "levels": [], "threshold": []}

    # Initialize hpaData if empty
    if not hpaData:
        hpaData = {
            "genes": [],
            "tissues": [],
            "celltypes": [],
            "levels": [],
            "types": [],
            "reliabilities": [],
            "gene2Level": [],
            "gene2Type": [],
            "gene2Reliability": [],
        }

    # Check that the tissue exists (case-insensitive)
    tissue_upper = tissue.upper()

    in_hpa = any(t.upper() == tissue_upper for t in hpaData["tissues"])
    in_array = any(t.upper() == tissue_upper for t in arrayData["tissues"])

    if not (in_hpa or in_array):
        raise ValueError("The tissue name does not match")  # equivalent to dispEM

    if any(celltype):
        # Check that both data types has cell type defined if that is to be used
        if "celltypes" not in hpaData.keys() or "celltypes" not in arrayData.keys():
            raise ValueError("Both hpaData and arrayData must contain cell type information if cell type is to be used")
        if celltype.upper() not in hpaData["celltypes"].upper() and celltype.upper() not in arrayData["celltypes"].keys().upper():
            raise ValueError("The cell type ma,e does not match")

    # Some preprocessing of the structures to increase efficiency
    # Remove all tissues that are not the correct one
    J = [t.lower() != tissue.lower() for t in hpaData["tissues"]]

    # If cell type is supplied, then only keep that cell type
    if celltype:
        J_celltype = [ct.lower() != celltype.lower() for ct in hpaData.get("celltypes", [])]
        J = [j or j_ct for j, j_ct in zip(J, J_celltype)]

    J_array = np.array(J)

    hpaData["tissues"] = [t for t, keep in zip(hpaData["tissues"], ~J_array)]
    if "celltypes" in hpaData:
        hpaData["celltypes"] = [ct for ct, keep in zip(hpaData["celltypes"], ~J_array)]

    for key in ["gene2Level", "gene2Type", "gene2Reliability"]:
        if key in hpaData and isinstance(hpaData[key], np.ndarray):
            hpaData[key] = hpaData[key][:, ~J_array]

    # Remove all genes from the structures that are not in model or that aren't measured in the tissue
    if hpaData["genes"]:  # Check if it is not empty
        genes_array = np.array(hpaData["genes"])
        gene_mask = np.isin(genes_array, model["genes"], invert=True)  # genes NOT in model.genes

        level_mask = np.sum(hpaData["gene2Level"], axis=1) == 0  # genes with all zeros
        I = gene_mask | level_mask  # combine both masks
    else:
        I = np.array([], dtype=bool)

    # Remove from hpaData fields
    if I.size > 0:
        hpaData["genes"] = [g for g, keep in zip(hpaData["genes"], ~I)]

        if "gene2Level" in hpaData:
            hpaData["gene2Level"] = hpaData["gene2Level"][~I, :]
        if "gene2Type" in hpaData:
            hpaData["gene2Type"] = hpaData["gene2Type"][~I, :]
        if "gene2Reliability" in hpaData:
            hpaData["gene2Reliability"] = hpaData["gene2Reliability"][~I, :]

    # Fine matching tissue (and optionally celltype)
    tissue_mask = [t.lower() == tissue.lower() for t in arrayData["tissues"]]

    if celltype:  # If specific celltype is supplied
        tissue_mask = [tm and (ct.lower() == celltype.lower()) for tm, ct in zip(tissue_mask, arrayData["celltypes"])]

    tissue_mask_np = np.array(tissue_mask)
    I_coles = tissue_mask_np

    # Create row mask: genes not in model or all NaN in matching tissue(s)
    genes_array = np.array(arrayData["genes"])
    not_in_model = np.isin(genes_array, model["genes"], invert=True)
    no_expression = np.all(np.isnan(arrayData["levels"][:, I_coles]), axis=1)

    J = not_in_model | no_expression

    # Filter arrayData
    arrayData["genes"] = [g for g, keep in zip(arrayData["genes"], ~J)]
    arrayData["levels"] = arrayData["levels"][:, ~J]

    if "threshold" in arrayData:
        arrayData["threshold"] = [th for th, keep in zip(arrayData["threshold"], ~J)]

    # Calculate the scores for the arrayData. These scores are calculated for each genes from its fold change between
    # the tissue/celltype(s) in question and all other celltypes, or the threshold if supplied. This is a lower quality data
    # than protein abundance, since gene abundance is an indirect estimate of protein level. These scores are therefore
    # only used for genes for which there is no HPA data available. The fold changes are transformed as min(5*log(x), 10) ]
    # for x > 1 and max(5*log(x),-5) for x < 1 in order to have negative scores for lower expressed genes and to scale
    # the scores to have somewhat lower weights than the HPA scores

    ###----------------Helper functions----------------###

    def scoreSimpleRule(rule, genes, gScores, isozymeScoring, complexScoring):
        ruleGenes = re.findall(r"[^&|() ]+", rule)
        geneInd = [genes.index(g) for g in ruleGenes if g in genes]
        if not geneInd:
            return np.nan
        scores = gScores[geneInd]
        method = isozymeScoring if "|" in rule else complexScoring

        if method == "min":
            return np.nanmin(scores)
        elif method == "max":
            return np.nanmax(scores)
        elif method == "median":
            return np.nanmedian(scores)
        elif method == "average":
            return np.nanmean(scores)
        else:
            raise ValueError("Invalid scoring method.")

    def scoreComplexRule(rule, genes, gScores, isozymeScoring, comlexScoring):
        search_phrases = [r"\([^&|() ]+( & [^&|() ]+)+\)", r"\([^&|() ]+( \| [^&|() ]+)+\)"]
        subsets = []
        c = 1
        r_orig = rule
        for _ in range(100):
            for phrase in search_phrases:
                new_subset = re.findall(phrase, rule)
                if new_subset:
                    subsets.extend(new_subset)
                    for i, s in enumerate(new_subset):
                        tag = f"#{c + i}#"
                        rule = rule.replace(s, tag, 1)
                    c += len(new_subset)
            if rule == r_orig:
                break
            r_orig = rule
        subsets.append(rule)

        for i, sub in enumerate(subsets):
            score = scoreSimpleRule(sub, genes, gScores, isozymeScoring, complexScoring)
            gScores = np.append(gScores, score)
            genes = genes + [f"#{i  + 1}#"]
        return gScores[-1]

    def find_nonzero(martrix):
        coo = csr_matrix(martrix).tocoo()
        return coo.row, coo.col, np.arange(len(coo.data))

    def unique_with_index(lst):
        seen = {}
        out = []
        index = []
        for i, val in enumerate(lst):
            if val not in seen:
                seen[val] = len(out)
                out.append(val)
            index.append(seen[val])
        return out, index

    ###----------------Start of scoring logic----------------###

    # Handle array data
    tempArrayLevels = np.copy(arrayData["levels"])
    tempArrayLevels[np.isnan(tempArrayLevels)] = 0

    if "threshold" in arrayData and len(arrayData["threshold"]) > 0:
        average = np.array(arrayData["threshold"])
    else:
        average = np.sum(tempArrayLevels, axis=1) / np.sum(~np.isnan(tempArrayLevels), axis=1)

    I = arrayData["tissue_mask"]  # Precomputed mask for tissue/celltype of interest

    if multipleCellScoring.lower() == "max":
        current = np.max(tempArrayLevels[:, I], axis=1)
    else:
        current = np.sum(tempArrayLevels[:, I], axis=1) / np.sum(~np.isnan(arrayData["levels"][:, I]), axis=1)

    if current.size > 0:
        with np.errstate(divide="ignore", invalid="ignore"):
            aScores = 5 * np.log(current / average)
            aScores[np.isnan(aScores)] = -5
        aScores[aScores > 0] = np.minimum(aScores[aScores > 0], 10)
        aScores[aScores < 0] = np.maximum(aScores[aScores < 0], -5)
    else:
        aScores = np.array([])

    # Map HPA levels to scores
    level_names_upper = [s.upper() for s in hpaLevelScores["names"]]
    hpa_levels_upper = [s.upper() for s in hpaData["levels"]]
    J = [level_names_upper.index(x) if x in level_names_upper else -1 for x in hpa_levels_upper]

    if any(j == -1 for j in J):
        raise ValueError("There are expression level categories that do not match to hpaLevelScores.")

    scores = np.array([hpaLevelScores["scores"][j] for j in J])
    K, L, M = find_nonzero(hpaData["gene2Level"])

    shape = (len(hpaData["genes"]), len(hpaData["tissues"]))
    sparse_scores = csr_matrix((scores[M], (K,L)), shape=shape)

    if multipleCellScoring.lower() == "max":
        hScores = np.max(sparse_scores.toarray(), axis=1)
    else:
        hScores = np.mean(sparse_scores.toarray(), axis=1)

    # Assign scores to model.gees
    geneScores = np.full(len(model["genes"]), np.nan)
    hpaScores = np.full(len(model["genes"]), -np.inf)
    arrayScores = np.full(len(model["genes"]), -np.inf)

    for i, g in enumerate(model["genes"]):
        if g in arrayData["genes"]:
            j = arrayData["genes"].index(g)
            geneScores[i] = aScores[j]
            arrayScores[i] = aScores[j]
        if g in hpaData["genes"]:
            j = hpaData["genes"].index(g)
            geneScores[i] = hScores[j]
            hpaScores[i] = hScores[j]

    # Score unique grRules
    uRules, rule_ind = unique_with_index(model["grRules"])
    uScores = np.full(len(uRules), np.nan)

    uRules = [r.replace(" and ", " & ").replace(" or ", " | ") for r in uRules]

    for i, rule in enumerate(uRules):
        if not rule:
            uScores[i] = noGeneScore
        elif "&" in rule and "|" in rule:
            uScores[i] = scoreComplexRule(rule, model["genes"], geneScores, isozymeScoring, complexScoring)
        else:
            uScores[i] = scoreSimpleRule(rule, model["genes"], geneScores, isozymeScoring, complexScoring)

    uScores[np.isnan(uScores)] = noGeneScore
    rxnScores = [uScores[i] for i in rule_ind]

    ###----------------End of scoring logic----------------###