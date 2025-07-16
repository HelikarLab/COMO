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
import numpy as np
from jupyter_server.auth import passwd

from networkx.classes import is_empty
from numpy.ma.extras import unique
from plotly.utils import node_generator
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
        arrayData = {
            "genes": [],
            "tissues": [],
            "celltypes": [],
            "levels": [],
            "threshold": []
        }

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
            "gene2Reliability": []
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
        if celltype.upper() not in hpaData["celltypes"].keys().upper() and celltype.upper() not in arrayData["celltypes"].keys().upper():
            raise ValueError("The cell type ma,e does not match")

    # Some preprocessing of the structures to increase efficiency
    # Remove all tissues that are not the correct one
    J = [t.lower() != tissue.lower() for t in hpaData["tissues"]]

    # If cell type is supplied, then only keep that cell type
    if celltype:
        J_celltype = [ct.lower() != celltype.lower() for ct in hpaData.get("celltypes",[])]
        J = [j or j_ct for j, j_ct in zip(J, J_celltype)]

    J_array = np.array(J)

    hpaData["tissues"] = [t for t, keep in zip(hpaData["tissues"], ~J_array)]
    if "celltypes" in hpaData:
        hpaData["celltypes"] = [ct for ct, keep in zip(hpaData["celltypes"], ~J_array)]

    for key in ["gene2Level", "gene2Type", "gene2Reliability"]:
        if key in hpaData and isinstance(hpaData[key], np.ndarray):
            hpaData[key] = hpaData[key][:,~J_array]

    # Remove all genes from the structures that are not in model or that aren't measured in the tissue
    if hpaData["genes"]: # This should not be necessary, but the summation is a 0x1 matrix and the other is []
        pass



    # If cell type is supplied, then only keep that cell type
    if any(celltype):
        pass

    # Remove all genes from the structures that are not in model of that aren't measured in the tissue
    pass


    # Calculate the scores for the arrayData. These scores are calculated for each genes from its fold change between
    # the tissue/celltype(s) in question and all other celltypes, or the threshold if supplied. This is a lower quality data
    # than protein abundance, since gene abundance is an indirect estimate of protein level. These scores are therefore
    # only used for genes for which there is no HPA data available. The fold changes are transformed as min(5*log(x), 10) ]
    # for x > 1 and max(5*log(x),-5) for x < 1 in order to have negative scores for lower expressed genes and to scale
    # the scores to have somewhat lower weights than the HPA scores

    tempArrayLevels = arrayData["levels"]
    #tempArrayLevels(isnan(tempArrayLevels)) = 0

    if "threshold" in arrayData and arrayData["threshold"]:
        pass
    else:
        pass
    """
    more code comes here
    """
    aScores(aScores > 0) = min(aScores(aScores > 0),10)
    aScores(aScores < 0) = max(aScores(aScores < 0), -5)
    aScores(isnan(aScores)) = -5 # NaNs occur when gene expression is zero across all tissues

    # Map the HPA levels to scores
    """
    more code comes here
    """

    # Assign gene scores, prioitizing HPA (protein) data over arrayData (RNA)

    # To speed things up, only need to score each unique grRule once

    # convert logic operators to symbols

    # score based on prescence/combination of & and | operators

    # NaN reaction scores should be changed to the no-gene score

    # re-map unique rule scores to model


    # Score reactions with simple gene rules (those with all ANDs or all ORs)


    # Score reactions with complex gene rules (those with both ANDs and ORs)

    # Specify phrases to search for in hte grRule. These phrases will fine genes grouped by all ANDs (first phrase) or all ORs (second phrase)

    # initialize some variables


    # score each subset and append to gene list and gene scores

    # the final subset score in the overall reaction score