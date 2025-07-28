# This is the Python version of removeLowScoreGenes (used in ftinit) originally written in MATLAB.
# Source code: https://github.com/SysBioChalmers/RAVEN/blob/cd82b5c6ba26a88bed2a9ca179c682a6434fd541/INIT/removeLowScoreGenes.m

# remove_low_score_genes: remove low-scoring genes from model
# This function removes genes from a model based on their scores, a step used by the tINIT package. The function recognizes and differentiates
# between isozymes and subunits of an enzyme complex. Genes are removed from each grRule, subject to the following conditions:
# 1) At least one gene must remain associated with the reaction
# 2) Genes involved in a complex (joined by ANDs) are not removed

## INPUT
# model: Model structure from which genes are to be removed
# gene_score: A vector of scores associated with the model genes. Genes with a positive score will remain in the model, whereas
#             genes with a negative score will try to be removed.
#             If all genes associated with a reaction have a negative score, then the least-negative gene will remain; if there is a tie,
#             one will be selected at random.
#             If a negative-scoring gene is a subunit in a complex, it will not be removed; however, the entire complex may be removed.
#             See the following example cases:
#              Original: G1 or (G2 and G3 and G4)
#              Negative: G1, G2
#              New: G2 and G3 and G4
#
#              Original