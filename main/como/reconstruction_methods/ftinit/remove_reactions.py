# This is the Python version of removeReactions (used in ftinit) originally written in MATLAB.
# Source code: https://github.com/SysBioChalmers/RAVEN/blob/cd82b5c6ba26a88bed2a9ca179c682a6434fd541/core/removeReactions.m

# remove_reactions: deletes a set of reacgtions from a model

## INPUT
# model: a model structure
# rxns_to_remove: either a cell array of reaction IDs, a logical vector with the same number of elements as reactions in the model, or a vector of indexes to remove
# remove_unused_mets: remove metabolites that are no longer in use (optional, default false)
# remove_unused_genes: remove genes that are no longer in use (optional, default false)
# remove_unused_comps: remove compartments that are no longer in use (optional, default false)

## OUTPUT
# reduced_model: an updated model structure

