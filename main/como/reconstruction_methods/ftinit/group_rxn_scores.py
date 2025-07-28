# This is the Python version of groupRxnScores (used in ftinit) originally written in MATLAB.
# Source code: https://github.com/SysBioChalmers/RAVEN/blob/cd82b5c6ba26a88bed2a9ca179c682a6434fd541/INIT/groupRxnScores.m

# group_rxn_scores: this function sums up the reaction scores for all reactions that were merged into one by the linear merge

## INPUT
# model: the model with linearly merged rxns
# orig_rxn_scores: the rxn_scores from the model before the linear merge
# orig_rxn_ids: the rxn ids of the model before the linear merge
# group_ids: the group_ids vector output from linear_merge. there is one integer for each rxn in orig_rxn_ids. 0 means the reaction was
# merged. a non-zero integer means that the reactions was merged with all other rxns having the same integer.
# orig_rxns_to_zero: a logical vector saying which of the original rxns that should not be part of the problem. the way this is solved
# is that all such reactions have a rxn_score of 0. if any original rxn_score value should be zero (which is very unlikely) it is changed to 0.001.
# If the sum of the rxn_scores for a merged rxn becomes zero while some of them are non-zero, the new value will also be 0.01, to distinguish the rxn
# from rxns with only rxns to zero. there are two reasons why we don't want zeros in the reaction scores unless these reactions should be ignored:
#   1) we want to be able to separate those
#   2) it is difficult to handle a zero value in the MILP - the on/off of such a reaction can be random, so better to fix it in one direction

import numpy as np


def group_rxn_scores(model, orig_rxn_scores, orig_rxn_ids, group_ids, orig_rxns_to_zero):
    new_rxn_scores = np.zeros(len(model.reactions))

    # find indices of original reactions in the model
    ia, ib = np.intersect1d(model.reactions, orig_rxn_ids, return_indices=True)

    group_ids_merged = np.full(len(model.reactions), np.nan)
    group_ids_merged[ia] = group_ids[ib]

    # Check if any of the original scores are 0, change them to 0.01
    orig_rxn_scores[orig_rxn_scores == 0] = 0.01
    # Set the reaction scores for reactions to zero to 0
    orig_rxn_scores[orig_rxns_to_zero] = 0

    # Fill in original scores
    new_rxn_scores[ia] = orig_rxn_scores[ib]

    for i in range(len(model.reactions)):
        # For reactions that are not merged with anything, just keep scores as it is
        if group_ids_merged[i] != 0:
            # Find all original reactions in the group
            sel = group_ids == group_ids_merged[i]
            new_rxn_scores[i] = np.sum(orig_rxn_scores[sel])

            if new_rxn_scores[i] == 0 and np.any(orig_rxn_scores[sel] != 0):
                # Special care where the reactions sum to 0 while some are non-zero
                new_rxn_scores[i] = 0.01
    return new_rxn_scores
