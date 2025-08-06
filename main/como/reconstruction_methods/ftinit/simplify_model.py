# This is the Python version of ftINITFillGaps (used in ftinit_fill_gaps and prep_init_model) originally written in MATLAB.
# Source code: https://github.com/SysBioChalmers/RAVEN/blob/cd82b5c6ba26a88bed2a9ca179c682a6434fd541/core/simplifyModel.m

# Simplify model: Simplifies a model by deleting reactions/metabolites

## INPUT
# model: a model structure
# delete_unconstrained: delete metabolites marked as unconstrained (optional, default true)
# delete_duplicates: deleted all but one of duplicate reactions (optional, default false)
# delete_zero_interval: delete reactions that are constrained to zero flux (optional, default false)
# delete_inaccessible: delete dead end reactions (optional, default false)
# delete_min_max: deleted reactions that cannot carry flux by trying to minimize/maximize the flux through that reaction.
#                 May be time consuming (optional, default false)
# group_linear: group linearly dependent pathways (optional, default false)
# constrain_reversible: check if there are reversible reactions which can only carry flux in one direction, and if so constrain them
#                       to be irreversible. This tends to allow for more reaction grouped when using group_linear (optional, default false)
# reverse_rxns: cell array with reaction IDs that are not allowed to be removed (optional)
# supress_warnings: true if warnings should be suppressed (optional, default false)

## OUTPUT
# reduce_model: an updated model structure
# deleted_reactions: a cell array with the IDs of all deleted reactions
# deleted_metabolites: a cell array with the IDs of all deleted metabolites

# This function is for reducing the model size by removing reactions and associated metabolites that cannot carry flux.
# It can also be used for identifying different types of gaps.

import cobra
from cobra.flux_analysis import flux_variability_analysis
from cobra.util.solver import linear_reaction_coefficients
import warnings

def simplify_model(model, delete_unconstrained=True, delete_duplicates=True, delete_zero_interval=True, delete_inaccessible=False,
                   delete_min_max=False, group_linear=False, constrain_reversible=False, reserved_rxns=None, supress_warnings=False):
    if reserved_rxns is None:
        reserved_rxns = set()
    else:
        reserved_rxns = set(reserved_rxns)

    reduced_model = model.copy()
    deleted_reactions = set()
    deleted_metabolites = set()

    if delete_unconstrained and hasattr(reduced_model, "unconstrained"):
        unconstrained_flags = getattr(reduced_model,"unconstrained")
        unconstrained_mets = [m for m, f in zip(reduced_model.metabolutes, unconstrained_flags) if f != 0]
        deleted_metabolites.update([met.id for met in unconstrained_mets])
        reduced_model.remove_metabolites(unconstrained_mets, destructive=True)
        delattr(reduced_model, "unconstrained")
