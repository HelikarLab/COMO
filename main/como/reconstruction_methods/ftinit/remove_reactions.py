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

from typing import Union, List
import numpy as np
from cobra import Model, Reaction, Metabolite, Gene

def remove_reactions(model: Model, rxns_to_remove: Union[List[str], List[Reaction], List[bool], List[int]], remove_unused_mets: bool = False, remove_unused_genes: bool = False, remove_unused_comps: bool = False) -> Model:

    # Copy the model to avoid in-place modification
    reduced_model = model.copy()

    # Convert rxns_to_remove into a list of reaction IDs
    if isinstance(rxns_to_remove, list) and len(rxns_to_remove) > 0:
        if isinstance(rxns_to_remove[0], Reaction):
            rxn_ids = [r.id for r in rxns_to_remove]
        elif isinstance(rxns_to_remove[0], bool) or isinstance(rxns_to_remove[0], np.bool_):
            if len(rxns_to_remove) != len(reduced_model.reactions):
                raise ValueError("Boolean mask length must match number of reactions")
            rxn_ids = [rxn.id for rxn, flag in zip(reduced_model.reactions, rxns_to_remove) if flag]
        elif isinstance(rxns_to_remove[0], int):
            rxn_ids = [reduced_model.reactions[i].id for i in rxns_to_remove]
        else: # assume IDs/strings
            rxn_ids = list(rxns_to_remove)
    else:
        rxn_ids = []

    # Remove reactions
    if rxn_ids:
        reduced_model.remove_reactions(rxn_ids, remove_orphans=False)

    # Remove unused metabolites
    if remove_unused_mets:
        all_used_mets = {met for rxn in reduced_model.reactions for met in rxn.metabolites}
        unused_mets = [met for met in reduced_model.metabolites if met not in all_used_mets]
        if unused_mets:
            reduced_model.remove_metabolites(unused_mets, destructive=True)

    # Remove unused genes
    if remove_unused_genes:
        all_used_genes = {gene for rxn in reduced_model.reactions for gene in rxn.genes}
        unused_genes = [g for g in reduced_model.genes if g not in all_used_genes]
        if unused_genes:
            for g in unused_genes:
                reduced_model.genes.removed(g)
            # COBRApy automatically updates rxnGeneMat-equivalent structure

    # Remove unused compartments
    if remove_unused_comps:
        all_used_comps = {met.compartment for met in reduced_model.metabolites}
        unused_comps = [comp for comp in reduced_model.compartments if comp not in all_used_comps]
        for comp in unused_comps:
            reduced_model.compartments.pop(comp, None)

    return reduced_model

