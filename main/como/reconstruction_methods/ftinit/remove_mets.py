# This is the Python version of removeMets (used in ftinit) originally written in MATLAB.
# Source code: https://github.com/SysBioChalmers/RAVEN/blob/cd82b5c6ba26a88bed2a9ca179c682a6434fd541/core/removeMets.m

# remove_mets: deletes a set of metabolites from a model

## INPUT
# model: a model structure
# mets_to_remove: either a cell array of metabolite IDs, a logical vector with the same number of elements as metabolites
#                 in the model, of the vector of indexes to remove
# is_names: true if the supplied mets represent metabolite names (as opposed to IDs). This is a way delete metabolites in
#                several compartments at once without knowing the exact IDs. This only works if mets_to_remove is a cell array (optional, default false)
# remove_unused_rxns: remove reactions that are no longer in use (optional, default false)
# remove_unused_genes: remove genes that are no longer in use (optional, default false)
# remove_unused_comps: remove compartments that are no longer in use (optional, default false)

## OUTPUT
# reduced_model: an updated model structure

from typing import Union, List
import numpy as np
from cobra import Model, Metabolite

from reconstruction_methods.ftinit.remove_reactions import remove_reactions


def remove_mets(model: Model, mets_to_remove: Union[List[str], List[Metabolite], List[bool], List[int]], is_names: bool = False, remove_unused_rxns: bool = False, remove_unused_genes: bool = False, remove_unused_comps: bool = False) -> Model:
    reduced_model = model.copy()

    # Convert mets_to_remove to a list of metabolite IDs
    if isinstance(mets_to_remove, list) and len(mets_to_remove) > 0:
        if isinstance(mets_to_remove[0], Metabolite):
            met_ids = [m.id for m in mets_to_remove]
        elif isinstance(mets_to_remove[0], bool) or isinstance(mets_to_remove[0], np.bool_):
            if len(mets_to_remove) != len(reduced_model.metabolites):
                raise ValueError("Boolean mask length must match number of metabolites")
            met_ids = [met.id for met, flag in zip(reduced_model.metabolites, mets_to_remove) if flag]
        elif isinstance(mets_to_remove[0], int):
            met_ids = [reduced_model.metabolites[i].id for i in mets_to_remove]
        else: # assume strings
            if is_names:
                met_ids = [m.id for m in reduced_model.metabolites if m.name in mets_to_remove]
            else:
                met_ids = list(mets_to_remove)
    else:
        met_ids = []

    # Remove metabolites
    if met_ids:
        reduced_model.remove_metabolites(met_ids, destructive=True)

    # Remove unused reactions
    if remove_unused_rxns:
        unused_rxns = [rxn.id for rxn in reduced_model.reactions if len(rxn.metabolites) == 0]
        if unused_rxns:
            reduced_model = remove_reactions(reduced_model, remove_unused_rxns, remove_unused_mets=False, remove_unused_genes=remove_unused_genes, remove_unused_comps=remove_unused_comps)

    # Remove unused compartments
    if remove_unused_comps:
        all_used_comps = {met.compartment for met in reduced_model.metabolites}
        unused_comps = [comp for comp in reduced_model.compartments if comp not in all_used_comps]
        for comp in unused_comps:
            reduced_model.compartments.pop(comp,None)

    return reduced_model

