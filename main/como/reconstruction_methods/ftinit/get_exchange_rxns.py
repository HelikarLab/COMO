# This is the Python version of getExchangeRxns (used in ftinit) originally written in MATLAB.
# Source code: https://github.com/SysBioChalmers/RAVEN/blob/cd82b5c6ba26a88bed2a9ca179c682a6434fd541/core/getExchangeRxns.m

# get_exchange_rxns: Retreives the exchange reactions from a model. Exchange reactions are identified by having either no substrates or products.

## INPUT
# model: a model structure
# reaction_type: which exchange reactions should be returned
                # 'all': all reactions, irrespective of reaction bounds
                # 'uptake': reactions with bounds that imply that only uptake are allowed. Reaction direction, upper and lower bounds are all considered.
                # 'excrete': reactions with bounds that imply that only excretion are allowed. Reaction direction, upper and lower bounds are all considered.
                # 'reverse': reactions with non-zero upper and lower bounds that imply that both uptake and excretion are allowed.
                # 'blocked': reactions that have zero upper and lower bounds, not allowing any flux
                # 'in': reactions where the boundary metabolite is the substrate of the reavtion, a positie flux value would  imply uptake, but reaction bounds are not considered.
                # 'out': reactions where the boundary metabolite is the substrate of the reaction, a positive flux value would imply uptake but reaction bounds are not considered.

## OUTPUT
# exchange_rxns: cell array with the IDs of the exchange reactions
# exchane_rxns_indexes: vector with the indexes of the exchange reactions

# Note: The union of 'in' and 'out' equals 'all'. Also, the union of 'uptake', 'excrete', 'reverse' and 'blocked' equals all.

from typing import Tuple, List, Union
import numpy as np
from cobra import Model, Reaction
from sympy.stats import Expectation


def get_exchange_rxns(model: Model, reaction_type: str = "all") -> Tuple[List[str], List[int]]:
    reaction_type = str(reaction_type).lower()

    # Precompute stochiometric sign flags per reactions
    # For each reaction, inspect metabolite coefficients:
    n_rxns = len(model.reactions)
    has_no_prod = np.zeros(n_rxns, dtype=bool) # no positive coefficient (>0)
    has_no_subs = np.zeros(n_rxns, dtype=bool) # no negative coefficient (<0)

    # If the model has a custom 'unconstrained' attribute (MATLAB) variant, try to mimic:
    # In COBRApy, this is uncommon; we default to full S-check if not present
    if hasattr(model, "unconstrained") and model.unconstrained is not None:
        # Assume model.unconstrained is a boolean vector on metabolites (MATLAB-style)
        # Map metab unconstrained -> use only those rows. If shape mismatch, fallback.
        try:
            unconstrained_mask = np.asarray(model.unconstrained).astype(bool)
            # Build per-reaction lists from metabolite list order
            mets = list(model.metabolites)
            if len(unconstrained_mask) == len(mets):
                for j, rxn in enumerate(model.reactions):
                    # get stoich for unconstrained metanolites only
                    coeffs = []
                    for i, met in enumerate(mets):
                        if not unconstrained_mask[i]:
                            continue
                        coeff = rxn.metabolites.get(met, 0.0)
                        coeffs.append(coeff)
                    coeffs = np.array(coeffs) if coeffs else np.array([0.0])
                    has_no_prod[j] = np.all(coeffs <= 0)
                    has_no_subs[j] = np.all(coeffs >= 0)
            else:
                # fallback to full-model check
                raise ValueError("unconstrained length mismatch, falling back to full S check")
        except Exception:
            # fallback
            for j, rxn in enumerate(model.reactions):
                coeffs = np.array(list(rxn.metabolites.values())) if rxn.metabolites else np.array([0.0])
                has_no_prod[j] = np.all(coeffs <= 0)
                has_no_subs[j] = np.all(coeffs >= 0)
    else:
        # Standard path: inspect all stochiometric coefficients for each reaction
        for j, rxn in enumerate(model.reactions):
            coeffs = np.array(list(rxn.metabolites.values())) if rxn.metabolites else np.array([0.0])
            # Note: 'np product' means no coeffs > 0
            has_no_prod[j] = np.all(coeffs <= 0)
            # 'no substrate' means no coeffs < 0
            has_no_subs[j] = np.all(coeffs >= 0)

    # union (MATLAB used vertical concatenation of indexes)
    all_exch_mask = np.logical_or(has_no_prod, has_no_subs)
    all_exch_idx = np.nonzero(all_exch_mask[0].tolist)
    has_no_prod_idx = np.nonzero(has_no_prod)[0].tolist()
    has_no_subs_idx = np.nonzero(has_no_subs)[0].tolist()

    # Choose according to reaction_type (MATLAB switch)
    if reaction_type in ("both", "all"):
        sel_idx = all_exch_idx
    elif reaction_type == "in":
        sel_idx = has_no_subs_idx
    elif reaction_type == "out":
        sel_idx = has_no_prod_idx
    elif reaction_type == "blocked":
        # lb == 0 and ub == 0
        sel_idx = [i for i in all_exch_idx if float(model.reactions[i].lower_bound) == 0 and float(model.reactions[i].upper_bound) == 0]
    elif reaction_type == "uptake":
        # combine conditions from MATLAB
        # (model.lb(hasNoProd) < 0 & model.ub(hasNoProd) <= 0) OR (model.lb(hasNoSubs) >= 0 & model.ub(hasNoSubs) > 0)
        idx_a = [i for i in has_no_prod_idx if float(model.reactions[i].lower_bound) < 0 and float(model.reactions[i].upper_bound) <= 0]
        idx_b = [i for i in has_no_subs_idx if float(model.reactions[i].lower_bound) >= 0 and float(model.reactions[i].upper_bound) >0]
        sel_idx = list(dict.fromkeys(idx_a + idx_b)) # preserve uniqueness
    elif reaction_type == "excrete":
        # (model.lb(hasNoProd) >= 0 & model.ub(hasNoProd) > 0) OR (model.lb(hasNoSubs) < 0 & model.ub(hasNoSubs) <= 0)
        idx_a = [i for i in has_no_prod_idx if float(model.reactions[i].lower_bound) >= 0 and float(model.reactions[i].upper_bound) > 0]
        idx_b = [i for i in has_no_subs_idx if float(model.reactions[i].lower_bound) < 0 and float(model.reactions[i].upper_bound) <= 0]
        sel_idx = list(dict.fromkeys(idx_a + idx_b))
    else:
        raise ValueError(f"Invalid reaction_type specified: {reaction_type}")
    # Sort indices to match MATLAB sort(exchange_rxn_indexes)
    sel_idx_sorted = sorted(sel_idx)

    # Retuen reaction IDs and 0-based indices
    sel_rxn_ids = [model.reactions[i].id for i in sel_idx_sorted]
    return sel_rxn_ids, sel_idx_sorted

