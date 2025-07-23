# This is the Python version of ftINITInternalAlg.m (used in ftinit) originally written in MATLAB.
# Source code: https://github.com/SysBioChalmers/RAVEN/blob/cd82b5c6ba26a88bed2a9ca179c682a6434fd541/INIT/ftINITInternalAlg.m

# ftINITInternalAlg: This function runs the MILP for a step in ftINIT

## INPUT
# model: a reference model structure
# rxn_scores: a vector of scores for the reactions in the model.
#             positive scores are reactions to keep and negative scores are reactions to exclude. rxns set to 0 are excluded from the problem.
# met_data: boolean matrix with mets as rows and rxns as columns saying which reaction produces each detected met (optional, default [])
# essential_rxns: cell array of reactions that are essential and that have to be in the resulting model. this is normally
#                 used when fitting a model to task (fee fitTasks)
# prod_weight: a score that determines the value of having net-production of metabolites. this is a way of having a more functional
#              network as it provides a reason on for including bad reactions for connectivity reasons. this score is for each metabolite,
#              and the sum of these weights and the scores for the reactions is what is optimized
# allow_excretion: true if excretion of all metabolites should be allowed. this results in fewer reactions being considered dead-ends,
#                  but all reactions in the resulting model may not be able to carry flux. if this is "false" then the equality constraints
#                  are taken from model.b. if the input model lacks exchange reactions then this should probably be "true", or a large proportion
#                  of the model would be excluded for connectivity reasons.
# rem_posrev: if true, the positive reversible reactions are removed from the problem. this is used in step 1 of ftINIT.
# params: parameters for the MILP, for example MIPGap and TimeLimit
# start_vals: Start values for the MILP, typically used when rerunning with a higher MIPGap, to use the results from the previous run
# fluxes: fluxes from the last run
# verbose: if true, the MILP progression will be shown.

## OUTPUT
# deleted_rxns: reactions which were deleted by the algorithm (only rxns included in the problem)
# met_production: array that indicates which of the metabolites in presentMets that could be produced
#                 0: metabolite could not be produced
#                 1: metabolite could be produced
# res: the result from the MILP
# turned_on_rxns: the reactions determined to be present (only rrxns included in the problem)
# fluxes: the fluxes from the MILP

# This function is the actual implementation of the algorithm. See ftINIT for a higher-level function for model reconstruction.

import numpy as np

def ftinit_internal_alg(model, rxn_scores, met_data, essential_rxns, prod_weight, allow_excretion, rem_pos_rev, params, start_vals, fluxes, verbose):
    essential_rxns = params.get("essential_rxns", [])
    essential_rxns = np.array(essential_rxns).flatten()
    prod_weight = params.get("prod_weight", 0.5)

    # The model should be in the reversible format and all relevant exchange reactions should be open
    if 'constrained' in model:
        print("Exchange metabolites are still present in the model. Use simplifyModel if this is not intended")
    essential = np.isin(model.reactions, essential_rxns)

    # Some nice to have numbers
    n_mets = len(model.metabolites)
    n_rxns = len(model.reactions)
    nrxns_with_on_off = n_rxns

    # Reactions with score 0 will just be left in the model, and will not be a part of the problem (but can carry flux).
    # It is possible to set score = 0 for e.g. spontaneous reactions, exchange rxns, etc., which may not be that interesting to remove

    # If makeIrrev is on, we can just as well skip all positive reversible rxns - they will be turned on without carrying
    # any flux since they can form a loop within themself (fwd-rev)
    rem_pos_rev = params.get("rem_pos_rev", False)
    if rem_pos_rev:
        rxn_scores[(rxn_scores > 0) & (model['rev'] != 0)] = 0

    # Handle metabolomics
    # A problem with the metabolomics is that some of the producer reactions for a metabolite could be excluded from the problem.
    # We solve this by adding them with score 0. They can still be seen as positive, since they either don't matter (there is another producer)
    # or they can contribute to an increased score.
    rev_rxns = model.rev != 0
    ess_rev_rxns = np.where(rev_rxns & essential)[0]
    ess_irrev_rxns = np.where(~rev_rxns & essential)[0]

    met_data = params.get("met_data", None)
    if met_data is not None:
        # Remove any metData rows that are connected to an essential rxn
        contains_essential = np.any(met_data[:,ess_rev_rxns], axis=1)
        met_data = met_data[~contains_essential, :]

    if met_data is not None:
        met_rxns = np.any(met_data, axis = 0)
        pos_rxns = (rxn_scores > 0) | ((rxn_scores == 0) & met_rxns)
    else:
        pos_rxns = rxn_scores > 0

    return deleted_rxns, met_production, res, turned_on_rxns, fluxes
