## This is a fast leak test originally written in MATLAB.
# source code: https://github.com/opencobra/cobratoolbox/blob/950f40cff4256d738fc8421a6ecfe91e084f46dd/src/reconstruction/modelGeneration/fastLeakTest.m

## INPUT
# model: model structure
# testRxns: List of exchange reactions to be tested for leaks
# demandTest (optional): if 'TRUE' is entered, demand reactions for all metabolites in the model are created
## OUTPUT
# LeakMets: List of exchange reactions for leaking metabolites
# modelClosed: Model structure that has been tested for leaks
# FLuxExV: Flux vector for computed exchange reactions in the closed model

import cobra
import numpy as np
from cobra.util.array import create_stoichiometric_matrix
from optlang.symbolics import Zero
from scipy.sparse import csc_matrix


def fast_leak_test_function(model, testRxns, demandTest):
    tol = 1e-6
    modelClosed = model.copy()

    # Get stoichiometric matrix in CSC format
    # Step 1: Get a dense NumPy matrix
    S_dense = create_stoichiometric_matrix(modelClosed)
    # Step 2: Convert it to a sparse CSC matrix
    S = csc_matrix(S_dense)

    rxns = modelClosed.reactions
    num_rxns = len(rxns)
    # Count non-zero stoichiometric entries per reaction (i.e., columns)
    num_nonzero_per_rxn = np.array((S != 0).sum(axis=0)).flatten()

    # Generate boolean masks (same as MATLAB logic)
    exp = np.array([(num_nonzero_per_rxn[i] == 1 and rxns[i].lower_bound >= 0) for i in range(num_rxns)])
    upt = np.array([(num_nonzero_per_rxn[i] == 1 and rxns[i].upper_bound <= 0) for i in range(num_rxns)])
    reversibleEX = np.array([(num_nonzero_per_rxn[i] == 1 and rxns[i].lower_bound < 0 and rxns[i].upper_bound > 0) for i in range(num_rxns)])

    count = exp | upt | reversibleEX
    # Get the reactions where count[i] is True
    ExR = [rxn for i, rxn in enumerate(modelClosed.reactions) if count[i]]
    # Set lower_bound = 0 for those reactions
    for i in range(num_rxns):
        if count[i]:
            rxns[i].lower_bound = 0
    modelexchangeAbbr = list(set(testRxns).union([rxn.id for rxn in ExR]))

    FluxEx = []
    cnt = 1
    rxnNames = []

    # test for all demand reactions is an option
    if demandTest:  # add demand reactions for all metabolites in model to check for those too
        rxnNames = []
        for metabolite in modelClosed.metabolites:
            demand_rxn = add_demand_reaction(modelClosed, metabolite)
            if demand_rxn is not None:
                rxnNames.append(demand_rxn.id)  # collect demand reaction IDs
        # return modelClosed (and reaction names in the model??)
    else:
        rxnNames = []

    modelexchangeAbbr = list(set(modelexchangeAbbr).union(rxnNames))
    TestRxnNum = len(modelexchangeAbbr)
    FluxExV = []

    ###-----------------------------------------------------original-----------------------------------------------------###
    # sets all reactions in modelexchangeAbbr as the objective at once
    # while modelexchangeAbbr:
    #     modelClosed.objective = Zero  # Reset objective
    #     for rxn_id in modelexchangeAbbr:
    #         rxn = modelClosed.reactions.get_by_id(rxn_id)
    #         rxn.objective_coefficient = 1
    #     FF2 = modelClosed.optimize()
    #
    #     print(f"Testing exchange: {modelexchangeAbbr[0]}")
    #     print(f"Objective value: {FF2.objective_value}")
    #     print(f"Remaining: {len(modelexchangeAbbr)}")
    #
    #     if FF2.status == 'infeasible':
    #         print('Trivial solution is not a solution of the mode.\n Check that you are not enforcing flux as '
    #               'Leak testing does not work with force fluxes.')
    #     elif FF2.status != 'optimal':
    #         print('Problems exist in the mode, which lead to the trivial problem being unbounded or otherwise '
    #               'problematic.\n If unbounded, one option could be to reduce the maximal upper/lower bounds to a specific value.')
    #     if FF2.objective_value >= tol:
    #         # Get reactions with flux > tol
    #         FluxR = [rxn.id for rxn in modelClosed.reactions if abs(FF2.fluxes[rxn.id]) > tol]
    #         # Get exchange reactions among them
    #         active_exchanges = list(set(modelexchangeAbbr).intersection(FluxR))
    #         # Append to FluxEx
    #         FluxEx.extend(active_exchanges)
    #         # Append corresponding flux values
    #         FluxExV.extend([FF2.fluxes[rxn_id] for rxn_id in active_exchanges])
    #         # Remove identified ones from modelexchangeAbbr
    #         modelexchangeAbbr = list(set(modelexchangeAbbr) - set(active_exchanges))
    #     else:
    #         break
    ###-----------------------------------------------------original-----------------------------------------------------###

    # set only the first reaction on the list as the objective (like how it behaves in the original MATLAB code)
    while modelexchangeAbbr:
        rxn_id = modelexchangeAbbr[0]  # Take the first one
        modelClosed.objective = modelClosed.reactions.get_by_id(rxn_id)  # Set it as objective

        FF2 = modelClosed.optimize()
        print(f"Testing exchange: {rxn_id}")
        print(f"Objective value: {FF2.objective_value}")
        print(f"Remaining: {len(modelexchangeAbbr)}")

        if FF2.status == "infeasible":
            print("Trivial solution is not a solution of the model.\n Check that you are not enforcing flux...")
        elif FF2.status != "optimal":
            print("Model has issues (e.g., unbounded). Try adjusting bounds.")

        if FF2.objective_value >= tol:
            FluxEx.append(rxn_id)
            FluxExV.append(FF2.objective_value)

        # Remove from the list whether it leaked or not (just like MATLAB logic)
        modelexchangeAbbr.pop(0)

    LeakMets = FluxEx

    return [LeakMets, modelClosed, FluxExV]


def add_demand_reaction(model, metabolite, demand_id=None):
    demand_id = demand_id or f"DM_{metabolite.id}"

    # Check if demand reaction already exists
    if demand_id in model.reactions:
        print(f"Ignoring reaction '{demand_id}' since it already exists.")
        return None

    # Create a new reaction
    demand_rxn = cobra.Reaction(demand_id or f"DM_{metabolite.id}")
    demand_rxn.name = f"Demand for {metabolite.id}"
    demand_rxn.lower_bound = 0
    demand_rxn.upper_bound = 1000  # Or whatever upper limit you want

    # Add stoichiometry: metabolite → ∅ (consume 1 unit)
    demand_rxn.add_metabolites({metabolite: -1})

    # Add to model
    model.add_reactions([demand_rxn])
    return demand_rxn
