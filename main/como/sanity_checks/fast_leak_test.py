"""Perform a Fast Leak Test.

This is a fast leak test originally written in MATLAB.
Source: https://github.com/opencobra/cobratoolbox/blob/950f40cff4256d738fc8421a6ecfe91e084f46dd/src/reconstruction/modelGeneration/fastLeakTest.m

Input
    model: model structure
    testRxns: List of exchange reactions to be tested for leaks
    demandTest (optional): if 'TRUE' is entered, demand reactions for all metabolites in the model are created

Output
    LeakMets: List of exchange reactions for leaking metabolites
    modelClosed: Model structure that has been tested for leaks
    FLuxExV: Flux vector for computed exchange reactions in the closed model
"""

from collections.abc import Iterable

import cobra
import numpy as np
from cobra.util.array import create_stoichiometric_matrix
from loguru import logger
from scipy.sparse import csc_matrix


def fast_leak_test(
    model: cobra.Model,
    test_reactions: Iterable[str],
    demand_test,
    tolerance: float = 1e-6,
):
    """Run a Fast Leak Test.

    This is a fast leak test originally written in MATLAB.
    Source: https://github.com/opencobra/cobratoolbox/blob/950f40cff4256d738fc8421a6ecfe91e084f46dd/src/reconstruction/modelGeneration/fastLeakTest.m

    Args:
        model: The model to test
        test_reactions: List of exchange reactions to be tested for leaks
        demand_test: If 'TRUE', demand reactions for all metabolites in the model are created
        tolerance: Tolerance level for determining significant flux

    """
    model_closed = model.copy()

    # Get stoichiometric matrix in CSC format
    # Step 1: Get a dense NumPy matrix
    s_dense = create_stoichiometric_matrix(model_closed)
    # Step 2: Convert it to a sparse CSC matrix
    s_sparse = csc_matrix(s_dense)

    rxns: list[cobra.Reaction] = list(model_closed.reactions)
    num_rxns = len(rxns)
    # Count non-zero stoichiometric entries per reaction (i.e., columns)
    num_nonzero_per_rxn = np.array((s_sparse != 0).sum(axis=0)).flatten()

    # Generate boolean masks (same as MATLAB logic)
    exp = np.array([(num_nonzero_per_rxn[i] == 1 and rxns[i].lower_bound >= 0) for i in range(num_rxns)])
    upt = np.array([(num_nonzero_per_rxn[i] == 1 and rxns[i].upper_bound <= 0) for i in range(num_rxns)])
    reversible_exchanges = np.array([(num_nonzero_per_rxn[i] == 1 and rxns[i].lower_bound < 0 < rxns[i].upper_bound) for i in range(num_rxns)])

    count = exp | upt | reversible_exchanges
    # Get the reactions where count[i] is True
    exchange_rxns = [rxn for i, rxn in enumerate(model_closed.reactions) if count[i]]
    # Set lower_bound = 0 for those reactions
    for i in range(num_rxns):
        if count[i]:
            rxns[i].lower_bound = 0
    model_exchange_abbr = list(set(test_reactions).union([rxn.id for rxn in exchange_rxns]))

    # test for all demand reactions is an option
    if demand_test:  # add demand reactions for all metabolites in model to check for those too
        rxn_names = []
        for metabolite in model_closed.metabolites:
            demand_rxn = add_demand_reaction(model_closed, metabolite)
            if demand_rxn is not None:
                rxn_names.append(demand_rxn.id)  # collect demand reaction IDs
        # return model_closed (and reaction names in the model??)
    else:
        rxn_names = []

    model_exchange_abbr = list(set(model_exchange_abbr).union(rxn_names))
    flux_ex_v = []

    # set only the first reaction on the list as the objective (like how it behaves in the original MATLAB code)
    flux_exchange = []
    while model_exchange_abbr:
        rxn_id = model_exchange_abbr[0]  # Take the first one
        model_closed.objective = model_closed.reactions.get_by_id(rxn_id)  # Set it as objective

        solution: cobra.Solution = model_closed.optimize()
        logger.debug(f"Testing exchange: {rxn_id}")
        logger.debug(f"Objective value: {solution.objective_value}")
        logger.debug(f"Remaining: {len(model_exchange_abbr)}")

        if solution.status == "infeasible":
            logger.info("Trivial solution is not a solution of the model.\n Check that you are not enforcing flux...")
        elif solution.status != "optimal":
            logger.info("Model has issues (e.g., unbounded). Try adjusting bounds.")

        if solution.objective_value >= tolerance:
            flux_exchange.append(rxn_id)
            flux_ex_v.append(solution.objective_value)

        # Remove from the list whether it leaked or not (just like MATLAB logic)
        model_exchange_abbr.pop(0)

    leak_mets = flux_exchange

    return [leak_mets, model_closed, flux_ex_v]


def add_demand_reaction(
    model: cobra.Model,
    metabolite: cobra.Metabolite,
    demand_id: str | None = None,
    demand_name: str | None = None,
    lower_bound: int = 0,
    upper_bound: int = 1000,
    consumption_rate: int = -1,
) -> cobra.Reaction:
    """Add a demand reaction to a given model.

    Args:
        model: The model to which the demand reaction will be added.
        metabolite: The metabolite for the demand reaction.
        demand_id: Optional ID for the demand reaction. If None, it will be generated as 'DM_{metabolite.id}'.
        demand_name: Optional name for the demand reaction. If None, it will be 'Demand for {metabolite.id}'.
        lower_bound: Lower bound for the demand reaction.
        upper_bound: Upper bound for the demand reaction.
        consumption_rate: The rate at which the metabolite is consumed.

    Returns:
        The demand reaction from the model (if it exists), or the newly created demand reaction.

    """
    demand_id: str = demand_id or f"DM_{metabolite.id}"
    demand_name: str = demand_name or f"Demand for {metabolite.id}"

    # Check if demand reaction already exists
    if demand_id in model.reactions:
        logger.info(f"Ignoring reaction '{demand_id}' since it already exists.")
        return model.reactions.get_by_id(demand_id)  # type: ignore

    # Create a new reaction
    demand_rxn: cobra.Reaction = cobra.Reaction(
        id=demand_id,
        name=demand_name,
        lower_bound=lower_bound,
        upper_bound=upper_bound,
    )

    demand_rxn.add_metabolites({metabolite: consumption_rate})
    model.add_reactions([demand_rxn])
    return demand_rxn
