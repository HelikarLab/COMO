"""Sanity check from COBRA Toolbox.

The original code is available at: https://opencobra.github.io/cobratoolbox/stable/tutorials/tutorial_modelSanityChecks.html
Initially written in MATLAB, this is the translated version written in Python.

Checks performed:
1. leak test
2. production of protons from nothing as well as from water, and/or oxygen alone
3. production of matter when atp hydrolysis reaction is allowed to work but all uptakes are closed
4. production of too much ATP from glucose under aerobic condition
5. duplicated reactions
6. empty columns in the model.rxnGeneMat
7. the single gene deletion analysis runs smoothly
8. ATP yield fro different carbon sources
9. metabolic objective functions
10. flux consistency
11. demand reactions with negative lower bound (should not occur based on definition of demand reactions)
12. consistency of model.rev which defines reaction reversibility, and the set values for the lower bounds on reactions
"""

import re
from pathlib import Path
from typing import Literal

import cobra
import cobra.util.array
import numpy as np
import numpy.typing as npt
from loguru import logger
from scipy.sparse import csc_matrix

from como.sanity_checks.fast_leak_test import fast_leak_test

# from cobra import Metabolite, Reaction
# from cobra.flux_analysis import single_gene_deletion


def create_metabolite_if_not_exists(model: cobra.Model, met_id: str, name: str, compartment: str) -> cobra.Metabolite:
    """Validate that a metabolite exists in a given model.

    Args:
        model: The model to check
        met_id: The metabolite ID to check
        name: The name of the metabolite to add if it does not exist
        compartment: The compartment of the metabolite to add if it does not exist

    Returns:
        The metabolite object, either existing or newly created.

    """
    if met_id in model.metabolites:
        return model.metabolites.get_by_id(met_id)  # type: ignore
    else:
        met = cobra.Metabolite(id=met_id, name=name, compartment=compartment)
        model.add_metabolites([met])
        return met


def preprocess(model_path: Path, solver: Literal["glpk", "cplex", "gurobi"] = "glpk") -> tuple[cobra.Model, list[str]]:  # noqa: C901
    """Close the given model.

    Args:
        model_path: Path to the SBML model file
        solver: The solver to use. Options are "glpk", "cplex", or "gurobi". Default is "glpk".

    Returns:
        A tuple containing the closed model and a list of selected exchange reactions.

    Raises:
        FileNotFoundError: If the model file does not exist.

    """
    if not model_path.exists():
        raise FileNotFoundError(f"File {model_path} does not exist")

    model: cobra.Model = cobra.io.read_sbml_model(model_path)

    # Change solver if necessary
    model.solver = solver

    #  Replace reaction abbreviation for the ATP hydrolysis (DM_atp_c_) and Biomass reaction used differently in various models.
    if "ATPM" in model.reactions:
        model.reactions.get_by_id("ATPM").id = "DM_atp[c]"
    if "ATPhyd" in model.reactions:
        model.reactions.get_by_id("ATPhyd").id = "DM_atp[c]"
    if "DM_atp(c)" in model.reactions:
        model.reactions.get_by_id("DM_atp(c)").id = "DM_atp[c]"
    if "EX_biomass_reaction" in model.reactions:
        model.reactions.get_by_id("EX_biomass_reaction").id = "biomass_reaction"
    if "EX_biomass_maintenance" in model.reactions:
        model.reactions.get_by_id("EX_biomass_maintenance").id = "biomass_maintenance"
    if "EX_biomass_reaction_noTrTr" in model.reactions:
        model.reactions.get_by_id("EX_biomass_reaction_noTrTr").id = "biomass_maintenance_noTrTr"

    if "DM_atp[c]" not in model.reactions:
        # Make sure ATP, ADP, H2O, H+ and Pi exist or are added

        atp: cobra.Metabolite = create_metabolite_if_not_exists(model=model, met_id="atp_c", name="ATP", compartment="c")
        adp: cobra.Metabolite = create_metabolite_if_not_exists(model=model, met_id="adp_c", name="ADP", compartment="c")
        pi: cobra.Metabolite = create_metabolite_if_not_exists(model=model, met_id="pi_c", name="Phosphate", compartment="c")
        h2o: cobra.Metabolite = create_metabolite_if_not_exists(model=model, met_id="h2o_c", name="H2O", compartment="c")
        h: cobra.Metabolite = create_metabolite_if_not_exists(model=model, met_id="h_c", name="Proton", compartment="c")

        atp_hydrolysis = cobra.Reaction("DM_atp[c]")
        atp_hydrolysis.name = "ATP maintenance reaction"
        atp_hydrolysis.lower_bound = 0
        atp_hydrolysis.upper_bound = 1000
        atp_hydrolysis.add_metabolites({adp: 1.0, atp: -1.0, h: 1.0, h2o: -1.0, pi: 1.0})
        model.add_reactions([atp_hydrolysis])

    # === Add water exchange reaction if missing ===
    if "EX_h2o[e]" not in model.reactions:
        h2o_e = create_metabolite_if_not_exists(model=model, met_id="h2o_e", name="H2O", compartment="e")
        model.add_boundary(h2o_e, type="exchange", reaction_id="EX_h2o[e]")

    # Set lower bound of the biomass reaction to 0.
    if "biomass_maintenance" in model.reactions:
        model.reactions.get_by_id("biomass_maintenance").lower_bound = 0  # type: ignore

    # Harmonize different use of brackets.
    for rxn in model.reactions:
        new_id = rxn.id
        new_id = new_id.replace("(", "[")
        new_id = new_id.replace(")", "]")
        new_id = new_id.replace("-", "_")
        new_id = re.sub(r"^Ex_", "EX_", new_id)
        new_id = re.sub(r"^(Sink|sink|SK)_", "sink_", new_id)
        rxn.id = new_id

    model_closed = model.copy()

    # Step 1: Collect reaction IDs by regex name patterns
    # match 'Ex_', 'EX_', 'DM_', or 'sink_' in the reaction id and add matched reaction ids in a list
    modelexchanges_ids: list[str] = [rxn.id for rxn in model_closed.reactions if re.search(pattern=r"^(EX|DM|sink)_", string=rxn.id)]

    # Step 2: Matrix-based method to find selected_exchanges reactions (single metabolite, 1 entry)
    s_matrix: npt.NDArray[np.floating] = cobra.util.array.create_stoichiometric_matrix(model_closed)

    # Step 2: Convert it to a sparse CSC matrix to optimize memory usage and performance
    s_sparse = csc_matrix(s_matrix)
    abs_s = abs(s_sparse)

    selected_exchanges: list[str] = [
        model_closed.reactions[j].id for j in range(s_sparse.shape[1]) if np.sum(abs_s[:, j] == 1) == 1 and np.sum(s_sparse[:, j] != 0) == 1
    ]

    # Step 3: Combine and deduplicate all identified exchange-like reactions
    modelexchanges: list[str] = sorted(set(modelexchanges_ids + selected_exchanges))

    # Step 4: Set bounds for these reactions
    try:
        for rxn_id in modelexchanges:
            model_closed.reactions.get_by_id(rxn_id).bounds = (0, 1000)  # type: ignore
    except KeyError:
        logger.warning(f"Reaction {rxn_id} not found in model.")

    # Optional: Save original model copy
    model_closed_original = model_closed

    logger.success("Model preprocessed!")
    return model_closed_original, selected_exchanges


def run_sanity_checks(model_closed: cobra.Model, selected_exchanges: list[str]):  # noqa: C901
    """Run the MATLAB sanity checks.

    Args:
        model_closed: The closed model from `preprocess`
        selected_exchanges: The selected exchanges from `preprocess`

    """
    table_checks = []
    tol = 1e-6

    # Perform leak test, i.e., whether the closed model can produce any exchanged metabolite, as defined in the model, from nothing.
    with model_closed as model_copy:
        leak_reactions, _, _ = fast_leak_test(
            model_copy,
            [rxn.id for rxn in model_copy.reactions if rxn.id in selected_exchanges],
            demand_test=False,
        )
        table_check_row = ["fastLeakTest1"]
        if len(leak_reactions) > 0:
            logger.warning("Model leaks metabolites!")
            table_check_row.append("Model leaks metabolites!")
        else:
            table_check_row.append("Leak free!")
        table_checks.append(table_check_row)

    # Test if something leaks when demand reactions for each metabolite in the model are added. Note that this step is time consuming.
    with model_closed as model_copy:
        table_checks = []
        leak_rxn_demand, _, _ = fast_leak_test(model_copy, [rxn.id for rxn in model_copy.reactions if rxn.id in selected_exchanges], demand_test=True)
        table_check_row = ["fastLeakTest 2 - add demand reactions for each metabolite in the model"]
        if len(leak_rxn_demand) > 0:
            table_check_row.append("Model leaks metabolites when demand reactions are added!")
        else:
            table_check_row.append("Leak free when demand reactions are added!")
        table_checks.append(table_check_row)

    # Test if the model produces energy from water
    with model_closed as model_copy:
        table_checks = []
        model_copy.objective = "DM_atp[c]"
        model_copy.reactions.get_by_id("DM_atp[c]").lower_bound = 0  # type: ignore
        model_copy.reactions.get_by_id("EX_h2o[e]").lower_bound = -1  # type: ignore

        solution: cobra.Solution = model_copy.optimize()
        table_check_row = ["EExchanges, sinks, and demands have lb = 0, except h2o"]
        if abs(solution.objective_value) > tol:
            table_check_row.append("model produces energy from water!")
        else:
            table_check_row.append("model DOES NOT produce energy from water!")
        table_checks.append(table_check_row)

    # Test if the model produces energy from water and oxygen
    with model_closed as model_copy:
        table_checks = []
        model_copy.objective = "DM_atp[c]"
        model_copy.reactions.get_by_id("DM_atp[c]").lower_bound = 0  # type: ignore
        model_copy.reactions.get_by_id("EX_h2o[e]").lower_bound = -1  # allowing water to come in # type: ignore
        model_copy.reactions.get_by_id("EX_o2[e]").lower_bound = -1  # allowing oxygen to come in # type: ignore

        solution: cobra.Solution = model_copy.optimize()
        table_check_row = ["Exchanges,sinks, and demands have lb = 0, except h2o and o2"]
        if abs(solution.objective_value) > tol:
            table_check_row.append("model produces energy from water and oxygen!")
        else:
            table_check_row.append("model DOES NOT produce energy from water and oxygen!")
        table_checks.append(table_check_row)

    # Test if the model produces matter when atp demand is reversed
    with model_closed as model_copy:
        model_copy.objective = "DM_atp[c]"
        model_copy.reactions.get_by_id("DM_atp[c]").lower_bound = -1000  # type: ignore
        solution = model_copy.optimize()
        table_check_row = ["Exchanges, sinks, and demands have lb = 0, allow DM_atp_c_ to be reversible"]
        if abs(solution.objective_value) > tol:
            table_check_row.append("model produces no matter when atp demand is reversed!")
        else:
            table_check_row.append("model DOES NOT produce no matter when atp demand is reversed!")
        table_checks.append(table_check_row)

    # Test if the model has flux through h[m] demand
    with model_closed as model_copy:
        model_copy.add_reactions(
            [
                cobra.Reaction(id="DM_h[m]", name="DM_h[m]", upper_bound=1000)  # lower bound not specified in original code
            ]
        )
        solution = model_copy.optimize()
        table_check_row = ["Exchanges, sinks, and demands have lb = 0, test flux through DM_h[m] (max)"]
        if abs(solution.objective_value) > tol:
            table_check_row.append("model has flux through DM_h[m] (max)!")
        else:
            table_check_row.append("model has NO flux through DM_h[m] (max)!")
        table_checks.append(table_check_row)

    # Test if the  model has flux through h[c] demand
    with model_closed as model_copy:
        model_copy.add_reactions(
            [
                cobra.Reaction(id="DM_h[c]", name="DM_h[c]", upper_bound=1000)  # lower bound not specified in original code
            ]
        )
        solution: cobra.Solution = model_copy.optimize()
        table_check_row = ["Exchanges, sinks, and demands have lb = 0, test flux through DM_h[c] (max)"]
        if abs(solution.objective_value) > tol:
            table_check_row.append("model has flux through DM_h[c] (max)!")
        else:
            table_check_row.append("model has NO flux through DM_h[c] (max)!")
        table_checks.append(table_check_row)

    # Test if the  model produces too much atp demand from glucose under aerobic condition.
    # Also consider using the tutorial testModelATPYield to test if the correct ATP yield from different carbon sources can be realized by the model.
    with model_closed as model_copy:
        model_copy.objective = "DM_atp[c]"
        model_copy.reactions.get_by_id("EX_o2[e]").lower_bound = -1000  # type: ignore
        model_copy.reactions.get_by_id("EX_h2o[e]").lower_bound = -1000  # type: ignore
        model_copy.reactions.get_by_id("EX_h2o[e]").upper_bound = 1000  # type: ignore
        model_copy.reactions.get_by_id("EX_co2[e]").upper_bound = 1000  # type: ignore
        original_sol: cobra.Solution = model_closed.optimize()

        table_check_row = ["ATP yield"]
        if abs(original_sol.objective_value) > 31:  # this is the theoretical value
            table_check_row.append("model produces too much atp demand from glc!")
        else:
            table_check_row.append("model DOES NOT produce too much atp demand from glc!")
        table_checks.append(table_check_row)

    # Test metabolic objective functions with open sinks. Note this step is time consuming and may only work.

    # Test metabolic objective functions with closed sinks (lb). Note this step is time consuming and may only work.

    # Compute ATP yield. This test is identical to the material covered in the tutorial testModelATPYield.

    # Check for duplicated reactions in the model.
    with model_closed as model_copy:
        table_check_row = ["Check duplicated reactions"]
        logger.debug("Checking for duplicated reactions...")
        all_reactions = model_copy.reactions
        for rxn1 in range(len(all_reactions)):  # print only first 20 for brevity
            for rxn2 in range(rxn1 + 1, len(all_reactions)):
                if all_reactions[rxn1].id == all_reactions[rxn2].id:
                    table_check_row.append(f"There exists a duplicated reaction: {all_reactions[rxn1].id}")
        table_checks.append(table_check_row)
        logger.info("Finished checking for duplicated reactions!")

    # Check empty columns in 'model.genes'.
    with model_closed as model_copy:
        table_check_row = ["Check empty columns in model_copy.genes"]
        all_genes = model_copy.genes
        for gene in all_genes:
            if gene is None:
                table_check_row.append("There is an empty column in model_copy.genes")
        table_checks.append(table_check_row)

    # Check that demand reactions have a lb >= 0.
    with model_closed as model_copy:
        table_check_row = ["Check that demand reactions have lb >= 0"]
        for rxn in model_copy.reactions:
            reaction_id = rxn.id
            if (re.search(r"\bDM_", reaction_id)) and rxn.lower_bound < 0:
                table_check_row.append(f"Demand reaction can have flux in backward direction. {rxn.id}")
        table_checks.append(table_check_row)

    # Check whether singleGeneDeletion runs smoothly.
    with model_closed as model_copy:
        table_check_row = ["Check whether singleGeneDeletion runs smoothly"]
        try:
            cobra.flux_analysis.single_gene_deletion(model_copy)
            table_check_row.append("singleGeneDeletion finished without problems")
        except Exception:
            table_check_row.append("There are problems with singleGeneDeletion.")
        table_checks.append(table_check_row)

    # Check for flux consistency.
    with model_closed as model_copy:
        table_check_row = ["Check for flux consistency"]
        # fastcc requires a numeric stoichiometric matrix and reaction list
        consistent_submodel = cobra.flux_analysis.fastcc(model_copy, 1e-4)

        if len(consistent_submodel.reactions) == len(model_copy.reactions):
            table_check_row.append("Model is flux consistent.")
        else:
            table_check_row.append("Model is NOT flux consistent.")

        table_checks.append(table_check_row)
        for row in table_checks:
            if len(row) >= 2:
                logger.debug(f"{row[0]:<60} {row[1]}")
            else:
                logger.debug(row)
        original_rxns = {r.id for r in model_copy.reactions}
        consistent_rxns = {r.id for r in consistent_submodel.reactions}
        inconsistent_rxns = original_rxns - consistent_rxns

        logger.info(f"Inconsistent reactions ({len(inconsistent_rxns)}):")
        logger.info(f"original reactions: {len(original_rxns)}")


if __name__ == "__main__":
    model_path = Path("/Users/satominakamura/Desktop/Dr.Helikar Lab/HealthyNaiveB.xml")  # Path to XML file

    model_closed, selected_exchanges = preprocess(model_path)
    run_sanity_checks(model_closed, selected_exchanges)
