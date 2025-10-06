# ruff: noqa: T201
"""COMO integration tests.

The aim of these tests are to validate COMO's behavior on file changes.


Seven behaviors exist in this file:
1. Validate that glucose restriction has no impact on B cell function.
2. Glucose distributes through G6P, F16BP, G3P, and 3PG. Check if the reactions exist in the model.
3. Naive B cells produce lactate.
4. In naive B cells, citrate, glutamate, glutamine, å-ketoglutarate, succinate, fumarate, malate are active.
5. Naive B cells showed no significant increase in the ECAR in response to the glucose addition.
6. B cell survival requires proper maintenance of glycolytic and oxidative glucose pathways.
7. Naive B cells seem to rely on FA as their main fuel.
"""

import csv
from pathlib import Path

import cobra
from cobra.flux_analysis import moma, single_reaction_deletion

# The aim of this unit test is to check if the current COMO version (after transition from R to Python)
# There are seven cell behaviors in total in the unit test


def glucose_restriction(model: cobra.Model):
    """Validate that glucose restriction has no impact on B cell function.

    Args:
        model: the model to open

    """
    with Path("glucose_validation.csv").open("w", encoding="utf-8") as o_stream:
        csv_writer = csv.DictWriter(f=o_stream, fieldnames=["lower_bound", "biomass"])
        csv_writer.writeheader()

        modified_lb = -1000
        while modified_lb < 0:
            for reaction in model.reactions:
                model.reactions.get_by_id(reaction.id).lower_bound = -1000  # type: ignore
                model.reactions.get_by_id(reaction.id).upper_bound = 1000  # type: ignore
                if reaction.id == "EX_glc_D_e":
                    model.reactions.get_by_id(reaction.id).lower_bound = modified_lb  # type: ignore

                solution_after: cobra.Solution = model.optimize()
                csv_writer.writerows([{"lower_bound": modified_lb, "biomass": solution_after.objective_value}])
                modified_lb += 20
                print(f"lower_bound = {modified_lb}")
                break


def glucose_distribution(model: cobra.Model):
    """Validate that glucose distributes through G6P, F16BP, G3P, and 3PG."""
    r_g6p = ["pgi", "pgmt"]
    r_f16bp = ["galor", "galtt", "ex_galt"]
    r_g3p = ["tpi", "fba", "gapd"]
    r_3pg = ["dpgase", "acyp", "pgk", "pgm"]
    como_model = model
    glc_dist = [r_g6p, r_f16bp, r_g3p, r_3pg]
    total_reaction = sum(len(x) for x in glc_dist)
    num_reaction = 0
    for reaction in como_model.reactions:
        i = 0
        while i < len(glc_dist):
            for r in glc_dist[i]:
                if r in reaction.id.lower():
                    num_reaction += 1
                    glc_dist[i].remove(r)
                    # print(f'{r} found in model')
                else:
                    continue
            i += 1
    if num_reaction == total_reaction:
        print("All reactions found in model")
    else:
        print("Not all reactions found in model")


def lactate_production(model: cobra.Model):
    """Validate that naive B cells produce lactate.

    Naive B cells produce lactate. Disrupt the following reaction:
    Define a sink reaction for lactose and determine how much metabolite is flowing through the pathways.
    """
    with model as model_copy:
        initial_solution = model_copy.optimize()
        print(initial_solution.objective_value)

        model_copy.add_boundary(metabolite=model_copy.metabolites.get_by_id("lac_L"), type="sink")  # type: ignore
        rxn: cobra.Reaction = model.reactions.get_by_id("EX_lac_D[e]")  # type: ignore
        rxn.knock_out()

        model.objective = "SK_lac_L[c]"
        solution_after = model.optimize()
        print(solution_after.objective_value)


def citric_acid_cycle_activation(model: cobra.Model):
    """Validate that Naive B cells activate the citric acid cycle.

    In naive B cells, citrate, glutamate, glutamine, å-ketoglutarate, succinate, fumarate, malate are active
    Confirm with escher map
    """
    metabolites = [
        "cit_c",
        "glu_L_c",
        "gln_L_c",
        "akg_c",
        "succ_c",
        "fum_c",
        "mal_L_c",
    ]
    for met in metabolites:
        reaction_list: list[str] = []

        for reaction in model.reactions:
            reaction: cobra.Reaction

            for metabolite in reaction.metabolites:
                metabolite: cobra.Metabolite

                if metabolite.id == met:
                    reaction_list.append(reaction.id)

        print(f"Metabolite: {met}, Reactions: {reaction_list}")


def ecar_response_with_glucose_addition(model: cobra.Model):
    """Validate that Naive B cells show no significant response in Extracellular Acidification Rate with glucose addition.

    Untreated or anti-IgM treated B cells showed no significant increase in the ECAR in response to the glucose addition
    Define a sink reaction for lactose and maximize the exchange reaction of glucose. If the model agrees with the
    verification, the lactose present in the pathway should increase as well
    """
    with model as model_copy:
        model_copy.add_boundary(metabolite=model_copy.metabolites.get_by_id("lac_L"), type="sink")  # type: ignore

        model_copy.add_boundary(model_copy.metabolites.get_by_id("lac_L[c]"), type="sink")  # type: ignore
        rxn: cobra.Reaction = model_copy.reactions.get_by_id("EX_lac_D[e]")  # type: ignore
        rxn.knock_out()

        solution_before = model_copy.optimize()
        print("Before", solution_before.objective_value)

        glucose_rxn: cobra.Reaction = model_copy.reactions.get_by_id("EX_glc_D[e]")  # type: ignore
        glucose_rxn.lower_bound = -1000
        glucose_rxn.upper_bound = 1000

        model_copy.objective = "SK_lac_L[c]"
        solution_after = model_copy.optimize()
        print("After", solution_after.objective_value)


def glycolysis_and_oxidative_glucose_pathway_maintenance(model: cobra.Model):
    """Validate that Naive B cells require glycolytic and oxidative glucose pathway maintenance.

    B cell survival requires proper maintenance of glycolytic and oxidative glucose pathways.
    Use the moma method instead of flux balance analysis and delete one reaction at a time.
    The result reflects the change in the cell's growth rate. The biomass does not change.
    """
    model.objective = "biomass_maintenance"
    solution_before = model.optimize()

    single_rxn_del = single_reaction_deletion(model, ["GAPD", "PGI", "HEX1"], method="moma", solution=solution_before)
    print(single_rxn_del)

    with model as model_copy:
        model_copy.reactions.get_by_id("GAPD").knock_out()  # type: ignore
        moma_result: cobra.Solution = moma(model_copy, solution_before, linear=False)
        print("MOMA:", moma_result.objective_value)


def fatty_acid_fuel_source(model: cobra.Model):
    """Validate that naive B cells rely on fatty acids as their main fuel source."""
    model.objective = "biomass_maintenance"
    solution_before = model.optimize()
    with model as model_copy:
        for group in model_copy.groups:
            if "fatty acid oxidation" in group.name.lower():
                for reaction in group.members:
                    print(reaction.id)
                    model.reactions.get_by_id(reaction.id).lower_bound = 0  # type: ignore
                    model.reactions.get_by_id(reaction.id).upper_bound = 0  # type: ignore
                    print(reaction.lower_bound)
                    print(reaction.upper_bound)
        solution_after: cobra.Solution = model_copy.optimize()
    print(solution_before.objective_value)
    print(solution_after.objective_value)


if __name__ == "__main__":
    como_paper_model_path = "/Users/satominakamura/Desktop/Dr.Helikar Lab/HealthyNaiveB.xml"  # Path to XML file
    como_paper_model: cobra.Model = cobra.io.read_sbml_model(como_paper_model_path)

    glucose_restriction(model=como_paper_model)
    glucose_distribution(model=como_paper_model)
    lactate_production(model=como_paper_model)
    citric_acid_cycle_activation(model=como_paper_model)
    ecar_response_with_glucose_addition(model=como_paper_model)
    glycolysis_and_oxidative_glucose_pathway_maintenance(model=como_paper_model)
    fatty_acid_fuel_source(model=como_paper_model)
