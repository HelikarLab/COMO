import csv

import cobra
from cobra.flux_analysis import moma, single_reaction_deletion

# The aim of this unit test is to check if the current COMO version (after transition from R to Python)
# There are seven cell behaviors in total in the unit test


# MEMO: Use to search for reaction ID, name, etc.
# Baseline of the healthy B cell model
# como_paper_model_path = "/Users/satominakamura/Desktop/Dr.Helikar Lab/HealthyNaiveB.xml"  # Path to XML file
# como_paper_model: cobra.Model = cobra.io.read_sbml_model(como_paper_model_path)
# reaction: cobra.Reaction
# # for reaction in como_paper_model.reactions:
# #     if ("pfk") in reaction.id.lower():
# #         print(f'ID: {reaction.id}, Name: {reaction.name}')
# # for metabolite in como_paper_model.:
# #     if "lactate" in metabolite.id.lower():
# #         print(metabolite.name)
# #         break
# for group in como_paper_model.groups:
#     if ("fatty acid oxidation") in group.name.lower():
#         print(group.name)
#         # print(group.members)
#         for reaction in group.members:
#             print(reaction.id)
#         break
# 1. Glucose restriction has no impact on naive B cell function
def test_1(model):
    como_model = model
    with open("/Users/satominakamura/Desktop/Dr.Helikar Lab/data_1.csv", "w") as f:
        field_names = ["lower_bound", "biomass"]
        writer = csv.DictWriter(f, fieldnames=field_names)
        writer.writeheader()

    modified_lb = -1000
    while modified_lb < 0:
        for reaction in como_model.reactions:
            como_model.reactions.get_by_id(reaction.id).lower_bound = -1000
            como_model.reactions.get_by_id(reaction.id).upper_bound = 1000
            if reaction.id == "EX_glc_D_e":
                como_model.reactions.get_by_id(reaction.id).lower_bound = modified_lb

            solution_after: cobra.Solution = como_model.optimize()
            data = [{"lower_bound": modified_lb, "biomass": solution_after.objective_value}]
            with open("/Users/satominakamura/Desktop/Dr.Helikar Lab/data_1.csv", "a") as f:
                writer = csv.DictWriter(f, fieldnames=field_names)
                writer.writerows(data)
            print(data)
            modified_lb += 20
            print(f"lower_bound = {modified_lb}")
            break


# 2. Glucose distributes through G6P, F16BP, G3P, and 3PG. Check if the reactions exist in the model.
#   GP6 involved in PGI, PGMT
#   F16BP involved in GALOR,GALTt,EX_galt
#   G3P involved in TPI,FBA,GAPD
#   3PG involved in DPGase,ACYP,PGK,PGM
def test_2(model):
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


# 3. Naive B cells produce lactate. Disrupt the following reaction:
#    Define a sink reaction for lactose and determine how much metabolite is flowing through the pathways.
def test_3(model):
    lactate_sink = cobra.Reaction(
        id="lac_L_c",
        name="Cytoplasmic Lactate Sink",
        subsystem="c",
    )

    solution = model.optimize()
    print(solution.objective_value)
    model.add_boundary(model.metabolites.get_by_id("lac_L[c]"), type="sink")

    rxn: cobra.Reaction = model.reactions.get_by_id("EX_lac_D[e]")
    rxn.knock_out()

    model.objective = "SK_lac_L[c]"
    solution_after = model.optimize()
    print(solution_after.objective_value)


# 4. In naive B cells, citrate, glutamate, glutamine, å-ketoglutarate, succinate, fumarate, malate are active
# Confirm with escher map
def test_4(model):
    pass


# 5. Untreated or anti-IgM treated B cells showed no significant increase in the ECAR in response to the glucose addition
#    Define a sink reaction for lactose and maximize the exchange reaction of glucose. If the model agrees with the
#    verification, the lactose present in the pathway should increase as well
def test_5(model):
    lactate_sink = cobra.Reaction(
        id="lac_L_c",
        name="Cytoplasmic Lactate Sink",
        subsystem="c",
    )

    model.add_boundary(model.metabolites.get_by_id("lac_L[c]"), type="sink")
    rxn: cobra.Reaction = model.reactions.get_by_id("EX_lac_D[e]")
    rxn.knock_out()

    solution_before = model.optimize()
    print("Before", solution_before.objective_value)

    glucose_rxn: cobra.Reaction = model.reactions.get_by_id("EX_glc_D[e]")
    glucose_rxn.lower_bound = -1000
    glucose_rxn.upper_bound = 1000

    model.objective = "SK_lac_L[c]"
    solution_after = model.optimize()
    print("After", solution_after.objective_value)


# 6. B cell survival requires proper maintenance of glycolytic and oxidative glucose pathways.
#    Use the moma method instead of flux balance analysis and delete one reaction at a time.
#    The result reflects the change in the cell's growth rate. The biomass does not change.
def test_6(model):
    model.objective = "biomass_maintenance"
    solution_before = model.optimize()

    single_rxn_del = single_reaction_deletion(model, ["GAPD", "PGI", "HEX1"], method="moma", solution=solution_before)
    # single_rxn_del = single_reaction_deletion(model, ["GAPD", "PGI", "HEX1"], method="fba", solution=solution_before)
    print(single_rxn_del)

    with model as model_copy:
        model_copy.reactions.get_by_id("GAPD").knock_out()
        solution_after: cobra.Solution = model_copy.optimize()
        moma_result: cobra.Solution = moma(model_copy, solution_before, linear=False)
        print("MOMA:", moma_result.objective_value)


# 7. Naive B cells seem to rely on FA as their main fuel
def test_7(model):
    model.objective = "biomass_maintenance"
    solution_before = model.optimize()
    with model as model_copy:
        for group in model_copy.groups:
            if ("fatty acid oxidation") in group.name.lower():
                for reaction in group.members:
                    print(reaction.id)
                    model.reactions.get_by_id(reaction.id).lower_bound = 0
                    model.reactions.get_by_id(reaction.id).upper_bound = 0
                    print(reaction.lower_bound)
                    print(reaction.upper_bound)
        solution_after: cobra.Solution = model_copy.optimize()
    print(solution_before.objective_value)
    print(solution_after.objective_value)


if __name__ == "__main__":
    # Baseline of the healthy B cell model
    como_paper_model_path = "/Users/satominakamura/Desktop/Dr.Helikar Lab/HealthyNaiveB.xml"  # Path to XML file
    como_paper_model: cobra.Model = cobra.io.read_sbml_model(como_paper_model_path)

    # test_1(como_paper_model)
    # test_2(como_paper_model)
    # test_3(como_paper_model)
    # test_4(como_paper_model)
    # test_5(como_paper_model)
    test_6(como_paper_model)
    # test_7(como_paper_model)

### ------------------------------------------------- MEMO ------------------------------------------------- ###
# print(reaction_after.objective_value)
# if "txas" in reaction.id.lower()
# if "R_TXASr"==reaction.id:
# reaction_after: cobra.Solution = como_paper_model.optimize()
# print(reaction_after.objective_value)
# while modified_lb < 0:
# modified_lb += 100
# como_paper_model.reactions.get_by_id(reaction.id).lower_bound = modified_lb
# reaction_after: cobra.Solution = como_paper_model.optimize()
# print(reaction_after.objective_value)
# data = [{"lower_bound": modified_lb}, {"biomass": reaction_after.objective_value}]
# with open('/Users/satominakamura/Desktop/Dr.Helikar Lab/data.csv', 'a') as f:
# field_names = ['lower_bound', 'biomass']
# writer = csv.DictWriter(f, fieldnames=field_names)
# writer.writerows(data)
# print(f"lower_bound = {modified_lb}")

# print(como_paper_model.reactions())
# como_paper_solution: cobra.Solution = como_paper_model.optimize()
# model_summary = como_paper_model.summary
# print(f'model summary = {model_summary}')
# print(f'Biomass value: {como_paper_solution.objective_value}')
