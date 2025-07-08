# This is a sanity check for metabolic model implemented in COBRA Toolbox.
# The original code is available at https://opencobra.github.io/cobratoolbox/stable/tutorials/tutorial_modelSanityChecks.html
# written in MATLAB. This is the translated version written in Python.

# Content
# 1. leak test
# 2. production of protons from nothing as well as from water, and/or oxygen alone
# 3. production of matter when atp hydrolysis reaction is allowed to work but all uptakes are closed
# 4. production of too much ATP from glucose under aerobic condition
# 5. duplicated reactions
# 6. empty columns in the model.rxnGeneMat
# 7. the single gene deletion analysis runs smoothly
# 8. ATP yield fro different carbon sources
# 9. metabolic objective functions
# 10. flux consistency
# 11. demand reactions with negative lower bound (should not occur based on definition of demand reactions)
# 12. consistency of model.rev which defines reaction reversibility, and the set values for the lower bounds on reactions

import os
import re
import sys
import warnings

import cobra
import numpy as np
from cobra import Metabolite, Model, Reaction
from cobra.flux_analysis import single_gene_deletion
from cobra.util.array import create_stoichiometric_matrix
from scipy.sparse import csc_matrix

from como.sanity_checks import fastLeakTest


def preprocess():
    ## EQUIPMENT SETUP ##
    # Load data
    model_path = "/Users/satominakamura/Desktop/Dr.Helikar Lab/HealthyNaiveB.xml"  # Path to XML file
    model: cobra.Model = cobra.io.read_sbml_model(model_path)

    # Change solver if necessary
    model.solver = "glpk"  # solvers available: 'glpk', 'ibm_cplex','tomlab_cplex'

    ## PROCEDURE ##
    ####-------------------------------------------------------OK------------------------------------------------------####
    ## MODEL HARMONIZATION ##
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

    # === Add DM_atp[c] if still missing ===
    if "DM_atp[c]" not in model.reactions:
        print("Adding DM_atp[c] manually.")

        # Make sure ATP, ADP, H2O, H+ and Pi exist or are added
        def ensure_met(met_id, name, compartment):
            if met_id in model.metabolites:
                return model.metabolites.get_by_id(met_id)
            else:
                met = Metabolite(id=met_id, name=name, compartment=compartment)
                model.add_metabolites([met])
                return met

        atp = ensure_met("atp_c", "ATP", "c")
        adp = ensure_met("adp_c", "ADP", "c")
        pi = ensure_met("pi_c", "Phosphate", "c")
        h2o = ensure_met("h2o_c", "H2O", "c")
        h = ensure_met("h_c", "Proton", "c")

        atp_hydrolysis = Reaction("DM_atp[c]")
        atp_hydrolysis.name = "ATP maintenance reaction"
        atp_hydrolysis.lower_bound = 0
        atp_hydrolysis.upper_bound = 1000
        atp_hydrolysis.add_metabolites({adp: 1.0, atp: -1.0, h: 1.0, h2o: -1.0, pi: 1.0})
        model.add_reactions([atp_hydrolysis])

    # === Add water exchange reaction if missing ===
    if "EX_h2o[e]" not in model.reactions:
        print("Adding EX_h2o[e] manually.")
        if "h2o_e" not in model.metabolites:
            h2o_e = Metabolite("h2o_e", name="H2O", compartment="e")
            model.add_metabolites([h2o_e])
        else:
            h2o_e = model.metabolites.get_by_id("h2o_e")

        model.add_boundary(h2o_e, type="exchange", reaction_id="EX_h2o[e]")

    if "DM_atp[c]" in model.reactions:
        print("Confirmed: DM_atp[c] is now in the model.")
    else:
        print("Warning: DM_atp[c] is still missing after harmonization.")

    # Set lower bound of the biomass reaction to 0.
    if "biomass_maintenance" in model.reactions:
        model.reactions.get_by_id("biomass_maintenance").lower_bound = 0
    print("Finished setting lower bound of the biomass reactions to 0")

    # Harmonize different use of brackets.
    for rxn in model.reactions:
        new_id = rxn.id
        # Apply regex substitutions to clean up the ID
        new_id = re.sub(r"\(", "[", new_id)  # Replace "(" with "["
        new_id = re.sub(r"\)", "]", new_id)  # Replace ")" with "]"
        new_id = re.sub(r"Ex_", "EX_", new_id)  # Replace "Ex_" at beginning with "EX_"
        new_id = re.sub(r"Sink_", "sink_", new_id)  # Replace "Sink_" at beginning with "sink_"
        new_id = re.sub(r"-", "_", new_id)  # Replace "-" with "_"

        # Assign the new ID back to the reaction
        rxn.id = new_id
    print("Finished performing regex substitution")

    # Define some parameters that we will need.
    cnt = 1
    tol = 1e-6
    ####-------------------------------------------------------OK------------------------------------------------------####
    # Define the closed model
    model_closed = model.copy()

    # Step 1: Collect reaction IDs by regex name patterns
    modelexchanges_ids = []

    # match 'Ex_', 'EX_', 'DM_', or 'sink_' in the reaction id and add matched reaction ids in a list
    for rxn in model_closed.reactions:
        rxn_id = rxn.id
        if re.search(r"\bEx_", rxn_id) or re.search(r"\bEX_", rxn_id) or re.search(r"\bDM_", rxn_id) or re.search(r"\bsink_", rxn_id):
            modelexchanges_ids.append(rxn_id)

    # Step 2: Matrix-based method to find selExc reactions (single metabolite, 1 entry)

    # Step 1: Get a dense NumPy matrix
    S_dense = create_stoichiometric_matrix(model_closed)

    # Step 2: Convert it to a sparse CSC matrix
    S = csc_matrix(S_dense)
    abs_S = abs(S)
    # column_counts = (abs_S != 0).sum(axis=0)

    selExc = [model_closed.reactions[j].id for j in range(S.shape[1]) if np.sum(abs_S[:, j] == 1) == 1 and np.sum(S[:, j] != 0) == 1]

    # Step 3: Combine and deduplicate all identified exchange-like reactions
    modelexchanges = sorted(set(modelexchanges_ids + selExc))

    # Step 4: Set bounds for these reactions
    for rxn_id in modelexchanges:
        try:
            model_closed.reactions.get_by_id(rxn_id).bounds = (0, 1000)
        except KeyError:
            print(f"Warning: Reaction {rxn_id} not found in model.")

    # Optional: Save original model copy
    model_closed_original = model_closed

    print("Success in model preprocessing!")
    return model_closed_original, selExc


def main(model_closed, selExc):
    TableChecks = []
    tol = 1e-6

    ## START WITH TESTS ##

    # 1. leak test
    # [SUCCESS] Perform leak test, i.e., whether the closed model can produce any exchanged metabolite, as defined in the model, from nothing.
    with model_closed as model_copy:
        LeakRxns, modelTested, LeakRxnsFluxVector = fastLeakTest.fast_leak_test_function(
            model_copy, [rxn.id for rxn in model_copy.reactions if rxn.id in selExc], demandTest=False
        )
        table_check_row = ["fastLeakTest1"]
        if len(LeakRxns) > 0:
            warnings.warn("model leaks metabolites!")
            table_check_row.append("Model leaks metabolites!")
        else:
            table_check_row.append("Leak free!")
        TableChecks.append(table_check_row)

    ##  [SUCCESS] Test if something leaks when demand reactions for each metabolite in the model are added. Note that this step is time consuming.
    with model_closed as model_copy:
        TableChecks = []
        LeakRxnsDM, modelTestedDM, LeakRxnsFluxVectorDM = fastLeakTest.fast_leak_test_function(
            model_copy, [rxn.id for rxn in model_copy.reactions if rxn.id in selExc], demandTest=True
        )
        table_check_row = ["fastLeakTest 2 - add demand reactions for each metabolite in the model"]
        if len(LeakRxnsDM) > 0:
            table_check_row.append("Model leaks metabolites when demand reactions are added!")
        else:
            table_check_row.append("Leak free when demand reactions are added!")
        TableChecks.append(table_check_row)

    ## [SUCCESS] Test if the model produces energy from water!
    with model_closed as model_copy:
        TableChecks = []
        model_copy.objective = "DM_atp[c]"
        model_copy.reactions.get_by_id("DM_atp[c]").lower_bound = 0
        model_copy.reactions.get_by_id("EX_h2o[e]").lower_bound = -1

        FBA3 = model_copy.optimize()
        table_check_row = ["EExchanges, sinks, and demands have lb = 0, except h2o"]
        if abs(FBA3.objective_value) > tol:
            table_check_row.append("model produces energy from water!")
        else:
            table_check_row.append("model DOES NOT produce energy from water!")
        TableChecks.append(table_check_row)

    ## [SUCCESS] Test if the model produces energy from water and oxygen!
    with model_closed as model_copy:
        TableChecks = []
        model_copy.objective = "DM_atp[c]"
        model_copy.reactions.get_by_id("DM_atp[c]").lower_bound = 0
        model_copy.reactions.get_by_id("EX_h2o[e]").lower_bound = -1  # allowing water to come in
        model_copy.reactions.get_by_id("EX_o2[e]").lower_bound = -1  # allowing oxygen to come in

        FBA6 = model_copy.optimize()
        table_check_row = ["Exchanges,sinks, and demands have lb = 0, except h2o and o2"]
        if abs(FBA6.objective_value) > tol:
            table_check_row.append("model produces energy from water and oxygen!")
        else:
            table_check_row.append("model DOES NOT produce energy from water and oxygen!")
        TableChecks.append(table_check_row)

    # ## [SUCCESS] Test if the model produces matter when atp demand is reversed!
    with model_closed as model_copy:
        model_copy.objective = "DM_atp[c]"
        model_copy.reactions.get_by_id("DM_atp[c]").lower_bound = -1000
        FBA = model_copy.optimize()
        table_check_row = ["Exchanges, sinks, and demands have lb = 0, allow DM_atp_c_ to be reversible"]
        if abs(FBA.objective_value) > tol:
            table_check_row.append("model produces no matter when atp demand is reversed!")
        else:
            table_check_row.append("model DOES NOT produce no matter when atp demand is reversed!")
        TableChecks.append(table_check_row)

    # ## [SUCCESS] Test if the model has flux through h[m] demand !
    with model_closed as model_copy:
        model_copy.add_reactions(
            [
                cobra.Reaction(id="DM_h[m]", name="DM_h[m]", upper_bound=1000)  # lower bound not specified in original code
            ]
        )
        FBA = model_copy.optimize()
        table_check_row = ["Exchanges, sinks, and demands have lb = 0, test flux through DM_h[m] (max)"]
        if abs(FBA.objective_value) > tol:
            table_check_row.append("model has flux through DM_h[m] (max)!")
        else:
            table_check_row.append("model has NO flux through DM_h[m] (max)!")
        TableChecks.append(table_check_row)

    # ## [SUCCESS] Test if the  model has flux through h[c] demand !
    with model_closed as model_copy:
        model_copy.add_reactions(
            [
                cobra.Reaction(id="DM_h[c]", name="DM_h[c]", upper_bound=1000)  # lower bound not specified in original code
            ]
        )
        FBA = model_copy.optimize()
        table_check_row = ["Exchanges, sinks, and demands have lb = 0, test flux through DM_h[c] (max)"]
        if abs(FBA.objective_value) > tol:
            table_check_row.append("model has flux through DM_h[c] (max)!")
        else:
            table_check_row.append("model has NO flux through DM_h[c] (max)!")
        TableChecks.append(table_check_row)

    # ## [SUCCESS] Test if the  model produces too much atp demand from glucose under aerobic condition.
    # # Also consider using the tutorial testModelATPYield to test if the correct ATP yield from different carbon sources
    # # can be realized by the model.
    with model_closed as model_copy:
        model_copy.objective = "DM_atp[c]"
        model_copy.reactions.get_by_id("EX_o2[e]").lower_bound = -1000
        model_copy.reactions.get_by_id("EX_h2o[e]").lower_bound = -1000
        model_copy.reactions.get_by_id("EX_h2o[e]").upper_bound = 1000
        model_copy.reactions.get_by_id("EX_co2[e]").upper_bound = 1000
        FBAOri = model_closed.optimize()

        table_check_row = ["ATP yield"]
        if abs(FBAOri.objective_value) > 31:  # this is the theoretical value
            table_check_row.append("model produces too much atp demand from glc!")
        else:
            table_check_row.append("model DOES NOT produce too much atp demand from glc!")
        TableChecks.append(table_check_row)

    ## [SKIP] Test metabolic objective functions with open sinks. Note this step is time consuming and may only work.

    ## [SKIP] Test metabolic objective functions with closed sinks (lb). Note this step is time consuming and may only work.

    ## [SKIP] Compute ATP yield. This test is identical to the material covered in the tutorial testModelATPYield.

    ## [SUCCESS] Check for duplicated reactions in the model.
    # # Method: look for duplicated reaction ID in the model. (Original approach)
    with model_closed as model_copy:
        table_check_row = ["Check duplicated reactions"]
        print("Checking for duplicated reactions...")
        all_reactions = model_copy.reactions
        for rxn1 in range(len(all_reactions)):  # print only first 20 for brevity
            for rxn2 in range(rxn1 + 1, len(all_reactions)):
                if all_reactions[rxn1].id == all_reactions[rxn2].id:
                    table_check_row.append(f"There exists a duplicated reaction: {all_reactions[rxn1].id}")
        TableChecks.append(table_check_row)
        print("Finished checking for duplicated reactions!")

    #
    ## [SUCCESS] (check with Josh) Check empty columns in 'model.genes'.
    with model_closed as model_copy:
        table_check_row = ["Check empty columns in model_copy.genes"]
        all_genes = model_copy.genes
        for gene in all_genes:
            if gene is None:
                table_check_row.append("There is an empty column in model_copy.genes")
        TableChecks.append(table_check_row)

    ## [SUCCESS] Check that demand reactions have a lb >= 0.
    with model_closed as model_copy:
        table_check_row = ["Check that demand reactions have lb >= 0"]
        for rxn in model_copy.reactions:
            reaction_id = rxn.id
            if (re.search(r"\bDM_", reaction_id)) and rxn.lower_bound < 0:
                table_check_row.append(f"Demand reaction can have flux in backward direction. {rxn.id}")
        TableChecks.append(table_check_row)

    ## [SUCCESS] Check whether singleGeneDeletion runs smoothly.
    with model_closed as model_copy:
        table_check_row = ["Check whether singleGeneDeletion runs smoothly"]
        try:
            single_gene_deletion(model_copy)
            table_check_row.append("singleGeneDeletion finished without problems")
        except:
            table_check_row.append("There are problems with singleGeneDeletion.")
        TableChecks.append(table_check_row)

    ## [CODE WORKS BUT INCONSISTENCY IN CURRENT MODEL] Check for flux consistency.
    with model_closed as model_copy:
        table_check_row = ["Check for flux consistency"]
        # fastcc requires a numeric stoichiometric matrix and reaction list
        consistent_submodel = cobra.flux_analysis.fastcc(model_copy, 1e-4)

        if len(consistent_submodel.reactions) == len(model_copy.reactions):
            table_check_row.append("Model is flux consistent.")
        else:
            table_check_row.append("Model is NOT flux consistent.")

        TableChecks.append(table_check_row)
        for row in TableChecks:
            if len(row) >= 2:
                print(f"{row[0]:<60} {row[1]}")
            else:
                print(row)
        original_rxns = set(r.id for r in model_copy.reactions)
        consistent_rxns = set(r.id for r in consistent_submodel.reactions)
        inconsistent_rxns = original_rxns - consistent_rxns

        print(f"Inconsistent reactions ({len(inconsistent_rxns)}):")
        print(f"original reactions: {len(original_rxns)}")
        # for rxn_id in sorted(inconsistent_rxns):
        #     print(rxn_id)


if __name__ == "__main__":
    # Baseline of the healthy B cell model
    # model_path = "/Users/satominakamura/Desktop/Dr.Helikar Lab/Recon3DModel_301.mat"  # Path to mat file (healthy human metabolism) from Recon3D
    # model: cobra.Model = cobra.io.read_sbml_model(model_path)

    model_closed, selExc = preprocess()
    main(model_closed, selExc)
