# ruff: noqa: D103 T201
"""COMO integration tests.

The aim of these tests are to validate COMO's behavior on file changes.


Seven behaviors exist in this file:
1. Validate that glucose restriction has no impact on B cell function.
2. Glucose distributes through G6P, F16BP, G3P, and 3PG. Check if the reactions exist in the model.
3. Naive B cells produce lactate.
4. In naive B cells, citrate, glutamate, glutamine, alpha-ketoglutarate, succinate, fumarate, malate are active.
5. Naive B cells showed no significant increase in the ECAR in response to the glucose addition.
6. B cell survival requires proper maintenance of glycolytic and oxidative glucose pathways.
7. Naive B cells seem to rely on FA as their main fuel.
"""

from collections import defaultdict
from contextlib import redirect_stderr, redirect_stdout, suppress
from pathlib import Path
from typing import cast

import cobra
import numpy as np
import pandas as pd
from cobra import Metabolite, Model, Reaction
from cobra.flux_analysis import single_reaction_deletion

from como.sanity_checks.naiveB_requirement import GLYCOLYSIS_REACTIONS, add_exchanges, add_glycolysis_reactions

# The aim of this unit test is to check if the current COMO version (after transition from R to Python)
# There are seven cell behaviors in total in the unit test


OPEN_LOWER_BOUND: int = -1000
OPEN_UPPER_BOUND: int = 1000
CLOSED_BOUNDARY: int = 0


def constrain_model(model: Model) -> Model:
    cast(Reaction, model.reactions.get_by_id("EX_glc_D[e]")).lower_bound = -1.0  # type: ignore[bad-argument-count]

    if "EX_hdca[e]" in model.reactions:
        cast(Reaction, model.reactions.get_by_id("EX_hdca[e]")).lower_bound = -20.0  # type: ignore[bad-argument-count]
    if "EX_gln_L[e]" in model.reactions:
        cast(Reaction, model.reactions.get_by_id("EX_gln_L[e]")).lower_bound = -5.0  # type: ignore[bad-argument-count]
    if "EX_o2[e]" in model.reactions:
        cast(Reaction, model.reactions.get_by_id("EX_o2[e]")).lower_bound = -1000.0  # type: ignore[bad-argument-count]

    for fuel in ["EX_fru[e]", "EX_sucr[e]", "EX_gal[e]"]:  # close off high-energy carbon sources
        if fuel in model.reactions:
            cast(Reaction, model.reactions.get_by_id(fuel)).lower_bound = 0.0  # type: ignore[bad-argument-count]
    return model


def glucose_restriction(model: Model):
    """Validate that glucose restriction has no impact on B cell function."""
    df = pd.DataFrame({"lower_bound": list(range(-1000, 1, 20)), "biomass": np.nan})
    glucose_rxn = "EX_glc_D[e]"

    for idx in range(len(df)):
        lb = df.at[idx, "lower_bound"]
        with model:
            cast(Reaction, model.reactions.get_by_id(glucose_rxn)).lower_bound = lb  # type: ignore[bad-arg-count]
            solution = model.optimize()
            df.loc[idx, "biomass"] = solution.fluxes["biomass_maintenance"]
    print(df)


def glycolysis_reactions_present(model: Model):
    """Validate major glycolysis reactions are present in the model."""
    found_all = True
    for rxn in GLYCOLYSIS_REACTIONS:
        if rxn not in model.reactions:
            print(f"Missing glycolysis reaction: {rxn}")
            found_all = False

    if found_all:
        print("All major glycolysis reactions are present in the model.")


def glucose_distribution(model: Model):
    """Validate flux through glycolytic reactions is non-zero."""

    def validate_flux():
        """Show that these reactions are doing 'something', even if it is zero-flux."""
        with model:
            model.objective = "biomass_maintenance"
            solution = model.optimize()

            for rxn_id in GLYCOLYSIS_REACTIONS:
                # this will fail if the reaction doen't exist
                # we expect these reactions to be in the model, and if they're not, that's a problem.
                if rxn_id not in solution.fluxes.index:
                    print(f"WARNING: Missing reaction '{rxn_id}'.")
                    continue

                flux = solution.fluxes[rxn_id]
                if abs(flux) <= 1e-6:
                    print(f"WARNING: Glycolysis reaction '{rxn_id}' has a (near-)zero  flux: {flux}")

    def set_sinks():
        metabolites = [
            "glc_D[c]",  # brought into model by exchange reaction
            "g6p[c]",  # produced by HEX1
            "f6p[c]",  # produced by PGI
            "g3p[c]",  # produced by FBA
            "13dpg[c]",  # produced by GAPD
            "2pg[c]",  # produced by PGM
            "pep[c]",  # produced by ENO
            "pyr[c]",  # produced by PYK
        ]
        for met in metabolites:
            with model:
                metabolite = cast(Metabolite, model.metabolites.get_by_id(met))
                sink = model.add_boundary(metabolite, type="sink")
                model.objective = sink.id
                solution = model.optimize()
                print(f"{met}: {solution.fluxes[sink.id]}")

    validate_flux()
    set_sinks()


def lactate_production(model: Model):
    """Validate that naive B cells produce lactate.

    Naive B cells produce lactate. Disrupt the following reaction:
    Define a sink reaction for lactate and determine how much metabolite is flowing through the pathways.
    """
    with model:
        model.add_boundary(
            metabolite=cast(Metabolite, model.metabolites.get_by_id("lac_D[c]")),
            type="sink",
            # reaction_id="sink_lac_D[c]"
        )  # type: ignore

        model.objective = "SK_lac_D[c]"
        initial_solution = model.optimize()

        with suppress(KeyError):
            rxn: Reaction = model.reactions.get_by_id("EX_lac_D[e]")  # type: ignore
            rxn.knock_out()

        solution_after = model.optimize()
        print(f"Flux through lac_d sink before knocking out EX_lac_D[e]: {initial_solution.fluxes['SK_lac_D[c]']}")
        print(f"Flux through lac_D sink after knocking out EX_lac_D[e]: {solution_after.fluxes['SK_lac_D[c]']}")


def citric_acid_cycle_activation(model: Model):
    """Validate that Naive B cells activate the citric acid cycle.

    In naive B cells, citrate, glutamate, glutamine, å-ketoglutarate, succinate, fumarate, malate are active
    Confirm with escher map
    """
    tca_cycle_rxns = [
        "CSm",
        "ACONTm",
        ("ICDHxm", "ICDHyrm"),
        "AKGDm",
        ("SUCOAS1m", "SUCOASm"),
        "SUCD1m",
        "FUMm",
        "MDHm",
    ]

    with model:
        model.objective = "biomass_maintenance"
        sol = model.optimize()

    fluxes = sol.fluxes
    flux_sum: dict[str | tuple[str, str], float] = defaultdict(float)

    for rxn in tca_cycle_rxns:
        if isinstance(rxn, str):
            if rxn not in fluxes.index:
                print(f"WARNING: Reaction '{rxn}' not found in the model.")
                continue
            flux_sum[rxn] += fluxes[rxn]
        elif isinstance(rxn, tuple):
            sums = 0
            for rxn_id in rxn:
                if rxn_id not in fluxes.index:
                    print(f"WARNING: Reaction '{rxn_id}' not found in the model.")
                    continue
                sums += fluxes[rxn_id]
            flux_sum[rxn] += sums

    df = pd.DataFrame.from_dict(flux_sum, orient="index", columns=["sum"])
    df.index.name = "reactions"
    print(df)


def ecar_response_with_glucose_addition(model: Model):
    """Validate that Naive B cells show no significant response in ECAR with glucose addition.

    ECAR: Extracellular Acidification Rate

    Untreated or anti-IgM treated B cells showed no significant increase in the ECAR in response to the glucose addition
    Define a sink reaction for lactate and maximize the exchange reaction of glucose. If the model agrees with the
    verification, the lactate present in the pathway should increase as well
    """
    with model:
        model.add_boundary(
            metabolite=cast(Metabolite, model.metabolites.get_by_id("lac_L[c]")),
            type="sink",
            reaction_id="sink_lac_L[c]",
            lb=0,
            ub=1000,
        )

        model.objective, model.objective.id = ("sink_lac_L[c]",) * 2
        wt_solution = model.optimize()

        orig_lower = cast(Reaction, model.reactions.get_by_id("EX_glc_D[e]")).lower_bound
        orig_upper = cast(Reaction, model.reactions.get_by_id("EX_glc_D[e]")).upper_bound
        cast(Reaction, model.reactions.get_by_id("EX_glc_D[e]")).lower_bound = 0  # type: ignore[bad-argument-count]
        cast(Reaction, model.reactions.get_by_id("EX_glc_D[e]")).upper_bound = 0  # type: ignore[bad-argument-count]
        modified_lower = cast(Reaction, model.reactions.get_by_id("EX_glc_D[e]")).lower_bound
        modified_upper = cast(Reaction, model.reactions.get_by_id("EX_glc_D[e]")).upper_bound

        lactate_solution = model.optimize()

        print("Before")
        print(f"Objective ID: {model.objective.id}")
        print(f"Glucose uptake lower bound: {orig_lower}")
        print(f"Glucose uptake upper bound: {orig_upper}")
        print(f"Objective: {wt_solution.fluxes['sink_lac_L[c]']}")
        print(f"Lactate flux: {wt_solution.fluxes['LDH_L']}")

        print("\nAfter")
        print(f"Objective ID: {model.objective.id}")
        print(f"Glucose uptake lower bound: {modified_lower}")
        print(f"Glucose uptake upper bound: {modified_upper}")
        print(f"Objective: {lactate_solution.fluxes['sink_lac_L[c]']}")
        print(f"Lactate flux: {lactate_solution.fluxes['LDH_L']}")


def glycolysis_and_oxidative_glucose_pathway_maintenance(model: Model):
    """Validate that Naive B cells require glycolytic and oxidative glucose pathway maintenance.

    B cell survival requires proper maintenance of glycolytic and oxidative glucose pathways.
    Use the moma method instead of flux balance analysis and delete one reaction at a time.
    The result reflects the change in the cell's growth rate. The biomass does not change.
    """
    with model :
        model_copy.objective = "biomass_maintenance"
        solution_before = model_copy.optimize()

        with (target := Path("/dev/null").open("w")), redirect_stderr(target), redirect_stdout(target):
            single_rxn_del = single_reaction_deletion(
                model_copy,
                [
                    "HEX1",  # Hexokinase
                    "PGI",  # Glucose-6-Phosphate isomerase
                    "PFK",  # Phosphofructokinase
                    "FBA",  # Fructose-Bisphosphate Aldolase
                    "GAPD",  # Glyceraldehyde 3-phosphate dehydrogenase
                    "PGK",  # Phosphoglycerate kinase
                    "PGM",  # Phosphoglyceromutase
                    "ENO",  # Enolase
                    "PYK",  # Pyruvate kinase
                ],
                method="moma",
                solution=solution_before,
                processes=9,
            )
        print(solution_before.fluxes["biomass_maintenance"])
        print(single_rxn_del)


def fatty_acid_fuel_source(model: Model):
    """Validate that naive B cells rely on fatty acids as their main fuel source."""
    with model :
        .objective = "biomass_maintenance"
        solution_before = .optimize()
        for rxn in .reactions:
            rxn: Reaction
            if rxn.subsystem == "Fatty acid oxidation":
                cast(Reaction, .reactions.get_by_id(rxn.id)).lower_bound = 0  # type: ignore[bad-argument-count]
                cast(Reaction, .reactions.get_by_id(rxn.id)).upper_bound = 0  # type: ignore[bad-argument-count]
        solution_after: cobra.Solution = .optimize()

    print(f"Objective flux before limiting fatty acid reactions: {solution_before.fluxes['biomass_maintenance']}")
    print(f"Objective flux after limiting fatty acid reactions: {solution_after.fluxes['biomass_maintenance']}")


if __name__ == "__main__":
    # fp = Path("/Users/joshl/Projects/COMO/main/data/results/naiveB/naiveB_imat_model_tpm.json")
    # fp = Path("/Users/joshl/Projects/COMO/main/data/results/naiveB/naiveB_imat_model_zfpkm.json")
    # fp = Path("/Users/joshl/Downloads/como_supp/Supplementary_Data_2_bbad387.xml")

    _root = Path("/Users/joshl/Projects/ImmunoMetabolism/results/model_build")
    # fp = _root / "A01/A01_b_naive/A01_b_naive_model_imat.json"
    # fp = _root / "A02/A02_b_naive/A02_b_naive_model_imat.json"
    # fp = _root / "A03/A03_b_naive/A03_b_naive_model_imat.json"
    # fp = _root / "B01/B01_b_naive/B01_b_naive_model_imat.json"
    # fp = _root / "B02/B02_b_naive/B02_b_naive_model_imat.json"
    # fp = _root / "B03/B03_b_naive/B03_b_naive_model_imat.json"
    # fp = _root / "C01/C01_b_naive/C01_b_naive_model_imat.json"
    # fp = _root / "C02/C02_b_naive/C02_b_naive_model_imat.json"
    # fp = _root / "C04/C04_b_naive/C04_b_naive_model_imat.json"
    # fp = _root / "D01/D01_b_naive/D01_b_naive_model_imat.json"
    # fp = _root / "D02/D02_b_naive/D02_b_naive_model_imat.json"
    # fp = _root / "D03/D03_b_naive/D03_b_naive_model_imat.json"
    # fp = _root / "E01/E01_b_naive/E01_b_naive_model_imat.json"
    # fp = _root / "E02/E02_b_naive/E02_b_naive_model_imat.json"
    fp = _root / "E03/E03_b_naive/E03_b_naive_model_imat.json"

    if fp.suffix == ".json":
        model_ = cobra.io.load_json_model(fp)
    elif fp.suffix == ".xml":
        model_ = cobra.io.read_sbml_model(fp)
    else:
        raise TypeError(f"Unknown extension '{fp.suffix}': {fp}")

    print("\n")
    model_ = constrain_model(model_)
    model_ = add_glycolysis_reactions(model_)
    model_ = add_exchanges(model_)
    print("\n")

    print(f"Reaction count: {len(model_.reactions)}")
    print(f"Gene count: {len(model_.genes)}")
    print(f"Metabolite count: {len(model_.metabolites)}")

    print("\n\nGlucose Restriction")
    glucose_restriction(model=model_)
    print("\n\nGlycolysis Reactions Present")
    glycolysis_reactions_present(model=model_)
    print("\n\nGlucose Distribution")
    glucose_distribution(model=model_)
    print("\n\nLactate Production")
    lactate_production(model=model_)
    print("\n\nCitric Acid Cycle Activation")
    citric_acid_cycle_activation(model=model_)
    print("\n\nECAR Response")
    ecar_response_with_glucose_addition(model=model_)
    print("\n\nFatty Acid Fuel Source")
    fatty_acid_fuel_source(model=model_)
    print("\n\nGlycolysis + Oxidative Glucose Pathway Maintenance")
    glycolysis_and_oxidative_glucose_pathway_maintenance(model=model_)
