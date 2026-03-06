from __future__ import annotations

import collections
import os
import re
import sys
from collections.abc import Sequence
from io import TextIOWrapper
from math import isfinite
from multiprocessing import Value
from pathlib import Path
from typing import Literal, Never, TextIO, cast, reveal_type

import cobamp.core.optimization
import cobra
import cobra.util.array
import numpy as np
import numpy.typing as npt
import pandas as pd
from loguru import logger
from troppo.methods.reconstruction.corda import CORDA, CORDAProperties
from troppo.methods.reconstruction.fastcore import FASTcore, FastcoreProperties
from troppo.methods.reconstruction.gimme import GIMME, GIMMEProperties
from troppo.methods.reconstruction.imat import IMAT, IMATProperties
from troppo.methods.reconstruction.tINIT import tINIT, tINITProperties

from como.data_types import (
    Algorithm,
    CobraCompartments,
    LogLevel,
    ModelBuildSettings,
    Solver,
    _BoundaryReactions,
)
from como.utils import set_up_logging, split_gene_expression_data


def _reaction_indices_to_ids(
    ref_model: cobra.Model, reaction_indices: Sequence[int] | npt.NDArray[np.integer]
) -> list[str]:
    rxns = list(ref_model.reactions)
    return [rxns[int(i)].id for i in reaction_indices]


def _correct_bracket(rule: str, name: str) -> str:
    """Correct GPR rules to format readable by.

    Args:
        rule: GPR rule string from a COBRA model
        name: Gene name string from a COBRA model

    Returns:
        A corrected GPR rule string
    """
    rule_match = re.search(r"or|and", rule)
    name_match = re.search(r"or|and", name)
    if rule_match is None or name_match is None:
        left_rule = rule
        left_name = name.strip()
        right_rule = ""
        right_name = ""
        operator = ""
    else:
        left_rule = rule[: rule_match.span()[0]]
        left_name = name[: name_match.span()[0]].strip()
        right_rule = rule[rule_match.span()[1] :]
        right_name = name[name_match.span()[1] :]
        operator = rule_match.group()

    new_right_rule: list[str] = []
    for char in list(left_rule):
        if char.isspace() or char.isdigit():
            new_right_rule.append(char)
        elif len(left_name) > 0 and char == left_name[0]:
            new_right_rule.append(char)
            left_name = left_name[1:]
    new_left_rule = "".join(new_right_rule).strip()
    final_right_rule = "" if rule_match is None else _correct_bracket(right_rule, right_name)
    return " ".join([new_left_rule, operator, final_right_rule]).strip()


def _gene_rule_logical(gpr_expression: str, level: int = 0) -> str:
    """Create an expression from GPR rule which can be evaluated as true or false.

    Args:
        gpr_expression: GPR rule string from a COBRA model
        level: Current recursion level (used for debugging)

    Returns:
        An evaluable string where "and" is replaced with "min" and "or" is replaced with "max"
    """
    try:
        loc_r = gpr_expression.index(")")
    except ValueError:
        if "and" in gpr_expression:
            gpr_expression = gpr_expression.replace(" and ", ", ")
            return "min{" + gpr_expression + "}"
        elif "or" in gpr_expression:
            gpr_expression = gpr_expression.replace(" or ", ", ")
            return "max{" + gpr_expression + "}"
        else:
            gpr_expression = gpr_expression.replace("[", "")
            return gpr_expression.replace("]", "")

    loc_l = gpr_expression[:loc_r].rindex("(")
    inner_string = gpr_expression[loc_l : loc_r + 1]
    inner_string = inner_string.replace("(", "[")
    inner_string = inner_string.replace(")", "]")
    if "and" in inner_string:
        inner_string = inner_string.replace("and", ",")
        inner_string = "min{" + inner_string + "}"
    elif "or" in inner_string:
        inner_string = inner_string.replace("or", ",")
        inner_string = "max{" + inner_string + "}"
    else:
        inner_string = inner_string.replace("[", "")
        inner_string = inner_string.replace("]", "")

    expression_out = f"{gpr_expression[:loc_l]}{inner_string}{gpr_expression[loc_r + 1 :]}"
    expression_out = _gene_rule_logical(expression_out, level + 1)

    return expression_out


def _feasibility_test(model_cobra: cobra.Model, step: str):
    # check number of unsolvable reactions for reference model under media assumptions
    # create flux consistant model (rmemoves some reactions)
    model_cobra_rm = cobra.flux_analysis.fastcc(model_cobra, flux_threshold=15, zero_cutoff=1e-7)
    incon_rxns = set(model_cobra.reactions.list_attr("id")) - set(model_cobra_rm.reactions.list_attr("id"))
    incon_rxns_cnt = len(incon_rxns)

    if step == "before_seeding":
        logger.warning(
            f"Under given boundary assumptions, there are {incon_rxns_cnt} infeasible reactions in"
            f" the reference model. These reactions will not be considered active in "
            f"context specific model construction. If any infeasible reactions are found to be "
            f"active according to expression data, or are found in the force reactions list, "
            f"they can be found found in 'InfeasibleRxns.csv'. It is normal for this value to be quite large; "
            f"however, if many of these reactions are active according to your expression data, "
            f"it is likely that you are missing some critical exchange (media) reactions."
        )
    elif step == "after_seeding":
        logger.warning(
            f"Under given boundary assumptions, with infeasible reactions from the general model not "
            f"considered there are {incon_rxns_cnt} new infeasible reactions in the context-specific model. "
            f"These reactions will be removed from the output model to ensure the model is solvable. "
            f"Note that this value should be very low compared to the reference model."
        )

    return incon_rxns, model_cobra_rm


def _build_with_gimme(
    reference_model: cobra.Model,
    expression_vector: Sequence[float] | npt.NDArray[np.floating],
    idx_objective: int,
    lower_bounds: npt.NDArray[np.floating],
    upper_bounds: npt.NDArray[np.floating],
    solver: str,
    threshold_percentile: int = 30,
):
    model = reference_model
    if threshold_percentile < 1:
        logger.critical(
            f"GIMME quantile calculation is less than 1 ({threshold_percentile}), but expected an integer of ~25! "
            f"This will likely result in a model with many more reactions than expected."
        )

    expression: npt.NDArray[np.floating] = np.asarray(expression_vector, dtype=float)
    nonzero = expression[expression > 0]
    min_expr_percentile = np.percentile(nonzero, threshold_percentile).astype(float)
    logger.debug(f"GIMME: minimum expression is {expression.min():.4f}")
    logger.debug(f"GIMME: minimum non-zero expression is {nonzero.min():.4f}")
    logger.debug(f"GIMME: maximum expression is {expression.max():.4f}")

    expression[(expression < min_expr_percentile) & (expression > 0)] = 0
    expression[expression >= min_expr_percentile] = 1

    suffix: str = {1: "st", 2: "nd", 3: "rd"}.get(threshold_percentile, "th")
    logger.debug(f"GIMME: {threshold_percentile}{suffix} percentile is {min_expr_percentile:.4f}")

    s_matrix: list[float] = list(np.asarray(cobra.util.array.create_stoichiometric_matrix(model=model), dtype=float))
    reaction_ids: list[str] = [r.id for r in model.reactions]
    metabolite_ids: list[str] = [m.id for m in model.metabolites]

    # `Becker and Palsson (2008). Context-specific metabolic networks are
    # consistent with experiments. PLoS Comput. Biol. 4, e1000082.`
    properties = GIMMEProperties(
        exp_vector=list(expression),
        objectives=[{idx_objective: 1}],
        preprocess=True,
        reaction_ids=reaction_ids,
        metabolite_ids=metabolite_ids,
        solver=solver.upper(),
    )
    algorithm = GIMME(s_matrix, list(lower_bounds), list(upper_bounds), properties)
    active_rxn_indices: list[int] = algorithm.run()
    to_remove_ids: list[str | cobra.Reaction] = [
        reaction_ids[i] for i in range(len(model.reactions)) if i not in active_rxn_indices
    ]
    model.remove_reactions(to_remove_ids, True)

    # psol = pfba(model)
    # reaction_ids = [r.id for r in context_cobra_model.reactions]
    # psol = context_cobra_model.optimize()
    # to_remove_ids = [reaction_ids[r] for r in np.where(abs(psol.fluxes) < 1e-8)[0]]
    # context_cobra_model.remove_reactions(to_remove_ids, True)

    return model


def _build_with_fastcore(
    reference_model: cobra.Model,
    lower_bounds: npt.NDArray[np.floating],
    upper_bounds: npt.NDArray[np.floating],
    exp_idx_list: Sequence[int],
    solver: str,
):
    # 'Vlassis, Pacheco, Sauter (2014). Fast reconstruction of compact
    # context-specific metabolic network models. PLoS Comput. Biol. 10,
    # e1003424.'
    model = reference_model
    logger.warning(
        "Fastcore requires a flux consistant model is used as refererence. "
        "To achieve this, fastcc is required, which is NOT reproducible."
    )
    s_matrix = cast(npt.NDArray[np.floating], cobra.util.create_stoichiometric_matrix(model=model))
    if lower_bounds.shape[0] != upper_bounds.shape[0] != s_matrix.shape[1]:
        raise ValueError(
            "Lower bounds, upper bounds, and stoichiometric matrix must have the same number of reactions."
        )
    logger.debug("Creating feasible model")
    _, reference_model = _feasibility_test(reference_model, "other")
    properties = FastcoreProperties(core=list(exp_idx_list), solver=solver)
    algorithm = FASTcore(s_matrix, np.asarray(lower_bounds), np.asarray(upper_bounds), properties)
    context_rxns = algorithm.fastcore()
    r_ids = [r.id for r in model.reactions]
    remove_rxns = [r_ids[int(i)] for i in range(s_matrix.shape[1]) if i not in context_rxns]
    model.remove_reactions(remove_rxns, True)

    return model


def _build_with_imat(
    reference_model: cobra.Model,
    lower_bounds: npt.NDArray[np.floating],
    upper_bounds: npt.NDArray[np.floating],
    expr_vector: npt.NDArray,
    low_expression_threshold: float,
    high_expression_threshold: float,
    force_reaction_indices: npt.NDArray[np.integer],
    solver: str,
    build_settings: ModelBuildSettings,
) -> cobra.Model:
    model = reference_model.copy()

    # `list()` is required because imat uses `if core:`; this returns an error on numpy arrays
    force_rxn_list = list(np.asarray(force_reaction_indices, dtype=int))
    properties: IMATProperties = IMATProperties(
        exp_vector=expr_vector,
        exp_thresholds=(
            low_expression_threshold,
            high_expression_threshold,
        ),
        tolerance=build_settings.min_reaction_flux,
        core=force_rxn_list,
        epsilon=build_settings.troppo_epsilon,
        solver=solver,
    )

    # Using cobra to create the stoichiometry matrix means we have less work to do
    s_matrix: npt.NDArray[np.floating] = cast(
        npt.NDArray[np.floating], cobra.util.array.create_stoichiometric_matrix(model=model)
    )
    algorithm: IMAT = IMAT(S=s_matrix, lb=lower_bounds, ub=upper_bounds, properties=properties)

    # Match Troppo's iMAT high/low bin calculations
    high_index = np.where(expr_vector >= high_expression_threshold)[0].astype(int)
    low_index = np.where((expr_vector >= 0) & (expr_vector < low_expression_threshold))[0].astype(int)
    if force_reaction_indices.size > 0:
        high_index = np.union1d(high_index, force_reaction_indices)
    lso, lsystem = algorithm.generate_imat_problem(
        S=algorithm.S,
        lb=algorithm.lb,
        ub=algorithm.ub,
        high_idx=high_index,
        low_idx=low_index,
        epsilon=build_settings.troppo_epsilon,
    )
    # reset default Troppo/Cobamp thresholds
    cfg = lsystem.get_configuration()
    cfg.verbosity = build_settings.solver_verbosity
    cfg.timeout = build_settings.solver_timeout
    cfg.tolerances.feasibility = build_settings.solver_feasibility
    cfg.tolerances.optimality = build_settings.solver_optimality
    cfg.tolerances.integrality = build_settings.solver_integrality
    grb_model = getattr(getattr(cfg, "problem", None), "problem", None)  # `getattr` to handle non-gurobi solvers
    if grb_model is not None:
        grb_model.Params.MIPGap = build_settings.gurobi_mipgap
        grb_model.Params.IntegralityFocus = build_settings.gurobi_integrality_focus
        grb_model.Params.NumericFocus = build_settings.gurobi_numeric_focus

    # Keep reactions from imat + force-included reactions
    solution: cobamp.core.optimization.Solution = lso.optimize()
    x_sol: npt.NDArray[np.floating] = solution.x().astype(float)

    # partition the solution into the following components:
    #   1) High-bin, forward-direction reactions
    #   2) High-bin, reverse-direction reactions
    #   3) Activated low-bin reactions
    n_rxns = len(reference_model.reactions)
    expected_len = n_rxns + (2 * high_index.size) + low_index.size
    if expected_len != x_sol.size:
        raise ValueError(
            f"Unexpected solution length: got {x_sol.size}, expected {expected_len}"
            f"(n_rxns={n_rxns}, n_high_index={high_index.size}, n_low_index={low_index.size}"
        )

    part: int = 0
    flux_vector: npt.NDArray[np.floating] = x_sol[part : part + n_rxns]
    part += n_rxns
    high_forward: npt.NDArray[np.floating] = x_sol[part : part + high_index.size]
    part += high_index.size
    high_reverse: npt.NDArray[np.floating] = x_sol[part : part + high_index.size]
    part += high_index.size + low_index.size

    # iMAT binarizes at ~0 and ~1 for off/on reactions,
    #   so we can use 0.5 as the midpoint for defining inclusion/exclusion
    high_forward_rxns: npt.NDArray[np.integer] = high_index[high_forward > 0.5]
    high_reverse_rxns: npt.NDArray[np.integer] = high_index[high_reverse > 0.5]
    high_on_rxns = np.unique(np.concatenate([high_forward_rxns, high_reverse_rxns])).astype(int)

    # Low bin reactions may be included if it has a flux greater than `min_reaction_flux`
    flux_keep = np.where(np.abs(solution.x())[:n_rxns] >= build_settings.min_reaction_flux)[0].astype(int)
    to_keep = np.unique(np.concatenate([flux_keep, high_on_rxns, force_reaction_indices]))

    # Collect all reaction IDs and their associated index (e.g., HEX1 is at index 123)
    all_rxn_ids: npt.NDArray[np.str_] = np.asarray([r.id for r in model.reactions], dtype=str)
    to_remove = np.ones(shape=len(model.reactions), dtype=bool)
    to_remove[to_keep] = False
    remove_rxn_ids: list[str | cobra.Reaction] = list(all_rxn_ids[to_remove])
    model.remove_reactions(reactions=remove_rxn_ids, remove_orphans=True)
    return model


def _build_with_tinit(
    reference_model: cobra.Model,
    lower_bounds: npt.NDArray[np.floating],
    upper_bounds: npt.NDArray[np.floating],
    expr_vector: npt.NDArray[np.floating],
    solver: str,
    idx_force: npt.NDArray[np.integer],
) -> None:
    raise NotImplementedError("tINIT is not yet implemented.")
    model = reference_model
    properties = tINITProperties(
        reactions_scores=list(expr_vector),
        solver=solver,
        essential_reactions=list(idx_force),
        production_weight=0.0,
        allow_excretion=False,
        no_reverse_loops=True,
    )
    s_matrix: npt.NDArray[np.floating] = np.asarray(
        cobra.util.array.create_stoichiometric_matrix(model=model), dtype=float
    )
    algorithm = tINIT(s_matrix, np.asarray(lower_bounds), np.asarray(upper_bounds), properties)
    active_rxn_indices: npt.NDArray[np.integer] = algorithm.run()


def _build_with_corda(
    reference_model: cobra.Model,
    neg_expression_threshold: float,
    high_expression_threshold: float,
    lower_bounds: npt.NDArray[np.floating],
    upper_bounds: npt.NDArray[np.floating],
    expression_vector: Sequence[float] | npt.NDArray[np.floating],
):
    """Reconstruct a model using CORDA.

    :param neg_expression_threshold: Reactions expressed below this value will be placed in "negative" expression bin
    :param high_expression_threshold: Reactions expressed above this value will be placed in the "high" expression bin
    """
    raise NotImplementedError("CORDA is not yet implemented")
    model = reference_model
    properties = CORDAProperties(
        high_conf_rx=[],
        medium_conf_rx=[],
        neg_conf_rx=[],
        pr_to_np=2,
        constraint=1,
        constrainby="val",
        om=1e4,
        ntimes=5,
        nl=1e-2,
        solver="GUROBI",
        threads=5,
    )
    s_matrix: npt.NDArray[np.floating] = np.asarray(
        cobra.util.array.create_stoichiometric_matrix(model=model), dtype=float
    )
    algorithm = CORDA(S=s_matrix, lb=np.asarray(lower_bounds), ub=np.asarray(upper_bounds), properties=properties)
    active_rxn_indices: npt.NDArray[np.integer] = algorithm.run()


def _map_expression_to_reaction(
    reference_model: cobra.Model,
    gene_expression_file: Path,
    recon_algorithm: Algorithm,
    taxon: int | str,
    low_percentile: int,
    high_percentile: int,
) -> tuple[collections.OrderedDict[str, float], float, float, float]:
    """Map gene ids to a reaction based on GPR (gene to protein to reaction) association rules.

    These rules should be defined in the general genome-scale metabolic model

    Args:
        reference_model: A COBRA model object representing the general genome-scale metabolic model.
        gene_expression_file: Path to a gene expression file (.csv, .tsv, .xlsx, or .xls)
        recon_algorithm: Algorithm to use for reconstruction (GIMME, FASTCORE, iMAT, or tINIT)
        taxon: Taxon ID or Taxon object for gene ID conversion.
        low_percentile: Low percentile for threshold calculation
        high_percentile: High percentile for threshold calculation

    Returns:
        [0]: An ordered dictionary mapping reaction IDs to their corresponding expression values.
        [1]: The abs(minimum value) from the expression vector
        [2]: The adjusted `low_expr` value
        [3]: The adjusted `high_expr` value

    Raises:
        ValueError: If neither 'entrez_gene_id' nor 'ensembl_gene_id' columns are found in the gene expression file.
    """
    expression_data = pd.read_csv(gene_expression_file)
    gene_activity = split_gene_expression_data(
        expression_data, identifier_column="entrez_gene_id", recon_algorithm=recon_algorithm
    )

    # adjust expression to have a minimum value of 0; if the minimum value is already 0, this does nothing
    gene_arr = gene_activity["active"].to_numpy(dtype=float, copy=False)
    min_val = float(gene_arr[np.isfinite(gene_arr)].min())
    if min_val < 0:
        gene_arr[np.isfinite(gene_arr)] += abs(min_val)

    # We are computing percentiles on non-zero data because we need to make sure that
    #   zero values are definitely placed in the low bin
    low_expr, high_expr = np.nanpercentile(gene_arr[gene_arr > 0], [low_percentile, high_percentile])
    gene_activity["active"] = gene_arr
    gene_activity.index = gene_activity.index.astype(str)  # maintain strings for reference model gene_ids

    # fmt: off
    # Define a default expression values for the following conditions:
    # 1. `no_expression`: If a reaction has a gene reaction rule but is not found in the data, then
    #       we can be positive that the reaction should be marked as inactive
    #   For all algorithms, placing no-expression data at an expression value of 0.0 will
    #       mark it as inactive, because we have adjusted the data to be >= 0 above
    no_expression = 0
    # 2. `missing_gpr`: If a reaction has no gene reaction rule, we don't know if it should be marked as (in)active
    #   Thus, place it in a no-cost location where possible, otherwise mark it as inactive
    #   IMAT/TINIT/GIMME: Place in mid-bin, which is a "no-cost weight" value
    #   FASTCORE: Mark as inactive, which is the best we can do with FASTCORE
    missing_gpr = (
        # -1 if recon_algorithm in {Algorithm.IMAT, Algorithm.TINIT, Algorithm.GIMME}
        (low_expr + high_expr) / 2 if recon_algorithm in {Algorithm.IMAT, Algorithm.TINIT, Algorithm.GIMME}
        else 0 if recon_algorithm == Algorithm.FASTCORE
        else 1
    )
    # fmt: on

    reaction_expression: collections.OrderedDict[str, float] = collections.OrderedDict()
    activity_map = cast(dict[str, float], gene_activity.loc[:, "active"].to_dict())
    error_count = 0
    _base_warn = "could not generate evaluable gene rule."
    for rxn in reference_model.reactions:
        rxn: cobra.Reaction
        gene_reaction_rule = rxn.gene_reaction_rule

        if not gene_reaction_rule:
            reaction_expression[rxn.id] = missing_gpr
            continue

        gene_ids: set[str] = set(re.findall(r"\d+", gene_reaction_rule))
        for gene_id in gene_ids:
            activity = activity_map.get(gene_id, no_expression)
            if np.isposinf(activity):  # For positive infinity values, place in high bin
                activity = high_expr + 1
            elif not isfinite(activity):  # For non-finite values, place in low/non-active bin
                activity = no_expression

            # Matches the standalone gene ID only when it's not part of a larger number.
            # It ensures the ID is not immediately preceded by a digit or decimal point and not immediately followed by
            #   a digit, so if the gene id is '48', decimals like 1.48 or numbers like 1648 are left unchanged.
            gene_reaction_rule = re.sub(
                pattern=rf"(?<![\d.]){gene_id}(?!\d)",
                repl=str(activity),
                string=str(gene_reaction_rule),
            )

        evaluable_gene_rule = None
        try:
            # We are using eval here because ast.literal_eval is unable to process an evaluable such as `max(-4, 0, 1)`
            # This isn't ideal, but ultimately the only other option is writing and maintaining a custom parsing engine
            #   but this is too much work
            evaluable_gene_rule = _gene_rule_logical(gene_reaction_rule).replace("{", "(").replace("}", ")")
            reaction_expression[rxn.id] = eval(evaluable_gene_rule)  # noqa: S307

        except ValueError:
            logger.warning(f"ValueError for Reaction '{rxn.id}': {evaluable_gene_rule or _base_warn}")
            error_count += 1
        except NameError:
            logger.warning(f"NameError for Reaction '{rxn.id}': {evaluable_gene_rule or _base_warn}")
            error_count += 1
        except SyntaxError as e:
            raise SyntaxError(f"SyntaxError on Reaction '{rxn.id}': {evaluable_gene_rule or _base_warn}") from e
        except Exception as e:
            error_type = type(e).__name__
            raise RuntimeError(
                f"Unexpected error {error_type} on Reaction '{rxn.id}': {evaluable_gene_rule or _base_warn}"
            ) from e

    logger.debug(f"Mapped gene expression to reactions, found {error_count} error(s).")
    return reaction_expression, float(min_val), float(low_expr), float(high_expr)


def _read_reference_model(filepath: Path) -> cobra.Model:
    if filepath.suffix == ".json":
        return cobra.io.load_json_model(filepath)
    if filepath.suffix == ".mat":
        return cobra.io.load_matlab_model(filepath)
    if filepath.suffix in {".xml", ".sbml"}:
        return cobra.io.read_sbml_model(filepath)
    raise ValueError(f"Reference model format must be .xml, .mat, or .json; found '{filepath.suffix}'")


def _build_model(
    reference_model: cobra.Model,
    gene_expression_file: Path,
    recon_algorithm: Algorithm,
    objective: str,
    boundary_reactions: list[str],
    exclude_reactions: list[str],
    force_reactions: list[str],
    lower_bounds: list[float],
    upper_bounds: list[float],
    solver: str,
    low_percentile: int,
    high_percentile: int,
    output_flux_result_filepath: Path,
    taxon: int | str,
    objective_direction: Literal["min", "max"],
    build_settings: ModelBuildSettings,
    *,
    force_boundary_rxn_inclusion: bool,
) -> cobra.Model:
    """Seed a context specific reference_model.

    Core reactions are determined from GPR associations with gene expression logicals.
    Core reactions that do not necessarily meet GPR association requirements can be forced if in the force reaction
    file. Metabolite exchange (media), sinks, and demands are determined from exchanges file. Reactions can also be
    force excluded even if they meet GPR association requirements using the force exclude file.

    Args:
        reference_model: The reference model used in reconstruction
        gene_expression_file: Path to a gene expression file (.csv, .tsv, .xlsx, or .xls)
        recon_algorithm: Algorithm to use for reconstruction (GIMME, FASTCORE, iMAT, or tINIT)
        objective: Objective reaction ID in the general model
        boundary_reactions: List of boundary reactions to set in the model
        exclude_reactions: List of reactions to exclude from the model
        force_reactions: List of reactions to force include in the model
        lower_bounds: List of lower bounds corresponding to boundary reactions
        upper_bounds: List of upper bounds corresponding to boundary reactions
        solver: Solver to use (e.g., 'glpk', 'cplex', 'gurobi')
        output_flux_result_filepath: Path to save flux results (for iMAT only)
        taxon: Taxon ID or Taxon object for gene ID conversion.
        force_boundary_rxn_inclusion: If True, ensure that all boundary reactions are included in the final model.

    Returns:
        A _BuildResults object containing:
         the context-specific model
         list of expression indices used
        A DataFrame of infeasible reactions.
    """
    if objective not in force_reactions:
        force_reactions.append(objective)  # forcibly include the objective in the reconstructed model

    # Set reference model boundaries
    for rxn, lb, ub in zip(boundary_reactions, lower_bounds, upper_bounds, strict=True):
        cast(cobra.Reaction, reference_model.reactions.get_by_id(rxn)).lower_bound = lb  # type: ignore[bad-argument-count]
        cast(cobra.Reaction, reference_model.reactions.get_by_id(rxn)).upper_bound = ub  # type: ignore[bad-argument-count]
    reference_model.solver = solver.lower()  # type: ignore[bad-argument-count]
    reference_model.objective_direction = objective_direction  # type: ignore[bad-argument-count]

    # Collect the lower/upper bounds for reactions in the reference model to provide it to reconstruction methods
    ref_lb: npt.NDArray[np.floating] = np.full(shape=(len(reference_model.reactions)), fill_value=np.nan, dtype=float)
    ref_ub: npt.NDArray[np.floating] = np.full(shape=(len(reference_model.reactions)), fill_value=np.nan, dtype=float)
    reaction_ids: list[str] = []
    for i, rxn in enumerate(reference_model.reactions):
        rxn: cobra.Reaction
        ref_lb[i] = float(rxn.lower_bound)
        ref_ub[i] = float(rxn.upper_bound)
        reaction_ids.append(rxn.id)
    if ref_lb.shape[0] != ref_ub.shape[0] != len(reaction_ids):
        raise ValueError(
            "Lower bounds, upper bounds, and reaction IDs must have the same length.\n"
            f"Number of reactions: {len(reaction_ids)}\n"
            f"Number of upper bounds: {ref_ub.shape[0]}\n"
            f"Number of lower bounds: {ref_lb.shape[0]}"
        )
    if np.isnan(ref_lb).any():
        raise ValueError("Lower bounds contains unfilled values!")
    if np.isnan(ref_ub).any():
        raise ValueError("Upper bounds contains unfilled values!")

    reaction_expression, min_expr_val, low_expr, high_expr = _map_expression_to_reaction(
        reference_model,
        gene_expression_file,
        recon_algorithm,
        taxon=taxon,
        low_percentile=low_percentile,
        high_percentile=high_percentile,
    )
    expression_vector: npt.NDArray[np.floating] = np.asarray(list(reaction_expression.values()), dtype=float)
    for rxn_id in force_reactions:
        if rxn_id not in reaction_ids:
            logger.warning(
                f"The force reaction '{rxn_id}' was not found in the reference model. "
                f"Check BiGG, or the relevant database for your reference model, for synonyms."
            )

    # Set force-include reactions to have an expression higher than the upper threshold
    # Set force-exclude reactions to have an expression lower than the lower threshold
    rxn_ids = list(reaction_expression.keys())
    force_mask = np.isin(rxn_ids, force_reactions)
    exclude_mask = np.isin(rxn_ids, exclude_reactions)
    expression_vector[force_mask] = high_expr + 1 if recon_algorithm in {Algorithm.TINIT, Algorithm.IMAT} else 1
    expression_vector[exclude_mask] = (
        max(0.0, low_expr - 1) if recon_algorithm in {Algorithm.TINIT, Algorithm.IMAT} else 0
    )

    if force_boundary_rxn_inclusion:
        all_forced: set[str] = {*force_reactions, *boundary_reactions}
        force_reaction_indices: npt.NDArray[np.uint16] = np.array(
            [reaction_ids.index(rxn) for rxn in all_forced if rxn in reaction_ids], dtype=np.uint16
        )
    else:
        force_reaction_indices: npt.NDArray[np.uint16] = np.array(
            [reaction_ids.index(rxn) for rxn in force_reactions if rxn in reaction_ids], dtype=np.uint16
        )

    expression_vector_indices = [i for (i, val) in enumerate(expression_vector) if val > 0]
    expression_threshold = (low_thresh, high_thresh)

    if recon_algorithm == Algorithm.IMAT:
        context_model_cobra: cobra.Model = _build_with_imat(
            reference_model=reference_model,
            lower_bounds=ref_lb,
            upper_bounds=ref_ub,
            expr_vector=expression_vector,
            low_expression_threshold=updated_low_thresh,
            high_expression_threshold=updated_high_thresh,
            force_reaction_indices=force_reaction_indices,
            solver=solver,
            build_settings=build_settings,
        )
    elif recon_algorithm == Algorithm.GIMME:
        expressed_rxn_ids: list[str] = list(reaction_expression.keys())
        metabolite_ids: set[str] = set()
        for rxn_id in expressed_rxn_ids:
            cobra_rxn: cobra.Reaction = cast(cobra.Reaction, reference_model.reactions.get_by_id(rxn_id))
            metabolite_ids.update([m.id for m in cobra_rxn.metabolites])

        context_model_cobra: cobra.Model = _build_with_gimme(
            reference_model=reference_model,
            expression_vector=expression_vector,
            idx_objective=objective_index,
            lower_bounds=ref_lb,
            upper_bounds=ref_ub,
            solver=solver,
        )
    elif recon_algorithm == Algorithm.FASTCORE:
        context_model_cobra: cobra.Model = _build_with_fastcore(
            reference_model=reference_model,
            lower_bounds=ref_lb,
            upper_bounds=ref_ub,
            exp_idx_list=expression_vector_indices,
            solver=solver,
        )
    elif recon_algorithm == Algorithm.TINIT:
        context_model_cobra: cobra.Model = _build_with_tinit(
            reference_model=reference_model,
            lower_bounds=ref_lb,
            upper_bounds=ref_ub,
            expr_vector=expression_vector,
            solver=solver,
            idx_force=force_reaction_indices,
        )
    else:
        raise ValueError(
            f"Reconstruction algorithm must be {Algorithm.GIMME.value}, "
            f"{Algorithm.FASTCORE.value}, {Algorithm.IMAT.value}, or {Algorithm.TINIT.value}. "
            f"Got: {recon_algorithm.value}"
        )

    inconsistent_and_infeasible_reactions: pd.DataFrame = pd.concat(
        [
            pd.DataFrame({"infeasible_reactions": inconsistent_reactions}),
            pd.DataFrame({"expressed_infeasible_reactions": infeasible_expression_reactions}),
            pd.DataFrame({"infeasible_force_reactions": infeasible_force_reactions}),
            pd.DataFrame({"infeasible_context_reactions": []}),  # Included to maintain legacy support
        ],
        ignore_index=True,
        axis=0,
    )

    return _BuildResults(
        model=context_model_cobra,
        expression_index_list=expression_vector_indices,
        infeasible_reactions=inconsistent_and_infeasible_reactions,
    )



def _create_df(path: Path, *, lowercase_col_names: bool = False) -> pd.DataFrame:
    if path.suffix not in {".csv", ".tsv"}:
        raise ValueError(f"File must be a .csv or .tsv file, got '{path.suffix}'")
    df: pd.DataFrame = await _read_file(path=path, header=0, sep="," if path.suffix == ".csv" else "\t", h5ad_as_df=True)

    if not isinstance(df, pd.DataFrame):
        _log_and_raise_error(
            f"Expected a pandas.DataFrame, got {type(df)}",
            error=TypeError,
            level=LogLevel.ERROR,
        )

    if lowercase_col_names:
        df.columns = [c.lower() for c in df.columns]
    return df


def _collect_boundary_reactions(
    path: Path,
    reference_model: cobra.Model,
    *,
    close_unlisted_exchanges: bool,
) -> _BoundaryReactions:
    df: pd.DataFrame = _create_df(path, lowercase_col_names=True)

    for column in df.columns:
        if column not in [
            "reaction",
            "abbreviation",
            "compartment",
            "lower bound",
            "upper bound",
        ]:
            raise ValueError(
                f"Boundary reactions file must have columns named 'reaction', 'abbreviation', 'compartment', "
                f"'lower bound', and 'upper bound'. Found: '{column}'"
            )

    for compartment in df["compartment"].unique():
        as_shorthand = CobraCompartments.get_shorthand(compartment) if len(compartment) > 2 else compartment
        if as_shorthand not in reference_model.compartments:
            raise ValueError(f"Unknown compartment: '{compartment}'")

    boundary_map = {"exchange": "EX", "demand": "DM", "sink": "SK"}
    for exchange in df["reaction"].str.lower().unique():
        if boundary_map.get(exchange) is None:
            raise ValueError(
                f"Unable to determine reaction type for: '{exchange}'\n"
                f"Validate the boundary reactions file has a 'Reaction' value of: 'Exchange', 'Demand', or 'Sink'."
            )

    reactions: list[str] = []
    lower_bounds: list[float] = []
    upper_bounds: list[float] = []
    for ref_rxn in reference_model.boundary:
        ref_rxn: cobra.Reaction

        # Remove prefixes from reactions: EX_glc_D[e] -> glc_D[e]; DM_h2o[c] -> h2o[c]; sink_o2[m] -> o2[m]
        # Then remove suffixes from reactions: glc_D[e] -> glc_D; h2o[c] -> h2o; o2[m] -> o2
        base_rxn_id: str = ref_rxn.id.split("_", maxsplit=1)[1]
        base_rxn_id = base_rxn_id.split("[", maxsplit=1)[0]

        # If we are excluding unlisted exchange reactions from the model,
        #   then we need to set their lower/upper bounds to 0
        # Otherwise, keep the unlisted exchanges as their default
        if base_rxn_id not in df["abbreviation"].values:
            reactions.append(ref_rxn.id)
            lower_bounds.append(0 if close_unlisted_exchanges else ref_rxn.lower_bound)
            upper_bounds.append(0 if close_unlisted_exchanges else ref_rxn.upper_bound)
            continue

        rxn_index = df[df["abbreviation"] == base_rxn_id].index.item()
        shorthand_compartment = CobraCompartments.get_shorthand(str(df.at[rxn_index, "compartment"]))
        exch_type: str | None = boundary_map.get(str(df.at[rxn_index, "reaction"]).lower())
        lower_bound, upper_bound = df.loc[rxn_index, ["lower bound", "upper bound"]].astype(float).values.ravel()

        # The boundary reactions are built on the base reaction id;
        #   if the reference model boundary has a matching base ID, but a non-matching
        #   boundary type (exchange, sink, demand), this will check to make sure we are not
        #   mistakingly setting boundaries that should not exist
        formatted_boundary_reaction = f"{exch_type}_{base_rxn_id}[{shorthand_compartment}]"
        if ref_rxn.id == formatted_boundary_reaction:
            reactions.append(formatted_boundary_reaction)
            lower_bounds.append(lower_bound)
            upper_bounds.append(upper_bound)
        else:
            reactions.append(ref_rxn.id)
            lower_bounds.append(0 if close_unlisted_exchanges else ref_rxn.lower_bound)
            upper_bounds.append(0 if close_unlisted_exchanges else ref_rxn.upper_bound)

    return _BoundaryReactions(
        reactions=reactions,
        lower_bounds=lower_bounds,
        upper_bounds=upper_bounds,
    )


def _write_model_to_disk(
    context_name: str,
    model: cobra.Model,
    output_filepaths: list[Path],
    mat_suffix: set[str],
    json_suffix: set[str],
    xml_suffix: set[str],
) -> None:
    for path in output_filepaths:
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.suffix in mat_suffix:
            cobra.io.save_matlab_model(model=model, file_name=path)
        elif path.suffix in json_suffix:
            cobra.io.save_json_model(model=model, filename=path, pretty=True)
        elif path.suffix in xml_suffix:
            cobra.io.write_sbml_model(cobra_model=model, filename=path)
        else:
            raise ValueError(
                f"Invalid output model filetype. Should be one of .xml, .sbml, .mat, or .json. Got '{path.suffix}'"
            )
        logger.success(f"Model for context '{context_name}' saved to to: '{path}'")


def create_context_specific_model(  # noqa: C901
    context_name: str,
    taxon: int | str | Taxon,
    reference_model: Path,
    active_genes_filepath: Path,
    output_infeasible_reactions_filepath: Path,
    output_flux_result_filepath: Path,
    output_model_filepaths: Path | list[Path],
    output_fastcore_expression_index_filepath: Path | None = None,
    objective: str = "biomass_reaction",
    boundary_rxns_filepath: str | Path | None = None,
    exclude_rxns_filepath: str | Path | None = None,
    force_rxns_filepath: str | Path | None = None,
    algorithm: Algorithm = Algorithm.GIMME,
    low_threshold: float = -5,
    high_threshold: float = -3,
    solver: Solver = Solver.GLPK,
    log_level: LogLevel = LogLevel.INFO,
    log_location: str | TextIO | TextIOWrapper = sys.stderr,
    *,
    force_boundary_rxn_inclusion: bool = False,
    close_unlisted_exchanges: bool = False,
) -> cobra.Model:
    """Create a context-specific model using the provided data.

    Args:
        context_name: Name of the context-specific model.
        taxon: NCBI taxonomy ID or name for the organism of interest.
        reference_model: Path to the general genome-scale metabolic model file (.xml, .mat, or .json).
        active_genes_filepath: Path to the gene expression data file (csv, tsv, or Excel).
        output_infeasible_reactions_filepath: Path to save infeasible reactions (csv).
        output_flux_result_filepath: Path to save flux results (csv).
        output_model_filepaths: Path or list of paths to save the context-specific model (.xml, .mat, or .json).
        output_fastcore_expression_index_filepath: Path to save Fastcore expression indices (txt). Required if using Fastcore.
        objective: Objective function reaction ID.
        boundary_rxns_filepath: Optional path to boundary reactions file (csv, tsv, or Excel).
        exclude_rxns_filepath: Optional path to reactions to exclude file (csv, tsv, or Excel).
        force_rxns_filepath: Optional path to reactions to force include file (csv, tsv, or Excel).
        algorithm: Algorithm to use for reconstruction. One of Algorithm.GIMME, Algorithm.FASTCORE, Algorithm.IMAT, Algorithm.TINIT.
        low_threshold: Low expression threshold for algorithms that require it.
        high_threshold: High expression threshold for algorithms that require it.
        solver: Solver to use. One of Solver.GLPK, Solver.CPLEX, Solver.GUROBI
        log_level: Logging level. One of LogLevel.DEBUG, LogLevel.INFO, LogLevel.WARNING, LogLevel.ERROR, LogLevel.CRITICAL
        log_location: Location for log output. Can be a file path or sys.stderr/sys.stdout.
        force_boundary_rxn_inclusion: If True, ensure that all provided boundary reactions are included in the final model.

    Raises:
        ImportError: If Gurobi solver is selected but gurobipy is not installed.
    """
    set_up_logging(level=log_level, location=log_location)
    if low_percentile is None:
        raise ValueError("low_percentile must be provided")
    if high_percentile is None:
        raise ValueError("high_percentile must be provided")
    # TODO: set up zfpkm threshold defaults

    boundary_rxns_filepath: Path | None = Path(boundary_rxns_filepath) if boundary_rxns_filepath else None
    output_model_filepaths = [output_model_filepaths] if isinstance(output_model_filepaths, Path) else output_model_filepaths

    if not reference_model.exists():
        raise FileNotFoundError(f"Reference model not found at {reference_model}")
    if not active_genes_filepath.exists():
        raise FileNotFoundError(f"Active genes file not found at {active_genes_filepath}")
    if algorithm == Algorithm.FASTCORE and not output_fastcore_expression_index_filepath:
        raise ValueError("The fastcore expression index output filepath must be provided")
    if boundary_rxns_filepath and not boundary_rxns_filepath.exists():
        raise FileNotFoundError(f"Boundary reactions file not found at {boundary_rxns_filepath}")

    if algorithm not in Algorithm:
        raise ValueError(f"Algorithm {algorithm} not supported. Use one of {', '.join(a.value for a in Algorithm)}")

    if solver not in Solver:
        raise ValueError(f"Solver '{solver}' not supported. Use one of {', '.join(s.value for s in Solver)}")

    mat_suffix, json_suffix, xml_suffix = {".mat"}, {".json"}, {".sbml", ".xml"}
    if any(path.suffix not in {*mat_suffix, *json_suffix, *xml_suffix} for path in output_model_filepaths):
        invalid_suffix: str = "\n".join(
            path.as_posix()
            for path in output_model_filepaths
            if path.suffix not in {*mat_suffix, *json_suffix, *xml_suffix}
        )
        raise ValueError(f"Invalid output filetype. Should be 'xml', 'sbml', 'mat', or 'json'. Got:\n{invalid_suffix}'")

    reference_model = _read_reference_model(reference_model_filepath)
    boundary_reactions = (
        _collect_boundary_reactions(
            boundary_rxns_filepath,
            reference_model=reference_model,
            close_unlisted_exchanges=close_unlisted_exchanges,
        )

    boundary_reactions = None
    if boundary_rxns_filepath:
        boundary_reactions = await _collect_boundary_reactions(boundary_rxns_filepath)

    exclude_rxns: list[str] = []
    if exclude_rxns_filepath:
        exclude_rxns_filepath: Path = Path(exclude_rxns_filepath)
        df = await _create_df(exclude_rxns_filepath)
        if "abbreviation" not in df.columns:
            raise ValueError("The exclude reactions file should have a single column with a header named Abbreviation")
        exclude_rxns = df["abbreviation"].tolist()

    force_rxns: list[str] = []
    if force_rxns_filepath:
        force_rxns_filepath: Path = Path(force_rxns_filepath)
        df = await _create_df(force_rxns_filepath, lowercase_col_names=True)
        if "abbreviation" not in df.columns:
            raise ValueError("The force reactions file should have a single column with a header named Abbreviation")
        force_rxns = df["abbreviation"].tolist()

    # Test that gurobi is using a valid license file
    if solver == Solver.GUROBI:
        from importlib.util import find_spec

        gurobi_present = find_spec("gurobipy")
        if not gurobi_present:
            raise ImportError(
                "The gurobi solver requires the gurobipy package to be installed. "
                "Please install gurobipy and try again. "
                "This can be done by installing the 'gurobi' optional dependency."
            )

        if not Path(f"{os.environ['HOME']}/gurobi.lic").exists():
            logger.critical(
                "Gurobi solver requested, but license information cannot be found. "
                "COMO will continue, but it is HIGHLY unlikely the resulting model will be valid."
            )

    logger.info(f"context='{context_name}', reconstruction method='{algorithm.value}', solver='{solver.value}'")
    build_results: _BuildResults = await _build_model(
        general_model_file=reference_model,
        gene_expression_file=active_genes_filepath,
        recon_algorithm=algorithm,
        objective=objective,
        boundary_reactions=boundary_reactions.reactions if boundary_reactions else [],
        exclude_reactions=exclude_rxns,
        force_reactions=force_rxns,
        lower_bounds=boundary_reactions.lower_bounds if boundary_reactions else [],
        upper_bounds=boundary_reactions.upper_bounds if boundary_reactions else [],
        solver=solver.value.lower(),
        low_thresh=low_threshold,
        high_thresh=high_threshold,
        output_flux_result_filepath=output_flux_result_filepath,
        force_boundary_rxn_inclusion=force_boundary_rxn_inclusion,
        taxon=taxon,
    )

    build_results.infeasible_reactions.dropna(inplace=True)
    build_results.infeasible_reactions.to_csv(output_infeasible_reactions_filepath, index=False)

    if algorithm == Algorithm.FASTCORE:
        fastcore_df = pd.DataFrame(build_results.expression_index_list)
        fastcore_df.dropna(inplace=True)
        fastcore_df.to_csv(output_fastcore_expression_index_filepath, index=False)

    await _write_model_to_disk(
        context_name=context_name,
        model=build_results.model,
        output_filepaths=output_model_filepaths,
        mat_suffix=mat_suffix,
        json_suffix=json_suffix,
        xml_suffix=xml_suffix,
    )
    logger.info(
        f"context={context_name}, reactions={len(build_results.model.reactions)}, genes={len(build_results.model.genes)}, metabolites={len(build_results.model.metabolites)}"
    )
