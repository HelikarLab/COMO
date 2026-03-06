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
    model_reconstruction = reference_model.copy()
    s_matrix: list[float] = list(cobra.util.array.create_stoichiometric_matrix(model=model_reconstruction))
    # `Becker and Palsson (2008). Context-specific metabolic networks are
    # consistent with experiments. PLoS Comput. Biol. 4, e1000082.`
    properties = GIMMEProperties(
        exp_vector=expr_vector,  # np.array(gimme_data['0']),
        obj_frac=0.9,
        objectives=[{idx_objective: 1}],
        preprocess=True,
        flux_threshold=0.9,
    )
    algorithm = GIMME(s_matrix, list(lower_bounds), list(upper_bounds), properties)
    gene_activity = algorithm.run()
    reaction_ids = [r.id for r in model_reconstruction.reactions]
    to_remove_ids = [reaction_ids[r] for r in np.where(gene_activity == 0)[0]]

    model_reconstruction.remove_reactions(to_remove_ids, True)
    psol = pfba(model_reconstruction)  # noqa: F841
    # reaction_ids = [r.id for r in context_cobra_model.reactions]
    # psol = context_cobra_model.optimize()
    # to_remove_ids = [reaction_ids[r] for r in np.where(abs(psol.fluxes) < 1e-8)[0]]
    # context_cobra_model.remove_reactions(to_remove_ids, True)

    return model_reconstruction


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
    _, cobra_model = _feasibility_test(cobra_model, "other")
    properties = FastcoreProperties(core=exp_idx_list, solver=solver)
    algorithm = FASTcore(s_matrix, lower_bounds, upper_bounds, properties)
    context_rxns = algorithm.fastcore()
    context_cobra_model = cobra_model.copy()
    r_ids = [r.id for r in context_cobra_model.reactions]
    remove_rxns = [r_ids[int(i)] for i in range(s_matrix.shape[1]) if i not in context_rxns]
    context_cobra_model.remove_reactions(remove_rxns, True)

    return context_cobra_model


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
    properties: IMATProperties = IMATProperties(
        exp_vector=expr_vector,
        exp_thresholds=expr_thresh,
        # core=np.array(force_gene_indices).tolist(),
        core=list(force_reaction_indices),  # `list()` is required because imat uses `if core:`; this returns an error on numpy arrays
        epsilon=0.01,
        solver=solver.upper(),
    )

    # Creating a copy of the model ensures we don't make any in-place modifications by accident
    # Using cobra to create the stoichiometry matrix means we have less work to do
    force_reaction_indices = np.array(force_reaction_indices, dtype=np.uint16)
    model_reconstruction: cobra.Model = reference_model.copy()
    s_matrix: npt.NDArray[float] = cobra.util.array.create_stoichiometric_matrix(model=model_reconstruction)
    algorithm: IMAT = IMAT(S=s_matrix, lb=np.array(lower_bounds), ub=np.array(upper_bounds), properties=properties)
    rxns_from_imat: npt.NDArray[np.uint16] = algorithm.run().astype(np.uint16)

    # Collect all reaction IDs and their associated index (e.g., HEX1 is at index 123)
    all_rxn_ids: npt.NDArray[str] = np.array([r.id for r in model_reconstruction.reactions], dtype=object)
    all_rxn_indices: npt.NDArray[np.uint16] = np.array(range(len(model_reconstruction.reactions)), dtype=np.uint16)

    rxn_indices_to_keep: npt.NDArray[int] = np.unique(np.concatenate([rxns_from_imat, force_reaction_indices], dtype=int))

    # Reaction indices to exclude from the model are thus reactions that are not marked to be included in the model
    # Assume unique is false because every value that is in `rxn_indices_to_keep` is included in `all_rxn_indices`
    rxn_indices_to_remove: npt.NDArray[np.uint16] = np.setdiff1d(all_rxn_indices, rxn_indices_to_keep, assume_unique=False)
    model_reconstruction.remove_reactions(reactions=all_rxn_ids[rxn_indices_to_remove].tolist(), remove_orphans=True)

    return model_reconstruction


def _build_with_tinit(
    reference_model: cobra.Model,
    lower_bounds,
    upper_bounds,
    expr_vector,
    solver,
    idx_force,
) -> Model:
    raise NotImplementedError("tINIT is not yet implemented.")
    model = reference_model
    properties = tINITProperties(
        reactions_scores=expr_vector,
        solver=solver,
        essential_reactions=idx_force,
        production_weight=0.0,
        allow_excretion=False,
        no_reverse_loops=True,
    )
    model_reconstruction = reference_model.copy()
    s_matrix: npt.NDArray[float] = np.asarray(cobra.util.array.create_stoichiometric_matrix(model=model_reconstruction), dtype=float)
    algorithm = tINIT(s_matrix, lower_bounds, upper_bounds, properties)
    algorithm.preprocessing()
    algorithm.build_problem()
    _log_and_raise_error("tINIT is not yet implemented.", error=NotImplementedError, level=LogLevel.CRITICAL)

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
    match filepath.suffix:
        case ".mat":
            reference_model = cobra.io.load_matlab_model(filepath)
        case ".xml" | ".sbml":
            reference_model = cobra.io.read_sbml_model(filepath)
        case ".json":
            reference_model = cobra.io.load_json_model(filepath)
        case _:
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
    reference_model: cobra.Model = _read_reference_model(general_model_file)

    if objective not in force_reactions:
        force_reactions.append(objective)
    reference_model = _set_boundaries(reference_model, boundary_reactions, lower_bounds, upper_bounds)
    reference_model.solver = solver.lower()

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

    # get expressed reactions
    reaction_expression: collections.OrderedDict[str, int] = await _map_expression_to_reaction(
        reference_model, gene_expression_file, recon_algorithm, high_thresh=high_thresh, low_thresh=low_thresh, taxon=taxon
    )
    expression_vector: npt.NDArray[float] = np.array(list(reaction_expression.values()), dtype=float)

    for rxn in force_reactions:
        if rxn not in reaction_ids:
            logger.warning(
                f"The force reaction '{rxn}' was not found in the reference model. "
                f"Check BiGG, or the relevant database for your reference model, for synonyms."
            )

    # collect list of reactions that are infeasible but active in expression data or user defined
    infeasible_expression_reactions = []
    infeasible_force_reactions = []

    for i, rxn in enumerate(reaction_expression):
        # log reactions in expressed and force lists that are infeasible that the user may wish to review
        if rxn in inconsistent_reactions and expression_vector[i] == 1:
            infeasible_expression_reactions.append(rxn)
        if rxn in inconsistent_reactions and rxn in force_reactions:
            infeasible_force_reactions.append(rxn)

        if rxn in force_reactions:
            expression_vector[i] = high_thresh + 0.1 if recon_algorithm in {Algorithm.TINIT, Algorithm.IMAT} else 1
        if rxn in inconsistent_reactions or rxn in exclude_reactions:
            expression_vector[i] = low_thresh - 0.1 if recon_algorithm in {Algorithm.TINIT, Algorithm.IMAT} else 0

    objective_index = reaction_ids.index(objective)

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
            "minimum reaction rate",
            "maximum reaction rate",
        ]:
            raise ValueError(
                f"Boundary reactions file must have columns named 'Reaction', 'Abbreviation', 'Compartment', "
                f"'Minimum Reaction Rate', and 'Maximum Reaction Rate'. Found: {column}"
            )

    reactions: list[str] = [""] * len(df)
    boundary_type: list[str] = df["reaction"].tolist()
    reaction_abbreviation: list[str] = list(df["abbreviation"].astype(str))
    reaction_compartment: list[str] = list(df["compartment"].astype(str))
    boundary_map = {"exchange": "EX", "demand": "DM", "sink": "SK"}
    for i in range(len(boundary_type)):
        boundary: str = boundary_type[i].lower()
        if boundary not in boundary_map:
            raise ValueError(f"Boundary reaction type must be 'Exchange', 'Demand', or 'Sink'. Found: {boundary}")

        shorthand_compartment = CobraCompartments.get_shorthand(reaction_compartment[i])
        reactions[i] = f"{boundary_map.get(boundary)}_{reaction_abbreviation[i]}[{shorthand_compartment}]"

    return _BoundaryReactions(
        reactions=reactions,
        lower_bounds=df["minimum reaction rate"].tolist(),
        upper_bounds=df["maximum reaction rate"].tolist(),
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
    force_boundary_rxn_inclusion: bool = True,
):
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
        log_and_raise_error(
            f"Invalid output filetype. Should be 'xml', 'sbml', 'mat', or 'json'. Got:\n{invalid_suffix}'",
            error=ValueError,
            level=LogLevel.ERROR,
        raise ValueError(f"Invalid output filetype. Should be 'xml', 'sbml', 'mat', or 'json'. Got:\n{invalid_suffix}'")
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
