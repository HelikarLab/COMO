from __future__ import annotations

import collections
import re
import sys
from collections.abc import Sequence
from io import TextIOWrapper
from pathlib import Path
from typing import Literal, TextIO, cast

import cobra
import cobra.util.array
import numpy as np
import numpy.typing as npt
import pandas as pd
from cobra import Model
from cobra.flux_analysis import pfba
from loguru import logger
from troppo.methods.reconstruction.fastcore import FASTcore, FastcoreProperties
from troppo.methods.reconstruction.gimme import GIMME, GIMMEProperties
from troppo.methods.reconstruction.imat import IMAT, IMATProperties
from troppo.methods.reconstruction.tINIT import tINIT, tINITProperties

from como.data_types import Algorithm, BoundaryReactions, BuildResults, CobraCompartments, LogLevel, Solver
from como.utils import _log_and_raise_error, read_file, set_up_logging, split_gene_expression_data


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

    new_right_rule = []
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


def _set_boundaries(
    model: cobra.Model,
    boundary_reactions: list[str],
    lower_bounds: list[float],
    upper_bounds: list[float],
) -> cobra.Model:
    # get boundary reactions
    exchange_rxns = [rxn.id for rxn in model.reactions if "EX_" in rxn.id]
    sink_rxns = [rxn.id for rxn in model.reactions if "sink_" in rxn.id]
    demand_rxns = [rxn.id for rxn in model.reactions if "DM_" in rxn.id]

    # Allows all boundary reactions to be used if none are given
    allow_all_boundary_rxns = not boundary_reactions

    # close sinks and demands not in boundary reactions unless no boundary reactions were given
    if not allow_all_boundary_rxns:
        for i, rxn in enumerate(sink_rxns):  # set sinks to 0
            getattr(model.reactions, rxn).lower_bounds = lower_bounds[i] if rxn in boundary_reactions else 0
            getattr(model.reactions, rxn).upper_bounds = upper_bounds[i] if rxn in boundary_reactions else 1000

        for i, rxn in enumerate(demand_rxns):
            getattr(model.reactions, rxn).lower_bounds = 0
            getattr(model.reactions, rxn).upper_bounds = upper_bounds[i] if rxn in boundary_reactions else 0

    # Reaction media
    medium = model.medium
    for rxn in exchange_rxns:  # open exchanges from exchange file, close unspecified exchanges
        medium[rxn] = -float(lower_bounds[boundary_reactions.index(rxn)]) if rxn in boundary_reactions else 0.0
    model.medium = medium

    return model


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
    lower_bounds: Sequence[float | np.floating],
    upper_bounds: Sequence[float | np.floating],
    idx_objective: int,
    expr_vector: npt.NDArray[np.floating],
):
    model_reconstruction = reference_model.copy()
    s_matrix: npt.NDArray[np.floating] = cobra.util.array.create_stoichiometric_matrix(model=model_reconstruction)
    # `Becker and Palsson (2008). Context-specific metabolic networks are
    # consistent with experiments. PLoS Comput. Biol. 4, e1000082.`
    properties = GIMMEProperties(
        exp_vector=expr_vector,  # np.array(gimme_data['0']),
        obj_frac=0.9,
        objectives=[{idx_objective: 1}],
        preprocess=True,
        flux_threshold=0.9,
    )
    algorithm = GIMME(s_matrix, lower_bounds, upper_bounds, properties)
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


def _build_with_fastcore(cobra_model, s_matrix, lower_bounds, upper_bounds, exp_idx_list, solver):
    # 'Vlassis, Pacheco, Sauter (2014). Fast reconstruction of compact
    # context-specific metabolic network models. PLoS Comput. Biol. 10,
    # e1003424.'
    logger.warning("Fastcore requires a flux consistant model is used as refererence, to achieve this fastcc is required which is NOT reproducible.")
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
    lower_bounds: Sequence[float],
    upper_bounds: Sequence[float],
    expr_vector: npt.NDArray,
    expr_thresh: tuple[float, float],
    force_gene_indices: Sequence[int],
    solver: str,
) -> cobra.Model:
    properties: IMATProperties = IMATProperties(
        exp_vector=expr_vector,
        exp_thresholds=expr_thresh,
        core=force_gene_indices,
        epsilon=0.01,
        solver=solver.upper(),
    )

    # Creating a copy of the model ensures we don't make any in-place modifications by accident
    # Using cobra to create the stoichiometry matrix means we have less work to do
    force_gene_indices = np.array(force_gene_indices, dtype=np.uint16)
    model_reconstruction: cobra.Model = reference_model.copy()
    s_matrix: npt.NDArray[np.floating] = cobra.util.array.create_stoichiometric_matrix(model=model_reconstruction)
    algorithm: IMAT = IMAT(S=s_matrix, lb=np.array(lower_bounds), ub=np.array(upper_bounds), properties=properties)
    rxns_from_imat: npt.NDArray[np.uint16] = algorithm.run().astype(np.uint16)

    # Collect all reaction IDs and their associated index (e.g., HEX1 is at index 123)
    all_rxn_ids: npt.NDArray[str] = np.array([r.id for r in model_reconstruction.reactions], dtype=object)
    all_rxn_indices: npt.NDArray[np.uint16] = np.array(range(len(model_reconstruction.reactions)), dtype=np.uint16)

    # Collect reactions to keep by creating a unique set of reactions from the iMAT algorithm and force-include reactions
    # dtype is set to uint16 because indices will not be below 0 or be greater than 65,535 (max size of uint16),
    #   because only ~10,000 reactions exist in Recon3D
    # Unsafe casting is OK because of these facts.
    rxn_indices_to_keep: npt.NDArray[np.uint16] = np.unique(np.concatenate([rxns_from_imat, force_gene_indices], dtype=np.uint16))

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
    properties = tINITProperties(
        reactions_scores=expr_vector,
        solver=solver,
        essential_reactions=idx_force,
        production_weight=0.0,
        allow_excretion=False,
        no_reverse_loops=True,
    )
    model_reconstruction = reference_model.copy()
    s_matrix: npt.NDArray[np.floating] = cobra.util.array.create_stoichiometric_matrix(model=model_reconstruction).astype(np.float64)
    algorithm = tINIT(s_matrix, lower_bounds, upper_bounds, properties)
    algorithm.preprocessing()
    algorithm.build_problem()
    _log_and_raise_error("tINIT is not yet implemented.", error=NotImplementedError, level=LogLevel.CRITICAL)


async def _map_expression_to_reaction(
    reference_model,
    gene_expression_file,
    recon_algorithm: Algorithm,
    low_thresh: float,
    high_thresh: float,
) -> collections.OrderedDict[str, int]:
    """Map gene ids to a reaction based on GPR (gene to protein to reaction) association rules.

    These rules should be defined in the general genome-scale metabolic model

    Args:
        reference_model: A COBRA model object representing the general genome-scale metabolic model.
        gene_expression_file: Path to a gene expression file (.csv, .tsv, .xlsx, or .xls)
        recon_algorithm: Algorithm to use for reconstruction (GIMME, FASTCORE, iMAT, or tINIT)
        low_thresh: Low expression threshold for algorithms that require it (iMAT, tINIT)
        high_thresh: High expression threshold for algorithms that require it (iMAT, tINIT)

    Returns:
        An ordered dictionary mapping reaction IDs to their corresponding expression values.

    Raises:
        ValueError: If neither 'entrez_gene_id' nor 'ensembl_gene_id' columns are found in the gene expression file.
    """
    expression_data = await read_file(gene_expression_file)
    identifier_column = next((col for col in ("entrez_gene_id", "ensembl_gene_id") if col in expression_data.columns), "")

    if not identifier_column:
        raise ValueError(
            f"At least one column of 'entrez_gene_id' or 'ensembl_gene_id' could not be found in the gene expression file '{gene_expression_file}'"
        )
    gene_activity = split_gene_expression_data(
        expression_data,
        identifier_column=cast(Literal["ensembl_gene_id", "entrez_gene_id"], identifier_column),
        recon_algorithm=recon_algorithm,
    )
    reaction_expression = collections.OrderedDict()

    # fmt: off
    # Define a default expression value if a gene ID is not found
    default_expression = (
        np.mean([low_thresh, high_thresh]) if recon_algorithm in {Algorithm.IMAT, Algorithm.TINIT}
        else -1 if recon_algorithm == Algorithm.GIMME
        else 0 if recon_algorithm == Algorithm.FASTCORE
        else 1
    )
    # fmt: on

    error_count = 0
    for rxn in reference_model.reactions:
        rxn: cobra.Reaction

        gene_reaction_rule = rxn.gene_reaction_rule
        if not gene_reaction_rule:
            continue

        gene_ids = set(re.findall(r"\d+", gene_reaction_rule))
        reaction_expression[rxn.id] = default_expression
        for gene_id in gene_ids:
            activity = gene_activity.at[gene_id, "active"] if gene_id in gene_activity.index else f"{default_expression!s}"
            # replace gene_id with activity, using optional whitespace before and after the gene id
            # Do not replace the whitespace (if it exists) before and after the gene ID
            gene_reaction_rule = re.sub(pattern=rf"\b{gene_id}\b", repl=activity, string=gene_reaction_rule)

        try:
            # We are using eval here because ast.literal_eval is unable to process an evaluable such as `max(-4, 0, 1)`
            # This isn't ideal, but ultimately the only other option is writing and maintaining a custom parsing engine, which is too much work
            evaluable_gene_rule = _gene_rule_logical(gene_reaction_rule).replace("{", "(").replace("}", ")")
            reaction_expression[rxn.id] = eval(evaluable_gene_rule)  # noqa: S307
        except ValueError:
            error_count += 1

    logger.debug(f"Mapped gene expression to reactions, found {error_count} error(s).")
    # expr_vector = np.array(list(reaction_expression.values()), dtype=float)

    return reaction_expression


async def _build_model(  # noqa: C901
    general_model_file: Path,
    gene_expression_file: Path,
    recon_algorithm: Algorithm,
    objective: str,
    boundary_reactions: list[str],
    exclude_reactions: list[str],
    force_reactions: list[str],
    lower_bounds: list[float],
    upper_bounds: list[float],
    solver: str,
    low_thresh: float,
    high_thresh: float,
    output_flux_result_filepath: Path,
    *,
    force_boundary_rxn_inclusion: bool,
) -> BuildResults:
    """Seed a context specific reference_model.

    Core reactions are determined from GPR associations with gene expression logicals.
    Core reactions that do not necessarily meet GPR association requirements can be forced if in the force reaction
    file. Metabolite exchange (media), sinks, and demands are determined from exchanges file. Reactions can also be
    force excluded even if they meet GPR association requirements using the force exclude file.

    Args:
        general_model_file: Path to a COBRA model file (.xml, .mat, or .json)
        gene_expression_file: Path to a gene expression file (.csv, .tsv, .xlsx, or .xls)
        recon_algorithm: Algorithm to use for reconstruction (GIMME, FASTCORE, iMAT, or tINIT)
        objective: Objective reaction ID in the general model
        boundary_reactions: List of boundary reactions to set in the model
        exclude_reactions: List of reactions to exclude from the model
        force_reactions: List of reactions to force include in the model
        lower_bounds: List of lower bounds corresponding to boundary reactions
        upper_bounds: List of upper bounds corresponding to boundary reactions
        solver: Solver to use (e.g., 'glpk', 'cplex', 'gurobi')
        low_thresh: Low expression threshold for algorithms that require it (iMAT, tINIT)
        high_thresh: High expression threshold for algorithms that require it (iMAT, tINIT)
        output_flux_result_filepath: Path to save flux results (for iMAT only)
        force_boundary_rxn_inclusion: If True, ensure that all boundary reactions are included in the final model.

    Returns:
        A _BuildResults object containing the context-specific model, list of expression indices used, and a DataFrame of infeasible reactions.
    """
    reference_model: cobra.Model
    match general_model_file.suffix:
        case ".mat":
            reference_model = cobra.io.load_matlab_model(general_model_file)
        case (".xml", ".sbml"):
            reference_model = cobra.io.read_sbml_model(general_model_file)
        case ".json":
            reference_model = cobra.io.load_json_model(general_model_file)
        case _:
            _log_and_raise_error(
                f"Reference model format must be .xml, .mat, or .json; found '{general_model_file.suffix}'",
                error=ValueError,
                level=LogLevel.ERROR,
            )

    if objective not in force_reactions:
        force_reactions.append(objective)
    reference_model = _set_boundaries(reference_model, boundary_reactions, lower_bounds, upper_bounds)
    reference_model.solver = solver.lower()

    # check number of unsolvable reactions for reference model under media assumptions
    # inconsistent_reactions, cobra_model = _feasibility_test(cobra_model, "before_seeding")
    inconsistent_reactions = []
    s_matrix = cobra.util.array.create_stoichiometric_matrix(reference_model, array_type="dense")
    lower_bounds = []
    upper_bounds = []
    reaction_ids = []
    for reaction in reference_model.reactions:
        lower_bounds.append(reaction.lower_bound)
        upper_bounds.append(reaction.upper_bound)
        reaction_ids.append(reaction.id)

    # get expressed reactions
    reaction_expression: collections.OrderedDict[str, int] = await _map_expression_to_reaction(
        reference_model,
        gene_expression_file,
        recon_algorithm,
        high_thresh=high_thresh,
        low_thresh=low_thresh,
    )
    expression_vector: npt.NDArray[np.int32] = np.array(list(reaction_expression.values()), dtype=np.int32)

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

    match recon_algorithm:
        case Algorithm.GIMME:
            context_model_cobra: cobra.Model = _build_with_gimme(
                reference_model=reference_model,
                lower_bounds=lower_bounds,
                upper_bounds=upper_bounds,
                idx_objective=objective_index,
                expr_vector=expression_vector,
            )
        case Algorithm.FASTCORE:
            context_model_cobra: cobra.Model = _build_with_fastcore(
                cobra_model=reference_model,
                s_matrix=s_matrix,
                lower_bounds=lower_bounds,
                upper_bounds=upper_bounds,
                exp_idx_list=expression_vector_indices,
                solver=solver,
            )
        case Algorithm.IMAT:
            context_model_cobra: cobra.Model = _build_with_imat(
                reference_model=reference_model,
                lower_bounds=lower_bounds,
                upper_bounds=upper_bounds,
                expr_vector=expression_vector,
                expr_thresh=expression_threshold,
                force_gene_indices=force_reaction_indices,
                solver=solver,
            )
            context_model_cobra.objective = objective
            flux_sol: cobra.Solution = context_model_cobra.optimize()
            fluxes: pd.Series = flux_sol.fluxes
            model_reactions: list[str] = [reaction.id for reaction in context_model_cobra.reactions]
            reaction_intersections: set[str] = set(fluxes.index).intersection(model_reactions)
            flux_df: pd.DataFrame = fluxes[~fluxes.index.isin(reaction_intersections)]
            flux_df.dropna(inplace=True)
            flux_df.to_csv(output_flux_result_filepath)
        case Algorithm.TINIT:
            context_model_cobra: cobra.Model = _build_with_tinit(
                reference_model=reference_model,
                lower_bounds=lower_bounds,
                upper_bounds=upper_bounds,
                expr_vector=expression_vector,
                solver=solver,
                idx_force=force_reaction_indices,
            )
        case _:
            _log_and_raise_error(
                (
                    f"Reconstruction algorithm must be {Algorithm.GIMME.value}, "
                    f"{Algorithm.FASTCORE.value}, {Algorithm.IMAT.value}, or {Algorithm.TINIT.value}. "
                    f"Got: {recon_algorithm.value}"
                ),
                error=ValueError,
                level=LogLevel.ERROR,
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

    return BuildResults(
        model=context_model_cobra,
        expression_index_list=expression_vector_indices,
        infeasible_reactions=inconsistent_and_infeasible_reactions,
    )


async def _create_df(path: Path, *, lowercase_col_names: bool = False) -> pd.DataFrame:
    if path.suffix not in {".csv", ".tsv"}:
        raise ValueError(f"File must be a .csv or .tsv file, got '{path.suffix}'")
    df: pd.DataFrame = await read_file(path=path, header=0, sep="," if path.suffix == ".csv" else "\t", h5ad_as_df=True)

    if not isinstance(df, pd.DataFrame):
        _log_and_raise_error(
            f"Expected a pandas.DataFrame, got {type(df)}",
            error=TypeError,
            level=LogLevel.ERROR,
        )

    if lowercase_col_names:
        df.columns = [c.lower() for c in df.columns]
    return df


async def _collect_boundary_reactions(path: Path) -> BoundaryReactions:
    df: pd.DataFrame = await _create_df(path, lowercase_col_names=True)
    for column in df.columns:
        if column not in [
            "reaction",
            "abbreviation",
            "compartment",
            "minimum reaction rate",
            "maximum reaction rate",
        ]:
            _log_and_raise_error(
                (
                    f"Boundary reactions file must have columns named 'Reaction', 'Abbreviation', 'Compartment', "
                    f"'Minimum Reaction Rate', and 'Maximum Reaction Rate'. Found: {column}"
                ),
                error=ValueError,
                level=LogLevel.ERROR,
            )

    reactions: list[str] = [""] * len(df)
    boundary_type: list[str] = df["reaction"].tolist()
    reaction_abbreviation: list[str] = df["abbreviation"].astype(str).tolist()
    reaction_compartment: list[str] = df["compartment"].astype(str).tolist()
    boundary_map = {"exchange": "EX", "demand": "DM", "sink": "SK"}
    for i in range(len(boundary_type)):
        boundary: str = boundary_type[i].lower()
        if boundary not in boundary_map:
            _log_and_raise_error(
                f"Boundary reaction type must be 'Exchange', 'Demand', or 'Sink'. Found: {boundary}",
                error=ValueError,
                level=LogLevel.ERROR,
            )

        shorthand_compartment = CobraCompartments.get_shorthand(reaction_compartment[i])
        reactions[i] = f"{boundary_map.get(boundary)}_{reaction_abbreviation[i]}[{shorthand_compartment}]"

    return BoundaryReactions(
        reactions=reactions,
        lower_bounds=df["minimum reaction rate"].tolist(),
        upper_bounds=df["maximum reaction rate"].tolist(),
    )


def _write_model_to_disk(context_name: str, model: cobra.Model, output_filepaths: list[Path]) -> None:
    for path in output_filepaths:
        path.parent.mkdir(parents=True, exist_ok=True)
        if path.suffix == ".mat":
            cobra.io.save_matlab_model(model=model, file_name=path)
        elif path.suffix == ".json":
            cobra.io.save_json_model(model=model, filename=path, pretty=True)
        elif path.suffix in {".sbml", ".xml"}:
            cobra.io.write_sbml_model(cobra_model=model, filename=path)
        else:
            _log_and_raise_error(
                f"Invalid output model filetype. Should be one of .xml, .sbml, .mat, or .json. Got '{path.suffix}'",
                error=ValueError,
                level=LogLevel.ERROR,
            )
        logger.success(f"Saved metabolic model for context '{context_name}' to '{path}'")


async def create_context_specific_model(  # noqa: C901
    context_name: str,
    reference_model: Path,
    active_genes_filepath: Path,
    output_infeasible_reactions_filepath: Path,
    output_flux_result_filepath: Path,
    output_model_filepaths: Path | list[Path],
    output_filetypes: list[str] | None = None,
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
        reference_model: Path to the general genome-scale metabolic model file (.xml, .mat, or .json).
        active_genes_filepath: Path to the gene expression data file (csv, tsv, or Excel).
        output_infeasible_reactions_filepath: Path to save infeasible reactions (csv).
        output_flux_result_filepath: Path to save flux results (csv).
        output_model_filepaths: Path or list of paths to save the context-specific model (.xml, .mat, or .json).
        output_filetypes: List of file types to save the model as ('xml', 'mat', 'json').
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
    _set_up_logging(level=log_level, location=log_location)
    output_model_filepaths = [output_model_filepaths] if isinstance(output_model_filepaths, Path) else output_model_filepaths
    for path in output_model_filepaths:
        if path.suffix not in {".mat", ".xml", ".sbml", ".json"}:
            _log_and_raise_error(
                f"Invalid output model filetype. Should be one of .xml, .sbml, .mat, or .json. Got '{path.suffix}'",
                error=ValueError,
                level=LogLevel.ERROR,
            )
    if len(output_model_filepaths) != len(output_model_filepaths):
        _log_and_raise_error(
            "The number of output model filepaths must be the same as the number of output flux result filepaths",
            error=ValueError,
            level=LogLevel.ERROR,
        )

    if not reference_model.exists():
        _log_and_raise_error(
            f"Reference model not found at {reference_model}",
            error=FileNotFoundError,
            level=LogLevel.ERROR,
        )
    if not active_genes_filepath.exists():
        _log_and_raise_error(
            f"Active genes file not found at {active_genes_filepath}",
            error=FileNotFoundError,
            level=LogLevel.ERROR,
        )
    if algorithm == Algorithm.FASTCORE and not output_fastcore_expression_index_filepath:
        _log_and_raise_error(
            "The fastcore expression index output filepath must be provided",
            error=ValueError,
            level=LogLevel.ERROR,
        )
    if boundary_rxns_filepath and not boundary_rxns_filepath.exists():
        _log_and_raise_error(
            f"Boundary reactions file not found at {boundary_rxns_filepath}",
            error=FileNotFoundError,
            level=LogLevel.ERROR,
        )

    output_filetypes = ["mat"] if output_filetypes is None else output_filetypes
    for output_type in output_filetypes:
        if output_type not in {"xml", "mat", "json"}:
            _log_and_raise_error(
                f"Output file type {output_type} not recognized. Must be one of: 'xml', 'mat', 'json'",
                error=ValueError,
                level=LogLevel.ERROR,
            )

    if algorithm not in Algorithm:
        _log_and_raise_error(
            f"Algorithm {algorithm} not supported. Use one of {', '.join(a.value for a in Algorithm)}",
            error=ValueError,
            level=LogLevel.ERROR,
        )

    if solver not in Solver:
        _log_and_raise_error(
            f"Solver '{solver}' not supported. Use one of {', '.join(s.value for s in Solver)}",
            error=ValueError,
            level=LogLevel.ERROR,
        )

    if boundary_rxns_filepath:
        boundary_reactions = await _collect_boundary_reactions(boundary_rxns_filepath)

    exclude_rxns: list[str] = []
    if exclude_rxns_filepath:
        exclude_rxns_filepath: Path = Path(exclude_rxns_filepath)
        df = await _create_df(exclude_rxns_filepath)
        if "abbreviation" not in df.columns:
            _log_and_raise_error(
                "The exclude reactions file should have a single column with a header named Abbreviation",
                error=ValueError,
                level=LogLevel.ERROR,
            )
        exclude_rxns = df["abbreviation"].tolist()

    force_rxns: list[str] = []
    if force_rxns_filepath:
        force_rxns_filepath: Path = Path(force_rxns_filepath)
        df = await _create_df(force_rxns_filepath)
        if "abbreviation" not in df.columns:
            _log_and_raise_error(
                "The force reactions file should have a single column with a header named Abbreviation",
                error=ValueError,
                level=LogLevel.ERROR,
            )
        force_rxns = df["abbreviation"].tolist()

    # Test that gurobi is using a valid license file
    if solver == Solver.GUROBI:
        # test if gurobi is available
        try:
            import gurobipy as gp
        except ImportError as e:
            logger.error(
                "The gurobi solver requires the gurobipy package to be installed. "
                "Please install gurobipy and try again. "
                "This can be done by installing the 'gurobi' optional dependency."
            )
            raise ImportError from e

        env = gp.Env()
        if env.getParam("WLSACCESSID") == "" or env.getParam("WLSSECRET") == "":
            logger.critical(
                "Gurobi solver requested, but license information cannot be found. "
                "COMO will continue, but it is HIGHLY unlikely the resulting model will be valid."
            )
        # remove gurobi-related information, it is no longer required
        del env, gp

    logger.info(f"Creating '{context_name}' model using '{algorithm.value}' reconstruction and '{solver.value}' solver")
    build_results: BuildResults = await _build_model(
        general_model_file=reference_model,
        gene_expression_file=active_genes_filepath,
        recon_algorithm=algorithm,
        objective=objective,
        boundary_reactions=boundary_reactions.reactions,
        exclude_reactions=exclude_rxns,
        force_reactions=force_rxns,
        lower_bounds=boundary_reactions.lower_bounds,
        upper_bounds=boundary_reactions.upper_bounds,
        solver=solver.value.lower(),
        low_thresh=low_threshold,
        high_thresh=high_threshold,
        output_flux_result_filepath=output_flux_result_filepath,
        force_boundary_rxn_inclusion=force_boundary_rxn_inclusion,
    )

    build_results.infeasible_reactions.dropna(inplace=True)
    build_results.infeasible_reactions.to_csv(output_infeasible_reactions_filepath, index=False)

    if algorithm == Algorithm.FASTCORE:
        fastcore_df = pd.DataFrame(build_results.expression_index_list)
        fastcore_df.dropna(inplace=True)
        fastcore_df.to_csv(output_fastcore_expression_index_filepath, index=False)

    _write_model_to_disk(context_name=context_name, model=build_results.model, output_filepaths=output_model_filepaths)
    logger.debug(f"Number of Genes: {len(build_results.model.genes):,}")
    logger.debug(f"Number of Metabolites: {len(build_results.model.metabolites):,}")
    logger.debug(f"Number of Reactions: {len(build_results.model.reactions):,}")
