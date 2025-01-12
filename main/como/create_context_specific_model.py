from __future__ import annotations

import collections
import re
import sys
from collections.abc import Sequence
from io import TextIOWrapper
from pathlib import Path

import cobra
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

from como.data_types import Algorithm, CobraCompartments, LogLevel, Solver, _BoundaryReactions, _BuildResults
from como.utils import _log_and_raise_error, _read_file, _set_up_logging, split_gene_expression_data


def _correct_bracket(rule: str, name: str) -> str:
    """Correct GPR rules to format readable by."""
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
    """Create an expression from GPR rule which can be evaluated as true or false."""
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

    expression_out = f"{gpr_expression[:loc_l]}{inner_string}{gpr_expression[loc_r + 1:]}"
    expression_out = _gene_rule_logical(expression_out, level + 1)

    return expression_out


def _set_boundaries(
    model: cobra.Model,
    boundary_reactions: list[str],
    lower_bounds: list[float],
    upper_bounds: list[float],
) -> cobra.Model:
    # get boundary reactions
    exchange_rxns = [rxn.id for rxn in model.reactions if re.search("EX_", rxn.id)]
    sink_rxns = [rxn.id for rxn in model.reactions if re.search("sink_", rxn.id)]
    demand_rxns = [rxn.id for rxn in model.reactions if re.search("DM_", rxn.id)]

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


def _build_with_gimme(cobra_model, s_matrix, lower_bounds, upper_bounds, idx_objective, expr_vector):
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
    context_cobra_model = cobra_model.copy()
    reaction_ids = [r.id for r in context_cobra_model.reactions]
    to_remove_ids = [reaction_ids[r] for r in np.where(gene_activity == 0)[0]]

    context_cobra_model.remove_reactions(to_remove_ids, True)
    psol = pfba(context_cobra_model)  # noqa: F841
    # reaction_ids = [r.id for r in context_cobra_model.reactions]
    # psol = context_cobra_model.optimize()
    # to_remove_ids = [reaction_ids[r] for r in np.where(abs(psol.fluxes) < 1e-8)[0]]
    # context_cobra_model.remove_reactions(to_remove_ids, True)

    return context_cobra_model


def _build_with_fastcore(cobra_model, s_matrix, lower_bounds, upper_bounds, exp_idx_list, solver):
    # 'Vlassis, Pacheco, Sauter (2014). Fast reconstruction of compact
    # context-specific metabolic network models. PLoS Comput. Biol. 10,
    # e1003424.'
    logger.warning(
        "Fastcore requires a flux consistant model is used as refererence, "
        "to achieve this fastcc is required which is NOT reproducible."
    )
    logger.debug("Creating feasible model")
    incon_rxns, cobra_model = _feasibility_test(cobra_model, "other")
    properties = FastcoreProperties(core=exp_idx_list, solver=solver)
    algorithm = FASTcore(s_matrix, lower_bounds, upper_bounds, properties)
    context_rxns = algorithm.fastcore()
    context_cobra_model = cobra_model.copy()
    r_ids = [r.id for r in context_cobra_model.reactions]
    remove_rxns = [r_ids[int(i)] for i in range(s_matrix.shape[1]) if i not in context_rxns]
    context_cobra_model.remove_reactions(remove_rxns, True)

    return context_cobra_model


def _build_with_imat(
    cobra_model: cobra.Model,
    s_matrix: npt.NDArray,
    lower_bounds: Sequence[float],
    upper_bounds: Sequence[float],
    expr_vector: npt.NDArray,
    expr_thresh: tuple[float, float],
    force_gene_ids: Sequence[int],
    solver: str,
) -> (cobra.Model, pd.DataFrame):
    expr_vector = np.array(expr_vector)
    properties = IMATProperties(
        exp_vector=expr_vector,
        exp_thresholds=expr_thresh,
        core=force_gene_ids,
        epsilon=0.01,
        solver=solver.upper(),
    )
    algorithm = IMAT(s_matrix, np.array(lower_bounds), np.array(upper_bounds), properties)
    context_rxns: npt.NDArray = algorithm.run()
    fluxes: pd.Series = algorithm.sol.to_series()
    context_cobra_model = cobra_model.copy()
    reaction_ids = [r.id for r in context_cobra_model.reactions]

    remove_rxns = [reaction_ids[i] for i in range(s_matrix.shape[1]) if i not in context_rxns]
    flux_df = pd.DataFrame(columns=["rxn", "flux"])
    for idx, (_, val) in enumerate(fluxes.items()):
        if idx <= len(cobra_model.reactions) - 1:
            r_id = str(context_cobra_model.reactions.get_by_id(reaction_ids[idx])).split(":")[0]
            getattr(context_cobra_model.reactions, r_id).fluxes = val
            flux_df.loc[len(flux_df.index)] = [r_id, val]

    context_cobra_model.remove_reactions(remove_rxns, True)

    return context_cobra_model, flux_df


def _build_with_tinit(
    cobra_model: cobra.Model,
    s_matrix,
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
    algorithm = tINIT(s_matrix, lower_bounds, upper_bounds, properties)
    algorithm.preprocessing()
    algorithm.build_problem()
    _log_and_raise_error("tINIT is not yet implemented.", error=NotImplementedError, level=LogLevel.CRITICAL)


async def _map_expression_to_reaction(
    model_cobra,
    gene_expression_file,
    recon_algorithm: Algorithm,
    low_thresh: float,
    high_thresh: float,
) -> collections.OrderedDict[str, int]:
    """Map gene ids to a reaction based on GPR (gene to protein to reaction) association rules.

    These rules should be defined in the general genome-scale metabolic model
    """
    gene_activity = split_gene_expression_data(await _read_file(gene_expression_file), recon_algorithm=recon_algorithm)
    reaction_expression = collections.OrderedDict()

    # fmt: off
    # Define a default expression value if a gene ID is not found
    default_expression = (
        np.mean([low_thresh, high_thresh]) if recon_algorithm in {Algorithm.IMAT, Algorithm.TINIT}
        else -1 if recon_algorithm in {Algorithm.GIMME}
        else 0 if recon_algorithm in {Algorithm.FASTCORE}
        else 1
    )
    # fmt: on

    error_count = 0
    for rxn in model_cobra.reactions:
        rxn: cobra.Reaction
        gene_reaction_rule = _correct_bracket(rxn.gene_reaction_rule, rxn.gene_name_reaction_rule)
        if gene_reaction_rule == "":
            continue

        gene_ids = re.findall(r"\d+", gene_reaction_rule)
        reaction_expression[rxn.id] = default_expression
        for gene_id in gene_ids:
            activity = (
                f"{gene_activity.at[gene_id, 'active']}"
                if gene_id in gene_activity.index
                else f"{default_expression!s}"
            )
            # replace gene_id with activity, using optional whitespace before and after the gene id
            # Do not replace the whitespace (if it exists) before and after the gene ID
            gene_reaction_rule = re.sub(
                pattern=rf"(?<!\S){gene_id}(?!\S)",
                repl=activity,
                string=gene_reaction_rule,
                count=1,
            )

        try:
            # We are using eval here because ast.literal_eval is unable to process an evaluable such as `max(-4, -4)`
            # This isn't ideal, but ultimately the only other option is writing and maintaining a custom parsing engine,
            #   which is too much work to do.
            evaluable_gene_rule = _gene_rule_logical(gene_reaction_rule).replace("{", "(").replace("}", ")")
            reaction_expression[rxn.id] = eval(evaluable_gene_rule)  # noqa: S307
        except ValueError:
            error_count += 1

    logger.debug(f"Mapped gene expression to reactions, found {error_count} error(s).")
    expr_vector = np.array(list(reaction_expression.values()), dtype=float)

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
) -> _BuildResults:
    """Seed a context specific reference_model.

    Core reactions are determined from GPR associations with gene expression logicals.
    Core reactions that do not necessarily meet GPR association requirements can be forced if in the force reaction
    file. Metabolite exchange (media), sinks, and demands are determined from exchanges file. Reactions can also be
    force excluded even if they meet GPR association requirements using the force exclude file.
    """
    reference_model: cobra.Model
    match general_model_file.suffix:
        case ".mat":
            reference_model = cobra.io.load_matlab_model(general_model_file)
        case ".xml":
            reference_model = cobra.io.read_sbml_model(general_model_file)
        case ".json":
            reference_model = cobra.io.load_json_model(general_model_file)
        case _:
            _log_and_raise_error(
                f"Reference model format must be .xml, .mat, or .json; found '{general_model_file.suffix}'",
                error=ValueError,
                level=LogLevel.ERROR,
            )

    reference_model.objective = {getattr(reference_model.reactions, objective): 1}
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
                f"The force reaction '{rxn}' was not found in the general reference_model. "
                f"Check BiGG, or the relevant database for your general reference_model, for synonyms."
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
    force_reaction_indices = [reaction_ids.index(rxn) for rxn in force_reactions if rxn in reaction_ids]
    expression_vector_indices = [i for (i, val) in enumerate(expression_vector) if val > 0]  # type: ignore
    expression_threshold = (low_thresh, high_thresh)

    match recon_algorithm:
        case Algorithm.GIMME:
            context_model_cobra = _build_with_gimme(
                cobra_model=reference_model,
                s_matrix=s_matrix,
                lower_bounds=lower_bounds,
                upper_bounds=upper_bounds,
                idx_objective=objective_index,
                expr_vector=expression_vector,
            )
        case Algorithm.FASTCORE:
            context_model_cobra = _build_with_fastcore(
                cobra_model=reference_model,
                s_matrix=s_matrix,
                lower_bounds=lower_bounds,
                upper_bounds=upper_bounds,
                exp_idx_list=expression_vector_indices,
                solver=solver,
            )
        case Algorithm.IMAT:
            context_model_cobra: cobra.Model
            context_model_cobra, flux_df = _build_with_imat(
                cobra_model=reference_model,
                s_matrix=s_matrix,
                lower_bounds=lower_bounds,
                upper_bounds=upper_bounds,
                expr_vector=expression_vector,
                expr_thresh=expression_threshold,
                force_gene_ids=force_reaction_indices,
                solver=solver,
            )
            imat_reactions = flux_df.rxn
            model_reactions = [reaction.id for reaction in context_model_cobra.reactions]
            reaction_intersections = set(imat_reactions).intersection(model_reactions)
            flux_df: pd.DataFrame = flux_df[~flux_df["rxn"].isin(reaction_intersections)]
            flux_df.dropna(inplace=True)
            flux_df.to_csv(output_flux_result_filepath)
        case Algorithm.TINIT:
            context_model_cobra = _build_with_tinit(
                cobra_model=reference_model,
                s_matrix=s_matrix,
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

    return _BuildResults(
        model=context_model_cobra,
        expression_index_list=expression_vector_indices,
        infeasible_reactions=inconsistent_and_infeasible_reactions,
    )


async def _create_df(path: Path) -> pd.DataFrame:
    match path.suffix:
        case ".csv" | ".tsv":
            df = await _read_file(path, header=0, sep="," if path.suffix == ".csv" else "\t")
        case ".xlsx" | ".xls":
            df = await _read_file(path, header=0)
        case _:
            _log_and_raise_error(
                f"File not found! Must be a csv, tsv, or Excel file. Searching for: {path}",
                error=FileNotFoundError,
                level=LogLevel.ERROR,
            )
    df.columns = [c.lower() for c in df.columns]
    return df


async def _collect_boundary_reactions(path: Path) -> _BoundaryReactions:
    df: pd.DataFrame = await _create_df(path)
    for column in df.columns:
        if column not in [
            "boundary",
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
    boundary_type: list[str] = df["boundary"].tolist()
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

    return _BoundaryReactions(
        reactions=reactions,
        lower_bounds=df["minimum reaction rate"].tolist(),
        upper_bounds=df["maximum reaction rate"].tolist(),
    )


async def _write_model_to_disk(
    model: cobra.Model,
    output_directory: Path,
    context_name: str,
    output_filetypes: list[str],
    algorithm: Algorithm,
) -> None:
    output_directory.mkdir(parents=True, exist_ok=True)
    if "mat" in output_filetypes:
        cobra.io.save_matlab_model(model, (output_directory / f"{context_name}_SpecificModel_{algorithm.value}.mat"))
    if "json" in output_filetypes:
        cobra.io.save_json_model(model, (output_directory / f"{context_name}_SpecificModel_{algorithm.value}.json"))
    if "xml" in output_filetypes:
        cobra.io.write_sbml_model(model, (output_directory / f"{context_name}_SpecificModel_{algorithm.value}.xml"))


async def create_context_specific_model(  # noqa: C901
    context_name: str,
    reference_model: Path,
    active_genes_filepath: Path,
    output_infeasible_reactions_filepath: Path,
    output_model_dirpath: Path,
    output_flux_result_filepath: Path,
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
    log_location: str | TextIOWrapper = sys.stderr,
):
    """Create a context-specific model using the provided data."""
    _set_up_logging(level=log_level, location=log_location)
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

    logger.info(f"Creating '{context_name}' model using '{algorithm.value}' reconstruction and '{solver.value}' solver")
    build_results: _BuildResults = await _build_model(
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
    )

    build_results.infeasible_reactions.dropna(inplace=True)
    build_results.infeasible_reactions.to_csv(output_infeasible_reactions_filepath, index=False)

    if algorithm == Algorithm.FASTCORE:
        fastcore_df = pd.DataFrame(build_results.expression_index_list)
        fastcore_df.dropna(inplace=True)
        fastcore_df.to_csv(output_fastcore_expression_index_filepath, index=False)

    await _write_model_to_disk(
        model=build_results.model,
        output_directory=output_model_dirpath,
        context_name=context_name,
        output_filetypes=output_filetypes,
        algorithm=algorithm,
    )

    logger.success(f"Saved metabolic model for context '{context_name}' to {output_model_dirpath}")
    logger.debug(f"Number of Genes: {len(build_results.model.genes):,}")
    logger.debug(f"Number of Metabolites: {len(build_results.model.metabolites):,}")
    logger.debug(f"Number of Reactions: {len(build_results.model.reactions):,}")
