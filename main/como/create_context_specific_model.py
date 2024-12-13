from __future__ import annotations

import argparse
import ast
import collections
import re
from collections.abc import Sequence
from dataclasses import dataclass
from enum import Enum
from pathlib import Path
from typing import NamedTuple

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

from como.project import Config
from como.utils import Algorithm, Compartments, split_gene_expression_data, stringlist_to_list


class Solver(Enum):
    """Solver used to seed context specific model."""

    GLPK = "GLPK"
    GUROBI = "GUROBI"
    SCIPY = "SCIPY"
    GLPK_EXACT = "GLPK_EXACT"

    @staticmethod
    def from_string(value: str) -> Solver:
        """Convert string to Solver enum."""
        match value.lower():
            case "glpk":
                return Solver.GLPK
            case "gurobi":
                return Solver.GUROBI
            case "scipy":
                return Solver.SCIPY
            case "glpk_exact":
                return Solver.GLPK_EXACT
            case _:
                raise ValueError(f"Unknown solver: {value}")


class _BoundaryReactions(NamedTuple):
    """Boundary reactions to be used in the context specific model."""

    reactions: list[str]
    lower_bounds: list[float]
    upper_bounds: list[float]


class _BuildResults(NamedTuple):
    """Results of building a context specific model."""

    model: cobra.Model
    expression_index_list: list[int]
    infeasible_reactions: pd.DataFrame


@dataclass
class _Arguments:
    """Arguments for creating a context specific model."""

    context_name: str
    reference_model: Path
    active_genes_filepath: Path
    objective: str
    boundary_reactions_filepath: Path
    exclude_reactions_filepath: Path
    force_reactions_filepath: Path
    recon_algorithm: Algorithm
    solver: Solver
    low_threshold: int
    high_threshold: int
    output_filetypes: list[str]

    def __post_init__(self):
        self.reference_model = Path(self.reference_model)
        self.active_genes_filepath = Path(self.active_genes_filepath)
        self.boundary_reactions_filepath = (
            Path(self.boundary_reactions_filepath) if self.boundary_reactions_filepath is not None else None
        )
        self.exclude_reactions_filepath = (
            Path(self.exclude_reactions_filepath) if self.exclude_reactions_filepath is not None else None
        )
        self.force_reactions_filepath = (
            Path(self.force_reactions_filepath) if self.force_reactions_filepath is not None else None
        )

        if not self.reference_model.exists():
            raise FileNotFoundError(f"Reference model not found at {self.reference_model}")
        if not self.active_genes_filepath.exists():
            raise FileNotFoundError(f"Active genes file not found at {self.active_genes_filepath}")
        if self.boundary_reactions_filepath and not self.boundary_reactions_filepath.exists():
            raise FileNotFoundError(f"Boundary reactions file not found at {self.boundary_reactions_filepath}")
        if self.exclude_reactions_filepath and not self.exclude_reactions_filepath.exists():
            raise FileNotFoundError(f"Exclude reactions file not found at {self.exclude_reactions_filepath}")
        if self.force_reactions_filepath and not self.force_reactions_filepath.exists():
            raise FileNotFoundError(f"Force reactions file not found at {self.force_reactions_filepath}")

        if self.high_threshold < self.low_threshold:
            raise ValueError(
                f"Low threshold must be less than high threshold. "
                f"Received low threshold: {self.low_threshold}, high threshold: {self.high_threshold}"
            )


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
    new_left_rule = "".join(new_right_rule)
    final_right_rule = "" if rule_match is None else _correct_bracket(right_rule, right_name)
    return " ".join([new_left_rule, operator, final_right_rule])


def _gene_rule_logical(gpr_expression: str, level: int = 0) -> str:
    """Create an expression from GPR rule which can be evaluated as true or false."""
    try:
        loc_r = gpr_expression.index(")")
    except ValueError:
        if "and" in gpr_expression:
            gpr_expression = gpr_expression.replace("and", ",")
            return "min{" + gpr_expression + "}"
        elif "or" in gpr_expression:
            gpr_expression = gpr_expression.replace("or", ",")
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


def _gene_rule_evaluable(expression_in: str) -> str:
    """Make expression rule evaluable."""
    gene_reaction_by_rule = _gene_rule_logical(expression_in)
    gene_reaction_by_rule = gene_reaction_by_rule.replace("{", "(")
    gene_reaction_by_rule = gene_reaction_by_rule.replace("}", ")")

    return gene_reaction_by_rule


def _set_boundaries(model: cobra.Model, bound_rxns: list, bound_lb, bound_ub) -> tuple[cobra.Model, list]:
    all_rxns = model.reactions  # get all reactions
    bound_rm_rxns = []

    # get boundary reactions
    exchange_rxns = [rxn.id for rxn in all_rxns if re.search("EX_", rxn.id)]
    sink_rxns = [rxn.id for rxn in all_rxns if re.search("sink_", rxn.id)]
    demand_rxns = [rxn.id for rxn in all_rxns if re.search("DM_", rxn.id)]

    # Allows all boundary reactions to be used if none are given
    allow_all_boundary_rxns = not bound_rxns

    # close sinks and demands not in boundary reactions unless no boundary reactions were given
    if not allow_all_boundary_rxns:
        for i, rxn in enumerate(sink_rxns):  # set sinks to 0
            if rxn not in bound_rxns:  # only allow sink accumulation
                getattr(model.reactions, rxn).lower_bounds = 0
                getattr(model.reactions, rxn).upper_bounds = 1000
            else:  # set from file
                getattr(model.reactions, rxn).lower_bounds = bound_lb[i]
                getattr(model.reactions, rxn).upper_bounds = bound_ub[i]

        for i, rxn in enumerate(demand_rxns):
            getattr(model.reactions, rxn).lower_bounds = 0
            getattr(model.reactions, rxn).upper_bounds = bound_ub[i] if rxn in bound_rxns else 0

    # Reaction media
    medium = model.medium  # get reaction media to modify
    for rxn in exchange_rxns:  # open exchanges from exchange file, close unspecified exchanges
        if rxn not in bound_rxns:
            medium[rxn] = 0.0
        else:
            medium[rxn] = -float(bound_lb[bound_rxns.index(rxn)])
    model.medium = medium  # set new media

    return model, bound_rm_rxns


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


def _build_with_gimme(cobra_model, s_matrix, lb, ub, idx_objective, expr_vector):
    # `Becker and Palsson (2008). Context-specific metabolic networks are
    # consistent with experiments. PLoS Comput. Biol. 4, e1000082.`
    properties = GIMMEProperties(
        exp_vector=expr_vector,  # np.array(gimme_data['0']),
        obj_frac=0.9,
        objectives=[{idx_objective: 1}],
        preprocess=True,
        flux_threshold=0.9,
    )
    algorithm = GIMME(s_matrix, lb, ub, properties)
    gene_activity = algorithm.run()
    context_cobra_model = cobra_model.copy()
    r_ids = [r.id for r in context_cobra_model.reactions]
    to_remove_ids = [r_ids[r] for r in np.where(gene_activity == 0)[0]]

    context_cobra_model.remove_reactions(to_remove_ids, True)
    r_ids = [r.id for r in context_cobra_model.reactions]
    psol = pfba(context_cobra_model)  # noqa: F841
    # psol = context_cobra_model.optimize()
    # to_remove_ids = [r_ids[r] for r in np.where(abs(psol.fluxes) < 1e-8)[0]]
    # context_cobra_model.remove_reactions(to_remove_ids, True)

    return context_cobra_model


def _build_with_fastcore(cobra_model, s_matrix, lb, ub, exp_idx_list, solver):
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
    algorithm = FASTcore(s_matrix, lb, ub, properties)
    context_rxns = algorithm.fastcore()
    context_cobra_model = cobra_model.copy()
    r_ids = [r.id for r in context_cobra_model.reactions]
    remove_rxns = [r_ids[int(i)] for i in range(s_matrix.shape[1]) if i not in context_rxns]
    context_cobra_model.remove_reactions(remove_rxns, True)

    return context_cobra_model


def _build_with_imat(
    cobra_model: cobra.Model,
    s_matrix: npt.NDArray,
    lb: Sequence[float],
    ub: Sequence[float],
    expr_vector: npt.NDArray,
    expr_thesh: tuple[float, float],
    force_gene_ids: Sequence[int],
    solver: str,
) -> (cobra.Model, pd.DataFrame):
    expr_vector = np.array(expr_vector)
    properties = IMATProperties(
        exp_vector=expr_vector,
        exp_thresholds=expr_thesh,
        core=force_gene_ids,
        epsilon=0.01,
        solver=solver.upper(),
    )
    algorithm = IMAT(s_matrix, np.array(lb), np.array(ub), properties)
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


def _build_with_tinit(cobra_model: cobra.Model, s_matrix, lb, ub, expr_vector, solver, idx_force) -> Model:
    properties = tINITProperties(
        reactions_scores=expr_vector,
        solver=solver,
        essential_reactions=idx_force,
        production_weight=0.0,
        allow_excretion=False,
        no_reverse_loops=True,
    )
    algorithm = tINIT(s_matrix, lb, ub, properties)
    algorithm.preprocessing()
    algorithm.build_problem()
    raise NotImplementedError("tINIT is not yet implemented")


def _map_expression_to_reaction(
    model_cobra,
    gene_expression_file,
    recon_algorithm: Algorithm,
    low_thresh=None,
    high_thresh=None,
):
    """Map gene ids to a reaction based on GPR (gene to protein to reaction) association rules.

    These rules should be defined in the general genome-scale metabolic model
    """
    expression_data = pd.read_csv(gene_expression_file)
    gene_expressions = split_gene_expression_data(expression_data, recon_algorithm=recon_algorithm)
    expression_rxns = collections.OrderedDict()

    unknown_val = 1
    if recon_algorithm in {Algorithm.IMAT, Algorithm.TINIT}:
        unknown_val = np.mean([low_thresh, high_thresh])  # put unknowns in mid bin
    elif recon_algorithm == Algorithm.GIMME:
        unknown_val = -1
    elif recon_algorithm == Algorithm.FASTCORE:
        unknown_val = 0

    error_count = 0
    for rxn in model_cobra.reactions:
        gene_reaction_rule = _correct_bracket(rxn.gene_reaction_rule, rxn.gene_name_reaction_rule)
        gene_ids = re.findall(r"\d+", gene_reaction_rule)
        expression_rxns[rxn.id] = unknown_val
        if gene_reaction_rule.strip() == "":
            continue
        for gid in gene_ids:
            if gid in gene_expressions.index:
                rep_val = f' {gene_expressions.at[gid, "active"]} '
            else:
                rep_val = f" {unknown_val!s} "
            gene_reaction_rule = f" {gene_reaction_rule} "  # pad white space to prevent gene matches inside floats
            gene_reaction_rule = gene_reaction_rule.replace(f" {gid} ", rep_val, 1)
        try:
            gene_reaction_by_rule = _gene_rule_evaluable(gene_reaction_rule)
            gene_reaction_by_rule = gene_reaction_by_rule.strip()
            expression_rxns[rxn.id] = ast.literal_eval(gene_reaction_by_rule)

        except BaseException:
            error_count += 1

    logger.info(f"Mapped gene expression to reactions, found {error_count} error(s).")
    expr_vector = np.array(list(expression_rxns.values()), dtype=float)

    return expression_rxns, expr_vector


def _build_model(  # noqa: C901
    context_name: str,
    general_model_file: Path,
    gene_expression_file: Path,
    recon_algorithm: Algorithm,
    objective: str,
    bound_rxns: list[str],
    exclude_rxns: list[str],
    force_rxns: list[str],
    bound_lb: list[float],
    bound_ub: list[float],
    solver: str,
    low_thresh: float,
    high_thresh: float,
) -> _BuildResults:
    """Seed a context specific reference_model.

    Core reactions are determined from GPR associations with gene expression logicals.
    Core reactions that do not necessarily meet GPR association requirements can be forced if in the force reaction
    file. Metabolite exchange (media), sinks, and demands are determined from exchanges file. Reactions can also be
    force excluded even if they meet GPR association requirements using the force exclude file.
    """
    config = Config()
    reference_model: cobra.Model
    match general_model_file.suffix:
        case ".mat":
            reference_model = cobra.io.load_matlab_model(general_model_file)
        case ".xml":
            reference_model = cobra.io.read_sbml_model(general_model_file)
        case ".json":
            reference_model = cobra.io.load_json_model(general_model_file)
        case _:
            raise NameError(
                f"Reference reference_model format must be .xml, .mat, or .json, found '{general_model_file.suffix}'"
            )

    reference_model.objective = {getattr(reference_model.reactions, objective): 1}  # set objective

    if objective not in force_rxns:
        force_rxns.append(objective)

    # set boundaries
    reference_model, bound_rm_rxns = _set_boundaries(reference_model, bound_rxns, bound_lb, bound_ub)

    # set solver
    reference_model.solver = solver.lower()

    # check number of unsolvable reactions for reference model under media assumptions
    # incon_rxns, cobra_model = _feasibility_test(cobra_model, "before_seeding")
    incon_rxns = []

    s_matrix = cobra.util.array.create_stoichiometric_matrix(reference_model, array_type="dense")
    lb = []
    ub = []
    rx_names = []
    for reaction in reference_model.reactions:
        lb.append(reaction.lower_bound)
        ub.append(reaction.upper_bound)
        rx_names.append(reaction.id)

    # get expressed reactions
    expression_rxns, expr_vector = _map_expression_to_reaction(
        reference_model,
        gene_expression_file,
        recon_algorithm,
        high_thresh=high_thresh,
        low_thresh=low_thresh,
    )

    for rxn in force_rxns:
        if rxn not in rx_names:
            logger.warning(
                f"The force reaction '{rxn}' was not found in the general reference_model. "
                f"Check BiGG, or the relevant database for your general reference_model, for synonyms."
            )

    # collect list of reactions that are infeasible but active in expression data or user defined
    infeas_exp_rxns = []
    infeas_force_rxns = []
    infeas_exp_cnt = 0
    infeas_force_cnt = 0

    for idx, rxn in enumerate(expression_rxns):
        # log reactions in expressed and force lists that are infeasible that the user may wish to review
        if rxn in incon_rxns and expr_vector[idx] == 1:
            infeas_exp_cnt += 1
            infeas_exp_rxns.append(rxn)
        if rxn in incon_rxns and rxn in force_rxns:
            infeas_force_cnt += 1
            infeas_force_rxns.append(rxn)

        # make changes to expressed reactions base on user defined force/exclude reactions
        # TODO: if not using bound reactions file, add two sets of exchange reactions to be put in either low or mid bin

        if rxn in force_rxns:
            expr_vector[idx] = high_thresh + 0.1 if recon_algorithm.value in {"TINIT", "IMAT"} else 1
        if rxn in incon_rxns or rxn in exclude_rxns:
            expr_vector[idx] = low_thresh - 0.1 if recon_algorithm.value in {"TINIT", "IMAT"} else 0

    idx_obj = rx_names.index(objective)
    idx_force = [rx_names.index(rxn) for rxn in force_rxns if rxn in rx_names]
    exp_idx_list = [i for (i, val) in enumerate(expr_vector) if val > 0]  # type: ignore
    exp_thresh = (low_thresh, high_thresh)

    if recon_algorithm == Algorithm.GIMME:
        context_model_cobra = _build_with_gimme(reference_model, s_matrix, lb, ub, idx_obj, expr_vector)
    elif recon_algorithm == Algorithm.FASTCORE:
        context_model_cobra = _build_with_fastcore(reference_model, s_matrix, lb, ub, exp_idx_list, solver)
    elif recon_algorithm == Algorithm.IMAT:
        context_model_cobra: cobra.Model
        context_model_cobra, flux_df = _build_with_imat(
            reference_model,
            s_matrix,
            lb,
            ub,
            expr_vector,
            exp_thresh,
            idx_force,
            solver=solver,
        )
        imat_reactions = flux_df.rxn
        model_reactions = [reaction.id for reaction in context_model_cobra.reactions]
        reaction_intersections = set(imat_reactions).intersection(model_reactions)
        flux_df = flux_df[~flux_df["rxn"].isin(reaction_intersections)]
        flux_df.to_csv(config.data_dir / "results" / context_name / f"{recon_algorithm.value}_flux.csv")
    elif recon_algorithm == Algorithm.TINIT:
        context_model_cobra = _build_with_tinit(reference_model, s_matrix, lb, ub, expr_vector, solver, idx_force)

    incon_rxns_cs = []
    incon_df = pd.DataFrame({"general_infeasible_rxns": list(incon_rxns)})
    infeas_exp_df = pd.DataFrame({"expressed_infeasible_rxns": infeas_exp_rxns})
    infeas_force_df = pd.DataFrame({"infeasible_rxns_in_force_list": infeas_exp_rxns})
    incon_df_cs = pd.DataFrame({"infeasible_rxns_from_seeding": list(incon_rxns_cs)})
    infeasible_df = pd.concat([incon_df, infeas_exp_df, infeas_force_df, incon_df_cs], ignore_index=True, axis=1)
    infeasible_df.columns = ["InfeasRxns", "ExpressedInfeasRxns", "ForceInfeasRxns", "ContextInfeasRxns"]

    return _BuildResults(
        model=context_model_cobra,
        expression_index_list=exp_idx_list,
        infeasible_reactions=infeasible_df,
    )


def _create_df(path: Path) -> pd.DataFrame:
    match path.suffix:
        case ".csv":
            df = pd.read_csv(path, header=0, sep=",")
        case ".tsv":
            df = pd.read_csv(path, header=0, sep="\t")
        case ".xlsx" | ".xls":
            df = pd.read_excel(path, header=0)
        case _:
            raise FileNotFoundError(f"File not found! Must be a csv, tsv, or Excel file. Searching for: {path}")
    df.columns = [c.lower() for c in df.columns]
    return df


def _collect_boundary_reactions(path: Path) -> _BoundaryReactions:
    df = _create_df(path)
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
    reaction_abbreviation: list[str] = df["abbreviation"].tolist()
    reaction_compartment: list[str] = df["compartment"].tolist()
    lower_bound = df["minimum reaction rate"].tolist()
    upper_bound = df["maximum reaction rate"].tolist()
    boundary_map = {"exchange": "EX", "demand": "DM", "sink": "SK"}
    for i in range(len(boundary_type)):
        boundary: str = boundary_type[i].lower()
        if boundary not in boundary_map:
            raise ValueError(f"Boundary reaction type must be 'Exchange', 'Demand', or 'Sink'. Found: {boundary[i]}")

        shorthand_compartment = Compartments.get(reaction_compartment[i])
        reactions[i] = f"{boundary_map.get(boundary)}_{reaction_abbreviation[i]}[{shorthand_compartment}]"

    return _BoundaryReactions(
        reactions=reactions,
        lower_bounds=lower_bound,
        upper_bounds=upper_bound,
    )


def _write_model_to_disk(
    model: cobra.Model,
    output_directory: Path,
    context_name: str,
    output_filetypes: list[str],
    algorithm: Algorithm,
) -> None:
    if "mat" in output_filetypes:
        cobra.io.save_matlab_model(model, output_directory / f"{context_name}_SpecificModel_{algorithm.value}.mat")
    if "xml" in output_filetypes:
        cobra.io.write_sbml_model(model, output_directory / f"{context_name}_SpecificModel_{algorithm.value}.xml")
    if "json" in output_filetypes:
        cobra.io.save_json_model(model, output_directory / f"{context_name}_SpecificModel_{algorithm.value}.json")


def create_context_specific_model(  # noqa: C901
    context_name: str,
    reference_model: Path,
    genes_file: Path,
    objective: str = "biomass_reaction",
    boundary_rxns_filepath: str | Path | None = None,
    exclude_rxns_filepath: str | Path | None = None,
    force_rxns_filepath: str | Path | None = None,
    algorithm: Algorithm = Algorithm.GIMME,
    low_threshold: float = -5,
    high_threshold: float = -3,
    solver: Solver = Solver.GLPK,
    output_filetypes: list[str] | None = None,
):
    """Create a context-specific model using the provided data."""
    if not reference_model.exists():
        raise FileNotFoundError(f"Reference model not found at {reference_model}")
    if not genes_file.exists():
        raise FileNotFoundError(f"Active genes file not found at {genes_file}")
    if output_filetypes is None:
        output_filetypes = ["mat"]

    for output_type in output_filetypes:
        if output_type not in {"xml", "mat", "json"}:
            raise ValueError(f"Output file type {output_type} not recognized. Must be one of: 'xml', 'mat', 'json'")

    if algorithm not in Algorithm:
        raise ValueError(f"Algorithm {algorithm} not supported. Use one of {', '.join(a.value for a in Algorithm)}")

    if solver not in Solver:
        raise ValueError(f"Solver '{solver}' not supported. Use one of {', '.join(s.value for s in Solver)}")

    if boundary_rxns_filepath:
        boundary_reactions = _collect_boundary_reactions(boundary_rxns_filepath)

    exclude_rxns: list[str] = []
    if exclude_rxns_filepath:
        exclude_rxns_filepath: Path = Path(exclude_rxns_filepath)
        df = _create_df(exclude_rxns_filepath)
        if "abbreviation" not in df.columns:
            raise ValueError("The exclude reactions file should have a single column with a header named Abbreviation")
        exclude_rxns = df["abbreviation"].tolist()

    force_rxns: list[str] = []
    if force_rxns_filepath:
        force_rxns_filepath: Path = Path(force_rxns_filepath)
        df = _create_df(force_rxns_filepath)
        if "abbreviation" not in df.columns:
            raise ValueError("The force reactions file should have a single column with a header named Abbreviation")
        force_rxns = df["abbreviation"].tolist()

    logger.info(f"Creating '{context_name}' model using '{algorithm.value}' reconstruction and '{solver.value}' solver")
    build_results: _BuildResults = _build_model(
        context_name=context_name,
        general_model_file=reference_model,
        gene_expression_file=genes_file,
        recon_algorithm=algorithm,
        objective=objective,
        bound_rxns=boundary_reactions.reactions,
        bound_lb=boundary_reactions.lower_bounds,
        bound_ub=boundary_reactions.upper_bounds,
        exclude_rxns=exclude_rxns,
        force_rxns=force_rxns,
        solver=solver.value.lower(),
        low_thresh=low_threshold,
        high_thresh=high_threshold,
    )

    config = Config()
    build_results.infeasible_reactions.to_csv(
        config.result_dir / context_name / f"{context_name}_infeasible_rxns.csv", index=False
    )

    if algorithm == Algorithm.FASTCORE:
        pd.DataFrame(build_results.expression_index_list).to_csv(
            config.result_dir / context_name / f"{context_name}_core_rxns.csv", index=False
        )

    output_directory = config.result_dir / context_name
    _write_model_to_disk(
        model=build_results.model,
        output_directory=output_directory,
        context_name=context_name,
        output_filetypes=output_filetypes,
        algorithm=algorithm,
    )

    logger.success(f"Saved output file to {output_directory}")
    logger.info(f"Number of Genes: {len(build_results.model.genes):,}")
    logger.info(f"Number of Metabolites: {len(build_results.model.metabolites):,}")
    logger.info(f"Number of Reactions: {len(build_results.model.reactions):,}")


def _parse_args():
    parser = argparse.ArgumentParser(
        prog="create_context_specific_model.py",
        description="Seed a context-specific model from a list of expressed genes, a reference",
        epilog="For additional help, please post questions/issues in the MADRID GitHub repo at "
        "https://github.com/HelikarLab/MADRID or email babessell@gmail.com",
    )
    parser.add_argument(
        "-n",
        "--context-name",
        type=str,
        required=True,
        dest="context_name",
        help="Name of context or context used consistent with outputs of merge_xomics.py.",
    )
    parser.add_argument(
        "-m",
        "--reference-model-filepath",
        type=str,
        required=True,
        dest="reference_model",
        help="Name of Genome-scale metabolic model to seed the context model to. For example, the "
        "GeneralModelUpdatedV2.mat, is a modified Recon3D model. We also provide iMM_madrid for mouse."
        "OT can be .mat, .xml, or .json.",
    )
    parser.add_argument(
        "-g",
        "--active-genes-filepath",
        type=str,
        required=True,
        dest="active_genes_filepath",
        help="Path to logical table of active genes output from merge_xomics.py called "
        "ActiveGenes_contextName_Merged.csv. Should be in the corresponding context/context folder "
        "inside /main/data/results/contextName/. The json file output from the function using "
        "the context of interest as the key can be used here.",
    )
    parser.add_argument(
        "-o",
        "--objective",
        type=str,
        default="biomass_reaction",
        dest="objective",
        help="Reaction ID of the objective function in the model. Generally a biomass function.",
    )
    parser.add_argument(
        "-b",
        "--boundary-reactions-filepath",
        type=str,
        default=None,
        dest="boundary_reactions_filepath",
        help="Path to file contains the exchange (media), sink, and demand reactions which "
        "the model should use to fulfill the reactions governed by transcriptomic and proteomics "
        "data inputs. It must be a csv or xlsx with three columns: Rxn, Lowerbound, Upperbound. If not "
        "specified, MADRID will allow ALL BOUNDARY REACTIONS THAT ARE OPEN IN THE REFERENCE MODEL "
        "TO BE USED!",
    )
    parser.add_argument(
        "-x",
        "--exclude-reactions-filepath",
        type=str,
        default=None,
        dest="exclude_reactions_filepath",
        help="Filepath to file that contains reactions which will be removed from active reactions "
        "the model to use when seeding, even if considered active from transcriptomic and "
        "proteomics data inputs. It must be a csv or xlsx with one column of reaction IDs consistent with "
        "the reference model",
    )
    parser.add_argument(
        "-f",
        "--force-reactions-filepath",
        type=str,
        default=None,
        dest="force_reactions_filepath",
        help="Filepath to file that contains reactions which will be added to active reactions for "
        "the model to use when seeding (unless it causes the model to be unsolvable), regardless "
        "of results of transcriptomic and proteomics data inputs. It must be a csv or xlsx with one "
        "column of reaction IDs consistent with the reference model",
    )
    parser.add_argument(
        "-a",
        "--algorithm",
        type=str,
        default="GIMME",
        dest="recon_algorithm",
        help="Algorithm used to seed context specific model to the Genome-scale model. "
        "Can be either GIMME, FASTCORE, iMAT, or tINIT.",
    )
    parser.add_argument(
        "-lt",
        "--low-threshold",
        type=float,
        default=-5,
        dest="low_threshold",
        help="Low to mid bin cutoff for iMAT solution",
    )
    parser.add_argument(
        "-ht",
        "--high-threshold",
        type=float,
        default=-3,
        dest="high_threshold",
        help="Mid to high bin cutoff for iMAT solution",
    )
    parser.add_argument(
        "-s",
        "--solver",
        type=str,
        default="glpk",
        dest="solver",
        help="Solver used to seed model and attempt to solve objective. Default is GLPK, also takes "
        "GUROBI but you must mount a container license to the Docker to use. An academic license "
        "can be obtained for free. See the README on the Github or Dockerhub for information on "
        "mounting this license.",
    )
    parser.add_argument(
        "-t",
        "--output-filetypes",
        type=str,
        nargs="+",
        default="mat",
        dest="output_filetypes",
        help="Filetypes to save seeded model type. Can be either a string with one filetype such as "
        "'xml' or multiple in the format \"['extension1', 'extension2', ... etc]\". If you want "
        "to output in all 3 accepted formats,  would be: \"['mat', 'xml', 'json']\" "
        "Note the outer quotes required to be interpreted by cmd. This a string, not a python list",
    )
    args = parser.parse_args()
    args.output_filetypes = stringlist_to_list(args.output_filetypes)
    args.solver = Solver.from_string(args.solver)  # type: ignore
    args.recon_algorithm = Algorithm.from_string(args.recon_algorithm)  # type: ignore
    return _Arguments(**vars(args))


if __name__ == "__main__":
    args = _parse_args()
    create_context_specific_model(
        context_name=args.context_name,
        reference_model=args.reference_model,
        genes_file=args.active_genes_filepath,
        objective=args.objective,
        boundary_rxns_filepath=args.boundary_reactions_filepath,
        exclude_rxns_filepath=args.exclude_reactions_filepath,
        force_rxns_filepath=args.force_reactions_filepath,
        algorithm=args.recon_algorithm,
        low_threshold=args.low_threshold,
        high_threshold=args.high_threshold,
        solver=args.solver,
        output_filetypes=args.output_filetypes,
    )
