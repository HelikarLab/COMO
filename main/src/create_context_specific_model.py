#!/usr/bin/python3

import os
import re
import sys
import cobra
from cobra import Model
import argparse
import collections
import numpy as np
import pandas as pd
from pathlib import Path
from warnings import warn
from cobra.flux_analysis import pfba
from cobamp.wrappers import COBRAModelObjectReader
from troppo.methods.reconstruction.imat import IMAT, IMATProperties
from troppo.methods.reconstruction.tINIT import tINIT, tINITProperties
from troppo.methods.reconstruction.gimme import GIMME, GIMMEProperties
from troppo.methods.reconstruction.fastcore import FASTcore, FastcoreProperties

from project import Configs
from como_utilities import stringlist_to_list, split_gene_expression_data, Compartments

sys.setrecursionlimit(1500)  # for re.search

configs = Configs()


def correct_bracket(rule: str, name: str) -> str:
    """
    Correct GPR rules to format readable by
    """
    rmatch = re.search(r"or|and", rule)
    nmatch = re.search(r"or|and", name)
    if rmatch is None or nmatch is None:
        lrule = rule
        lname = name.strip()
        rrule = ""
        rname = ""
        operator = ""
    else:
        lrule = rule[0: rmatch.span()[0]]
        lname = name[0: nmatch.span()[0]].strip()
        rrule = rule[rmatch.span()[1]:]
        rname = name[nmatch.span()[1]:]
        operator = rmatch.group()
    
    rlist_new = []
    for ch in list(lrule):
        if ch.isspace() or ch.isdigit():  # re.match(r'\w+', ch)
            rlist_new.append(ch)
        elif len(lname) > 0:
            if ch == lname[0]:
                rlist_new.append(ch)
                lname = lname[1:]
    rule_left = "".join(rlist_new)
    
    if rmatch is None:
        rule_right = ""
    else:
        rule_right = correct_bracket(rrule, rname)
    
    return " ".join([rule_left, operator, rule_right])


def gene_rule_logical(expression_in: str, level: int = 0) -> str:
    """
    creates an expression from GPR rule which can be evaluated as true or false
    """
    try:
        loc_r = expression_in.index(")")
    except BaseException:
        if "and" in expression_in:
            expression_in = expression_in.replace("and", ",")
            expression_in = "min{" + expression_in + "}"
        elif "or" in expression_in:
            expression_in = expression_in.replace("or", ",")
            expression_in = "max{" + expression_in + "}"
        else:
            expression_in = expression_in.replace("[", "")
            expression_in = expression_in.replace("]", "")
        
        return expression_in
    
    loc_l = expression_in[0:loc_r].rindex("(")
    inner_string = expression_in[loc_l: loc_r + 1]
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
    
    expression_out = "{}{}{}".format(
        expression_in[0:loc_l], inner_string, expression_in[loc_r + 1:]
    )
    expression_out = gene_rule_logical(expression_out, level + 1)
    
    return expression_out


def gene_rule_evaluable(expression_in: str) -> str:
    """
    Make expression rule evaluable
    """
    gene_reaction_by_rule = gene_rule_logical(expression_in)
    gene_reaction_by_rule = gene_reaction_by_rule.replace("{", "(")
    gene_reaction_by_rule = gene_reaction_by_rule.replace("}", ")")
    
    return gene_reaction_by_rule


def set_boundaries(model_cobra: cobra.Model, bound_rxns: list, bound_lb, bound_ub) -> tuple[cobra.Model, list]:
    all_rxns = model_cobra.reactions  # get all reactions
    bound_rm_rxns = []
    
    # get boundary reactions
    exchange_rxns = [rxn.id for rxn in all_rxns if re.search("EX_", rxn.id)]
    sink_rxns = [rxn.id for rxn in all_rxns if re.search("sink_", rxn.id)]
    demand_rxns = [rxn.id for rxn in all_rxns if re.search("DM_", rxn.id)]
    
    # set flag that allows all boundary reactions to be used if none are given
    if len(bound_rxns) == 0:
        allow_all_boundary_rxns = True
    else:
        allow_all_boundary_rxns = False
    
    # close sinks and demands not in boundary reactions unless no boundary reactions were given
    if not allow_all_boundary_rxns:
        
        for i, rxn in enumerate(sink_rxns):  # set sinks to 0
            if rxn not in bound_rxns:  # only allow sink accumulation
                getattr(model_cobra.reactions, rxn).lower_bound = 0
                getattr(model_cobra.reactions, rxn).upper_bound = 1000
            else:  # set from file
                getattr(model_cobra.reactions, rxn).lower_bound = bound_lb[i]
                getattr(model_cobra.reactions, rxn).upper_bound = bound_ub[i]
        
        for i, rxn in enumerate(demand_rxns):
            if rxn not in bound_rxns:  # demand is one way - outside the system
                getattr(model_cobra.reactions, rxn).lower_bound = 0
                getattr(model_cobra.reactions, rxn).upper_bound = 1000
            else:
                getattr(model_cobra.reactions, rxn).lower_bound = 0
                getattr(model_cobra.reactions, rxn).upper_bound = bound_ub[i]
    
    # Reaction media
    medium = model_cobra.medium  # get reaction media to modify
    for rxn in exchange_rxns:  # open exchanges from exchange file, close unspecified exchanges
        if rxn not in bound_rxns:
            medium[rxn] = 0.0
        else:
            medium[rxn] = -float(bound_lb[bound_rxns.index(rxn)])
    model_cobra.medium = medium  # set new media
    
    return model_cobra, bound_rm_rxns


def feasibility_test(model_cobra: cobra.Model, step: str):
    # check number of unsolvable reactions for reference model under media assumptions
    model_cobra_rm = cobra.flux_analysis.fastcc(model_cobra, flux_threshold=15,
                                                zero_cutoff=1e-7)  # create flux consistant model (rmemoves some reactions)
    incon_rxns = set(model_cobra.reactions.list_attr("id")) - set(model_cobra_rm.reactions.list_attr("id"))
    incon_rxns_cnt = len(incon_rxns)
    
    if step == "before_seeding":
        print(
            f"Under given boundary assumptions, there are {str(incon_rxns_cnt)} infeasible reactions in the "
            "reference model.\n"
        )
        print(
            "These reactions will not be considered active in context specific model construction. If any infeasible "
            "reactions are found to be active according to expression data, or, are found in the force reactions "
            "list, they can be found found in 'InfeasibleRxns.csv'\n"
        )
        print(
            "It is normal for this value to be quite large, however, if many of these reactions are active "
            "according to your expression data, it is likely that you are missing some critical exchange (media) "
            "reactions.\n"
        )
    elif step == "after_seeding":
        print(
            f"Under given boundary assumptions, with infeasible reactions from the general model not considered "
            f"there are {incon_rxns_cnt} new infeasible reactions in the context-specific model.\n"
        )
        print(
            "These reactions will be removed from the output model to ensure the model is solvable"
        )
        print(
            "Note that this value should be very low compared to the reference model."
        )
    else:
        pass
    
    return incon_rxns, model_cobra_rm


def seed_gimme(cobra_model, s_matrix, lb, ub, idx_objective, expr_vector):
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
    psol = pfba(context_cobra_model)
    # psol = context_cobra_model.optimize()
    # to_remove_ids = [r_ids[r] for r in np.where(abs(psol.fluxes) < 1e-8)[0]]
    # context_cobra_model.remove_reactions(to_remove_ids, True)
    
    return context_cobra_model


def seed_fastcore(cobra_model, s_matrix, lb, ub, exp_idx_list, solver):
    # 'Vlassis, Pacheco, Sauter (2014). Fast reconstruction of compact
    # context-specific metabolic network models. PLoS Comput. Biol. 10,
    # e1003424.'
    warn(
        "Fastcore requires a flux consistant model is used as refererence, "
        "to achieve this fastcc is required which is NOT reproducible."
    )
    print("Creating feasible model...")
    incon_rxns, cobra_model = feasibility_test(cobra_model, "other")
    properties = FastcoreProperties(core=exp_idx_list, solver=solver)
    algorithm = FASTcore(s_matrix, lb, ub, properties)
    context_rxns = algorithm.fastcore()
    context_cobra_model = cobra_model.copy()
    r_ids = [r.id for r in context_cobra_model.reactions]
    remove_rxns = [
        r_ids[int(i)] for i in range(s_matrix.shape[1]) if i not in context_rxns
    ]
    context_cobra_model.remove_reactions(remove_rxns, True)
    
    return context_cobra_model


def seed_imat(
    cobra_model, s_matrix, lb, ub, expr_vector, expr_thesh, idx_force, context_name
):
    expr_vector = np.array(expr_vector)
    properties = IMATProperties(
        exp_vector=expr_vector,
        exp_thresholds=expr_thesh,
        core=idx_force,
        epsilon=0.01
    )
    print("Setting properties")
    algorithm = IMAT(s_matrix, lb, ub, properties)
    
    print("Setting algorithm")
    context_rxns = algorithm.run()
    
    print("Running")
    fluxes = algorithm.sol.to_series()
    
    print("Obtained flux values")
    context_cobra_model = cobra_model.copy()
    r_ids = [r.id for r in context_cobra_model.reactions]
    pd.DataFrame({"rxns": r_ids}).to_csv(os.path.join(configs.data_dir, "rxns_test.csv"))
    remove_rxns = [r_ids[int(i)] for i in range(s_matrix.shape[1]) if not np.isin(i, context_rxns)]
    flux_df = pd.DataFrame(columns=["rxn", "flux"])
    for idx, (_, val) in enumerate(fluxes.items()):
        if idx <= len(cobra_model.reactions) - 1:
            r_id = str(context_cobra_model.reactions.get_by_id(r_ids[idx])).split(":")[0]
            getattr(context_cobra_model.reactions, r_id).fluxes = val
            flux_df.loc[len(flux_df.index)] = [r_id, val]
    
    context_cobra_model.remove_reactions(remove_rxns, True)
    
    return context_cobra_model, flux_df


def seed_tinit(
    cobra_model: cobra.Model,
    s_matrix,
    lb,
    ub,
    expr_vector,
    solver,
    idx_force
) -> Model:
    expr_vector = np.array(expr_vector)
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
    return Model()


def map_expression_to_rxn(model_cobra, gene_expression_file, recon_algorithm, low_thresh=None, high_thresh=None):
    """
    Map gene ids to a reaction based on GPR (gene to protein to reaction) association rules
    which are defined in general genome-scale metabolic model
    """
    expression_data = pd.read_csv(gene_expression_file)
    gene_expressions = split_gene_expression_data(
        expression_data, recon_algorithm=recon_algorithm
    )
    expression_rxns = collections.OrderedDict()
    
    cnt = 0
    if recon_algorithm in ["IMAT", "TINIT"]:
        # unknown_val = min(gene_expressions["Data"].tolist())
        unknown_val = np.mean([low_thresh, high_thresh])  # put unknowns in mid bin
    elif recon_algorithm == "GIMME":
        unknown_val = -1
    elif recon_algorithm == "FASTCORE":
        unknown_val = 0
    else:
        unknown_val = 1
    
    test_counter = 0
    for rxn in model_cobra.reactions:
        test_counter += 1
        gene_reaction_rule = correct_bracket(
            rxn.gene_reaction_rule, rxn.gene_name_reaction_rule
        )
        gene_ids = re.findall(r"\d+", gene_reaction_rule)
        expression_rxns[rxn.id] = unknown_val
        if gene_reaction_rule.strip() == "":
            continue
        for gid in gene_ids:
            if gid in gene_expressions.index:
                rep_val = " {} ".format(gene_expressions.at[gid, "Data"])
            else:
                rep_val = f" {str(unknown_val)} "
            gene_reaction_rule = (
                " " + gene_reaction_rule + " "
            )  # pad white space to prevent gene matches inside floats
            gene_reaction_rule = gene_reaction_rule.replace(
                " {} ".format(gid), rep_val, 1
            )
        try:
            gene_reaction_by_rule = gene_rule_evaluable(gene_reaction_rule)
            gene_reaction_by_rule = gene_reaction_by_rule.strip()
            expression_rxns[rxn.id] = eval(gene_reaction_by_rule)
        
        except BaseException:
            cnt += 1
    
    print("Map gene expression to reactions, {} errors.".format(cnt))
    expr_vector = np.array(list(expression_rxns.values()), dtype=float)
    
    return expression_rxns, expr_vector


def create_context_specific_model(
    general_model_file,
    gene_expression_file,
    recon_algorithm,
    objective,
    exclude_rxns,
    force_rxns,
    bound_rxns,
    bound_lb,
    bound_ub,
    solver,
    context_name,
    low_thresh,
    high_thresh
):
    """
    Seed a context specific model. Core reactions are determined from GPR associations with gene expression logicals.
    Core reactions that do not necessarily meet GPR association requirements can be forced if in the force reaction
    file. Metabolite exchange (media), sinks, and demands are determined from exchanges file. Reactions can also be
    force excluded even if they meet GPR association requirements using the force exclude file.
    """
    if general_model_file[-4:] == ".mat":
        cobra_model = cobra.io.load_matlab_model(general_model_file)
    elif general_model_file[-4:] == ".xml":
        cobra_model = cobra.io.load_sbml_model(general_model_file)
    elif general_model_file[-5:] == ".json":
        cobra_model = cobra.io.load_json_model(general_model_file)
    else:
        raise NameError("reference model format must be .xml, .mat, or .json")
    
    cobra_model.objective = {getattr(cobra_model.reactions, objective): 1}  # set objective
    
    if objective not in force_rxns:
        force_rxns.append(objective)
    
    # set boundaries
    cobra_model, bound_rm_rxns = set_boundaries(
        cobra_model, bound_rxns, bound_lb, bound_ub
    )
    
    # set solver
    cobra_model.solver = solver.lower()
    
    # check number of unsolvable reactions for reference model under media assumptions
    # incon_rxns, cobra_model = feasibility_test(cobra_model, "before_seeding")
    incon_rxns = []
    
    # CoBAMP methods
    cobamp_model = COBRAModelObjectReader(cobra_model)  # load model in readable format for CoBAMP
    cobamp_model.get_irreversibilities(True)
    s_matrix = cobamp_model.get_stoichiometric_matrix()
    lb, ub = cobamp_model.get_model_bounds(False, True)
    rx_names = cobamp_model.get_reaction_and_metabolite_ids()[0]
    
    # get expressed reactions
    expression_rxns, expr_vector = map_expression_to_rxn(
        cobra_model,
        gene_expression_file,
        recon_algorithm,
        high_thresh=high_thresh,
        low_thresh=low_thresh
    )
    
    # find reactions in the force reactions file that are not in general model and warn user
    for rxn in force_rxns:
        if rxn not in rx_names:
            warn(
                f"{rxn} from force reactions not in general model, check BiGG (or relevant database used in your "
                "general model) for synonyms"
            )
    
    # collect list of reactions that are infeasible but active in expression data or user defined
    infeas_exp_rxns = []
    infeas_force_rxns = []
    infeas_exp_cnt = 0
    infeas_force_cnt = 0
    
    for idx, (rxn, exp) in enumerate(expression_rxns.items()):
        # log reactions in expressed and force lists that are infeasible that the user may wish to review
        if rxn in incon_rxns and expr_vector[idx] == 1:
            infeas_exp_cnt += 1
            infeas_exp_rxns.append(rxn)
        if rxn in incon_rxns and rxn in force_rxns:
            infeas_force_cnt += 1
            infeas_force_rxns.append(rxn)
        
        # make changes to expressed reactions base on user defined force/exclude reactions
        # TODO: if not using bound reactions file, add two sets of exchange reactoins to be put in either low or mid bin
        
        if rxn in force_rxns:
            expr_vector[idx] = (
                high_thresh + 0.1 if recon_algorithm in ["TINIT", "IMAT"] else 1
            )
        if rxn in incon_rxns or rxn in exclude_rxns:
            expr_vector[idx] = (
                low_thresh - 0.1 if recon_algorithm in ["TINIT", "IMAT"] else 0
            )
    
    idx_obj = rx_names.index(objective)
    idx_force = [rx_names.index(rxn) for rxn in force_rxns if rxn in rx_names]
    exp_idx_list = [idx for (idx, val) in enumerate(expr_vector) if val > 0]
    exp_thresh = (low_thresh, high_thresh)
    
    # switch case dictionary runs the functions making it too slow, better solution then elif ladder?
    if recon_algorithm == "GIMME":
        context_model_cobra = seed_gimme(
            cobra_model, s_matrix, lb, ub, idx_obj, expr_vector
        )
    elif recon_algorithm == "FASTCORE":
        context_model_cobra = seed_fastcore(
            cobra_model, s_matrix, lb, ub, exp_idx_list, solver
        )
    elif recon_algorithm == "IMAT":
        flux_df: pd.DataFrame
        context_model_cobra: cobra.Model
        context_model_cobra, flux_df = seed_imat(
            cobra_model,
            s_matrix,
            lb,
            ub,
            expr_vector,
            exp_thresh,
            idx_force,
            context_name,
        )
        imat_reactions = flux_df.rxn
        model_reactions = [reaction.id for reaction in context_model_cobra.reactions]
        reaction_intersections = set(imat_reactions).intersection(model_reactions)
        flux_df = flux_df[~flux_df["rxn"].isin(reaction_intersections)]
        flux_df.to_csv(str(os.path.join(configs.data_dir, "results", context_name, f"{recon_algorithm}_flux.csv")))
    
    elif recon_algorithm == "TINIT":
        context_model_cobra = seed_tinit(
            cobra_model, s_matrix, lb, ub, expr_vector, solver, idx_force
        )
    else:
        raise ValueError(
            f"Unsupported reconstruction algorithm: {recon_algorithm}. Must be 'IMAT', 'GIMME', or 'FASTCORE'")
    
    # check number of unsolvable reactions for seeded model under media assumptions
    # incon_rxns_cs, context_model_cobra = feasibility_test(context_model_cobra, "after_seeding")
    incon_rxns_cs = []
    
    # if recon_algorithm in ["IMAT"]:
    #     final_rxns = [rxn.id for rxn in context_model_cobra.reactions]
    #     imat_rxns = flux_df.rxn
    #     for rxn in imat_rxns:
    #         if rxn not in final_rxns:
    #             flux_df = flux_df[flux_df.rxn != rxn]
    #
    #     flux_df.to_csv(os.path.join(configs.datadir, "results", context_name, f"{recon_algorithm}_flux.csv"))
    
    incon_df = pd.DataFrame({"general_infeasible_rxns": list(incon_rxns)})
    infeas_exp_df = pd.DataFrame({"expressed_infeasible_rxns": infeas_exp_rxns})
    infeas_force_df = pd.DataFrame({"infeasible_rxns_in_force_list": infeas_exp_rxns})
    incon_df_cs = pd.DataFrame({"infeasible_rxns_from_seeding": list(incon_rxns_cs)})
    infeasible_df = pd.concat(
        [incon_df, infeas_exp_df, infeas_force_df, incon_df_cs],
        ignore_index=True,
        axis=1,
    )
    infeasible_df.columns = [
        "InfeasRxns",
        "ExpressedInfeasRxns",
        "ForceInfeasRxns",
        "ContextInfeasRxns",
    ]
    
    return context_model_cobra, exp_idx_list, infeasible_df


def print_filetype_help():
    print("Unsupported model format. Current support is for: 'xml', 'mat', and 'json'."
          "Or use multiple with: 'xml mat json'")


def parse_args(argv):
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
        dest="modelfile",
        help="Name of Genome-scale metabolic model to seed the context model to. For example, the "
             "GeneralModelUpdatedV2.mat, is a modified Recon3D model. We also provide iMM_madrid for mouse."
             "OT can be .mat, .xml, or .json.",
    )
    parser.add_argument(
        "-g",
        "--active-genes-filepath",
        type=str,
        required=True,
        dest="genefile",
        help="Path to logical table of active genes output from merge_xomics.py called "
             "ActiveGenes_contextName_Merged.csv. Should be in the corresponding context/context folder "
             "inside /main/data/results/contextName/. The json file output from the function using "
             "the context of interest as the key can be used here.",
    )
    parser.add_argument(
        "-o",
        "--objective",
        type=str,
        required=False,
        default="biomass_reaction",
        dest="objective",
        help="Reaction ID of the objective function in the model. Generally a biomass function.",
    )
    parser.add_argument(
        "-b",
        "--boundary-reactions-filepath",
        type=str,
        required=False,
        default=None,
        dest="boundary_rxns_filepath",
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
        required=False,
        default=None,
        dest="exclude_rxns_filepath",
        help="Filepath to file that contains reactions which will be removed from active reactions "
             "the model to use when seeding, even if considered active from transcriptomic and "
             "proteomics data inputs. It must be a csv or xlsx with one column of reaction IDs consistent with "
             "the reference model",
    )
    parser.add_argument(
        "-f",
        "--force-reactions-filepath",
        type=str,
        required=False,
        default=None,
        dest="force_rxns_filepath",
        help="Filepath to file that contains reactions which will be added to active reactions for "
             "the model to use when seeding (unless it causes the model to be unsolvable), regardless "
             "of results of transcriptomic and proteomics data inputs. It must be a csv or xlsx with one "
             "column of reaction IDs consistent with the reference model",
    )
    parser.add_argument(
        "-a",
        "--algorithm",
        type=str,
        required=False,
        default="GIMME",
        dest="algorithm",
        help="Algorithm used to seed context specific model to the Genome-scale model. Can be either "
             "GIMME, FASTCORE, iMAT, or tINIT.",
    )
    parser.add_argument(
        "-lt",
        "--low-threshold",
        type=float,
        required=False,
        default=-5,
        dest="low_threshold",
        help="Low to mid bin cutoff for iMAT solution"
    )
    parser.add_argument(
        "-ht",
        "--high-threshold",
        type=float,
        required=False,
        default=-3,
        dest="high_threshold",
        help="Mid to high bin cutoff for iMAT solution"
    )
    parser.add_argument(
        "-s",
        "--solver",
        type=str,
        required=False,
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
        required=False,
        default="mat",
        dest="output_filetypes",
        help="Filetypes to save seeded model type. Can be either a string with one filetype such as "
             "'xml' or multiple in the format \"['extension1', 'extension2', ... etc]\". If you want "
             "to output in all 3 accepted formats,  would be: \"['mat', 'xml', 'json']\" "
             "Note the outer quotes required to be interpreted by cmd. This a string, not a python list",
    )
    args = parser.parse_args()
    return args


def main(argv):
    """
    Seed a context-specific model from a list of expressed genes, a reference
    """
    args = parse_args(argv)
    
    context_name = args.context_name
    reference_model = args.modelfile
    genefile = args.genefile
    objective = args.objective
    boundary_rxns_filepath = args.boundary_rxns_filepath
    exclude_rxns_filepath = args.exclude_rxns_filepath
    force_rxns_filepath = args.force_rxns_filepath
    recon_alg = args.algorithm.upper()
    low_threshold = args.low_threshold
    high_threshold = args.high_threshold
    solver = args.solver.upper()
    output_filetypes = stringlist_to_list(args.output_filetypes)
    
    if not os.path.exists(reference_model):
        raise FileNotFoundError(f"Reference model not found at {reference_model}")
    
    if not os.path.exists(genefile):
        raise FileNotFoundError(f"Active genes file not found at {genefile}")
    
    print(f"Active Genes: {genefile}")
    
    boundary_rxns = []
    boundary_rxns_upper: list[float] = []
    boundary_rxns_lower: list[float] = []
    
    if boundary_rxns_filepath:
        boundary_rxns_filepath: Path = Path(boundary_rxns_filepath)
        
        print(f"Boundary Reactions: {str(boundary_rxns_filepath)}")
        if boundary_rxns_filepath.suffix == ".csv":
            df: pd.DataFrame = pd.read_csv(boundary_rxns_filepath, header=0, sep=",")
        elif boundary_rxns_filepath.suffix in [".xlsx", ".xls"]:
            df: pd.DataFrame = pd.read_excel(boundary_rxns_filepath, header=0)
        else:
            raise FileNotFoundError(
                f"Boundary reactions file not found! Must be a csv or Excel file. Searching for: {boundary_rxns_filepath}")
        
        # convert all columns to lowercase
        df.columns = [column.lower() for column in df.columns]
        
        # Make sure the columns are named correctly. They should be "Reaction", "Abbreviation", "Compartment", "Minimum Reaction Rate", and "Maximum Reaction Rate"
        for column in df.columns:
            if column not in ["reaction", "abbreviation", "compartment", "minimum reaction rate",
                              "maximum reaction rate"]:
                raise ValueError(
                    f"Boundary reactions file must have columns named 'Reaction', 'Abbreviation', 'Compartment', 'Minimum Reaction Rate', and 'Maximum Reaction Rate'. Found: {column}")
        
        reaction_type: list[str] = df["reaction"].tolist()
        reaction_abbreviation: list[str] = df["abbreviation"].tolist()
        reaction_compartment: list[str] = df["compartment"].tolist()
        boundary_rxns_lower = df["minimum reaction rate"].tolist()
        boundary_rxns_upper = df["maximum reaction rate"].tolist()
        
        reaction_formula: list[str] = []
        for i in range(len(reaction_type)):
            current_type: str = reaction_type[i]
            temp_reaction: str = ""
            
            match current_type.lower():
                case "exchange":
                    temp_reaction += "EX_"
                case "demand":
                    temp_reaction += "DM_"
                case "sink":
                    temp_reaction += "SK_"
            
            shorthand_compartment = Compartments.get(reaction_compartment[i])
            temp_reaction += f"{reaction_abbreviation[i]}[{shorthand_compartment}]"
            boundary_rxns.append(temp_reaction)
            # reaction_formula.append(temp_reaction)
        
        del df
    
    exclude_rxns = []
    if exclude_rxns_filepath:
        exclude_rxns_filepath: Path = Path(exclude_rxns_filepath)
        try:
            print(f"Reading {exclude_rxns_filepath} for exclude reactions")
            if exclude_rxns_filepath.suffix == ".csv":
                df = pd.read_csv(exclude_rxns_filepath, header=0)
            elif exclude_rxns_filepath.suffix == ".xlsx" or exclude_rxns_filepath.suffix == ".xls":
                df = pd.read_excel(exclude_rxns_filepath, header=0)
            else:
                print(
                    "--exclude-reactions-file must be a path to a csv or xlsx file with one column."
                )
                print(
                    "This column should be populated with reaction ids in the reference model which will not be "
                    + "included in the context specific model, regardless of omics expression"
                )
                sys.exit()
            
            exclude_rxns = df["Abbreviation"].tolist()
        
        except FileNotFoundError:
            print(
                "--exclude-reactions-file must be a path to a csv or xlsx file with one column."
            )
            print(
                "This column should be populated with reaction ids in the reference model which will not be "
                + "included in the context specific model, regardless of omics expression"
            )
            sys.exit()
        except BaseException:
            print("exclude reactions file must be a csv or xlsx with one column.")
            print(
                "This column should be populated with reaction ids in the reference model which will not be "
                + "included in the context specific model, regardless of omics expression"
            )
            sys.exit()
    
    force_rxns = []
    if force_rxns_filepath:
        force_rxns_filepath: Path = Path(force_rxns_filepath)
        try:
            print(f"Force Reactions: {force_rxns_filepath}")
            if force_rxns_filepath.suffix == ".csv":
                df = pd.read_csv(force_rxns_filepath, header=0)
            elif force_rxns_filepath.suffix == ".xlsx" or force_rxns_filepath.suffix == ".xls":
                df = pd.read_excel(force_rxns_filepath, header=0)
            else:
                print(
                    "--force-reactions-file must be a path to a csv or xlsx file with one column."
                )
                print(
                    "This column should be populated with reaction ids in the reference model which will not be "
                    + "included in the context specific model, regardless of omics expression"
                )
                sys.exit()
            
            force_rxns = df["Abbreviation"].tolist()
        
        except FileNotFoundError:
            print(
                "--force-reactions-file must be a path to a csv or xlsx file with one column."
            )
            print(
                "This column should be populated with reaction ids in the reference model which will be force "
                + "included in the context specific model, regardless of omics expression (unless the reaction "
                + "causes the model to be infeasible)."
            )
            sys.exit()
        except BaseException:
            print("force reactions file must be a csv or xlsx with one column.")
            print(
                "This column should be populated with reaction ids in the reference model which will be force "
                + "included in the context specific model, regardless of omics expression (unless the reaction "
                + "causes the model to be infeasible)."
            )
            sys.exit()
    
    # Assert output types are valid
    for output_type in output_filetypes:
        if output_type not in ["xml", "mat", "json"]:
            print(f"Output file type {output_type} not recognized.")
            print("Output file types must be one of the following: xml, mat, json")
            sys.exit(1)
    
    if recon_alg not in ["FASTCORE", "GIMME", "IMAT"]:
        print(
            f"Algorithm {recon_alg} not supported. Please use one of: GIMME, FASTCORE, or IMAT"
        )
        sys.exit(2)
    
    if solver not in ["GUROBI", "GLPK"]:
        print(f"Solver {solver} not supported. Please use 'GLPK' or 'GUROBI'")
        sys.exit(2)
    
    print(f"Creating '{context_name}' model using '{recon_alg}' reconstruction and '{solver}' solver")
    context_model, core_list, infeas_df = create_context_specific_model(
        reference_model,
        genefile,
        recon_alg,
        objective,
        exclude_rxns,
        force_rxns,
        boundary_rxns,
        boundary_rxns_lower,
        boundary_rxns_upper,
        solver,
        context_name,
        low_threshold,
        high_threshold
    )
    
    infeas_df.to_csv(
        os.path.join(
            configs.root_dir,
            "data",
            "results",
            context_name,
            context_name + "_infeasible_rxns.csv",
        ),
        index=False,
    )
    
    if recon_alg == "FASTCORE":
        pd.DataFrame(core_list).to_csv(
            os.path.join(
                configs.root_dir,
                "data",
                "results",
                context_name,
                context_name + "_core_rxns.csv",
            ),
            index=False,
        )
    
    output_directory = os.path.join(configs.data_dir, "results", context_name)
    if "mat" in output_filetypes:
        cobra.io.save_matlab_model(
            context_model,
            os.path.join(output_directory, f"{context_name}_SpecificModel_{recon_alg}.mat"),
        )
    if "xml" in output_filetypes:
        cobra.io.write_sbml_model(
            context_model,
            os.path.join(output_directory, f"{context_name}_SpecificModel_{recon_alg}.xml"),
        )
    if "json" in output_filetypes:
        cobra.io.save_json_model(
            context_model,
            os.path.join(output_directory, f"{context_name}_SpecificModel_{recon_alg}.json"),
        )
    
    # os.path.join(configs.datadir, "results", context_name, outputfile)
    print("")
    print(f"Saved output file to {output_directory}")
    print(f"Number of Genes: {len(context_model.genes):,}")
    print(f"Number of Metabolites: {len(context_model.metabolites):,}")
    print(f"Number of Reactions: {len(context_model.reactions):,}")
    print("\nModel successfully created!")


if __name__ == "__main__":
    main(sys.argv[1:])
