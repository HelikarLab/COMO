#!/usr/bin/python3
import sys
import argparse
import cobra
import pandas as pd
import numpy as np
import os
import re
import collections
from warnings import warn
from cobamp.wrappers import COBRAModelObjectReader
from troppo.methods.reconstruction.gimme import GIMME, GIMMEProperties
from troppo.methods.reconstruction.fastcore import FASTcore, FastcoreProperties
from project import configs


def correct_bracket(rule, name):
    """
    Correct GPR rules to format readable by
    """
    rmatch = re.search(r'or|and', rule)
    nmatch = re.search(r'or|and', name)
    if rmatch is None:
        lrule = rule
        lname = name.strip()
        rrule = ''
        rname = ''
        operator = ''
    else:
        lrule = rule[0:rmatch.span()[0]]
        lname = name[0:nmatch.span()[0]].strip()
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
        rule_right = ''
    else:
        rule_right = correct_bracket(rrule, rname)

    return ' '.join([rule_left, operator, rule_right])


def gene_rule_logical(expression_in, level=0):
    """
    creates an expression from GPR rule which can be evauated as true or false
    """
    try:
        loc_r = expression_in.index(')')
    except BaseException:
        if 'and' in expression_in:
            expression_in = expression_in.replace('and', ',')
            expression_in = 'min{' + expression_in + '}'
        elif 'or' in expression_in:
            expression_in = expression_in.replace('or', ',')
            expression_in = 'max{' + expression_in + '}'
        else:
            expression_in = expression_in.replace('[', '')
            expression_in = expression_in.replace(']', '')
        return expression_in
    loc_l = expression_in[0:loc_r].rindex('(')
    inner_string = expression_in[loc_l:loc_r + 1]
    inner_string = inner_string.replace('(', '[')
    inner_string = inner_string.replace(')', ']')
    if 'and' in inner_string:
        inner_string = inner_string.replace('and', ',')
        inner_string = 'min{' + inner_string + '}'
    elif 'or' in inner_string:
        inner_string = inner_string.replace('or', ',')
        inner_string = 'max{' + inner_string + '}'
    else:
        inner_string = inner_string.replace('[', '')
        inner_string = inner_string.replace(']', '')

    expression_out = '{}{}{}'.format(expression_in[0:loc_l], inner_string, expression_in[loc_r + 1:])
    expression_out = gene_rule_logical(expression_out, level + 1)

    return expression_out


def gene_rule_evaluable(expression_in):
    """
    Make expression rule evaluable
    """
    gene_reaction_by_rule = gene_rule_logical(expression_in)
    gene_reaction_by_rule = gene_reaction_by_rule.replace('{', '(')
    gene_reaction_by_rule = gene_reaction_by_rule.replace('}', ')')

    return gene_reaction_by_rule


def create_context_specific_model(general_model_file, gene_expression_file, recon_algorithm, objective,
                                 exclude_rxns, force_rxns, bound_rxns, bound_lb, bound_ub, solver):
    """
    Seed a context specific model. Core reactions are determined from GPR associations with gene expression logicals.
    Core reactions that do not necessarily meet GPR association requirements can be forced if in the force reaction
    file. Metabolite exchange (media), sinks, and demands are determined from exchanges file. Reactions can also be
    force excluded even if they meet GPR association requirements using the force exclude file.
    """
    # Check for correct format of reference model
    if general_model_file[-4:] == ".mat":
        model_cobra = cobra.io.load_matlab_model(general_model_file)  # load matlab model
    elif general_model_file[-4:] == ".xml":
        model_cobra = cobra.io.load_sbml_model(general_model_file)  # load reference model
    elif general_model_file[-5:] == ".json":
        model_cobra = cobra.io.load_json_model(general_model_file)  # load json model
    else:
        raise NameError("reference model format must be .xml, .mat, or .json")

    model_cobra.objective = {getattr(model_cobra.reactions, objective): 1}  # set objective
    all_rxns = model_cobra.reactions  # get all reactions

    # get boundary reactions
    exchange_rxns = [rxn.id for rxn in all_rxns if re.search("EX_", rxn.id)]
    sink_rxns = [rxn.id for rxn in all_rxns if re.search("SK_", rxn.id)]
    demand_rxns = [rxn.id for rxn in all_rxns if re.search("DM_", rxn.id)]

    # set flag that allows all boundary reactions to be used if none are given
    if len(bound_rxns) == 0:
        allow_all_boundary_rxns = True
    else:
        allow_all_boundary_rxns = False

    # close sinks and demands not in boundary reactions unless no boundary reactions were given
    if not allow_all_boundary_rxns:

        for i, rxn in enumerate(sink_rxns):  # set sinks to 0
            if rxn not in bound_rxns:  # exchange file
                getattr(model_cobra.reactions, rxn).lower_bound = 0.
            else:
                getattr(model_cobra.reactions, rxn).lower_bound = bound_lb[i]
                getattr(model_cobra.reactions, rxn).upper_bound = bound_ub[i]

        for i, rxn in enumerate(demand_rxns):
            if rxn not in bound_rxns:  # exchange file
                getattr(model_cobra.reactions, rxn).upper_bound = 0.
            else:
                getattr(model_cobra.reactions, rxn).lower_bound = bound_lb[i]

        # Reaction media
        medium = model_cobra.medium  # get reaction media to modify
        for rxn in exchange_rxns:  # open exchanges from exchange file, close unspecified exchanges
            if rxn not in bound_rxns:
                medium[rxn] = 0.0
            else:
                medium[rxn] = -float(bound_lb[bound_rxns.index(rxn)])

        model_cobra.medium = medium  # set new media

    # check number of unsolvable reactions for reference model under media assumptions
    model_cobra_rm = cobra.flux_analysis.fastcc(model_cobra)  # create flux consistant model (rmemoves some reactions)
    incon_rxns = set(model_cobra.reactions.list_attr("id")) - set(model_cobra_rm.reactions.list_attr("id"))
    incon_rxns_cnt = len(incon_rxns)
    del model_cobra_rm
    print(('Under given boundary assumptions, there are "{}" infeasible reactions in the ' +
           'reference model.\n').format(str(incon_rxns_cnt)))
    print('These reactions will not be considered active in context specific model construction. If any infeasible ' +
          'reactions are found to be active according to expression data, or, are found in the force reactions ' +
          'list, they can be found found in "InfeasibleRxns.csv\n')
    print('It is normal for this value to be quite large, however, if many of these reactions are active ' +
          'according to your expression data, it is likely that you are missing some critical exchange (media) ' +
          'reactions.\n')

    # CoBAMP methods
    cobamp_model = COBRAModelObjectReader(model_cobra)  # load model in readable format for CoBAMP
    cobamp_model.get_irreversibilities(True)
    s_matrix = cobamp_model.get_stoichiometric_matrix()
    lb, ub = cobamp_model.get_model_bounds(False, True)
    rx_names = cobamp_model.get_reaction_and_metabolite_ids()[0]

    # get expressed reactions
    expression_rxns, expr_vector = map_expression_to_rxn(model_cobra, gene_expression_file)
    # expressed_rxns = list({k: v for (k, v) in expression_rxns.items() if v > 0}.keys())

    # find reactions in the force reactions file that are not in general model and warn user
    for rxn in force_rxns:
        if rxn not in rx_names:
            warn(f"{rxn} from force reactions not in general model, check BiGG (or relevant database used in your "
                 "general model) for synonyms")

    # collect list of reactions that are infeasible but active in expression data or user defined
    infeas_exp_rxns = []
    infeas_force_rxns = []
    infeas_exp_cnt = 0
    infeas_force_cnt = 0
    for idx, (rxn, exp) in enumerate(expression_rxns.items()):
        if rxn in incon_rxns and expr_vector[idx] == 1:
            infeas_exp_cnt += 1
            infeas_exp_rxns.append(rxn)
        if rxn in incon_rxns and rxn in force_rxns:
            infeas_force_cnt += 1
            infeas_force_rxns.append(rxn)

        # make changes to expressed reactions base on user defined force/exclude reactions
        if rxn in force_rxns:
            expr_vector[idx] = 1  # set force reactions to 1
        if rxn in incon_rxns or rxn in exclude_rxns:
            expr_vector[idx] = 0  # set exclude and infeasible reactions to 0
        # if rxn in incon_rxns and expr_vector[idx] == 1:
        #    do_nothing="nothing"
        if rxn == objective:
            expr_vector[idx] = 1

    # get index of objective
    idx_objective = rx_names.index(objective)
    exp_idx_list = [idx for (idx, val) in enumerate(expr_vector) if val > 0]

    # seed contextual model with GIMME
    # `Becker and Palsson (2008). Context-specific metabolic networks are
    # consistent with experiments. PLoS Comput. Biol. 4, e1000082.`
    if recon_algorithm == "GIMME":
        properties = GIMMEProperties(
            exp_vector=expr_vector,  # np.array(gimme_data['0']),
            obj_frac=0.9,
            objectives=[{idx_objective: 1}],
            preprocess=True,
            flux_threshold=0.9
        )
        algorithm = GIMME(s_matrix, lb, ub, properties)
        model_seeded = algorithm.run()
        model_seeded_final = model_cobra.copy()
        r_ids = [r.id for r in model_seeded_final.reactions]
        to_remove_ids = [r_ids[r] for r in np.where(model_seeded == 0)[0]]
        model_seeded_final.remove_reactions(to_remove_ids, True)

    # seed contextual model with FastCORE
    # 'Vlassis, Pacheco, Sauter (2014). Fast reconstruction of compact
    # context-specific metbolic network models. PLoS Comput. Biol. 10,
    # e1003424.'
    elif recon_algorithm == "FASTCORE":
        properties = FastcoreProperties(core=exp_idx_list, solver=solver)
        algorithm = FASTcore(s_matrix, lb, ub, properties)
        context_rxns = algorithm.fastcore()
        model_seeded_final = model_cobra.copy()
        r_ids = [r.id for r in model_seeded_final.reactions]
        remove_rxns = [r_ids[int(i)] for i in range(s_matrix.shape[1]) if i not in context_rxns]
        model_seeded_final.remove_reactions(remove_rxns, True)

    else:
        print("Invalid reconstruction algorithm")
        return None

    # check number of unsolvable reactions for seeded model under media assumptions
    model_cobra_rm = cobra.flux_analysis.fastcc(model_seeded_final)  # create flux consistant model
    incon_rxns_cs = set(model_seeded_final.reactions.list_attr("id")) - set(model_cobra_rm.reactions.list_attr("id"))
    incon_rxns_cnt_cs = len(incon_rxns_cs)
    del model_cobra_rm
    print('Under given boundary assumptions, with infeasible reactions from the general model not considered there ' +
          'are "{}" new infeasible reactions in the context-specific model.\n'.format(str(incon_rxns_cnt_cs)))
    print('These reactions will be removed from the output model to ensure the model is solvable')
    print('Note that this value should be very low compared to the reference model.')

    incon_df = pd.DataFrame({'general_infeasible_rxns': list(incon_rxns)})
    infeas_exp_df = pd.DataFrame({'expressed_infeasible_rxns': infeas_exp_rxns})
    infeas_force_df = pd.DataFrame({'infeasible_rxns_in_force_list': infeas_exp_rxns})
    incon_df_cs = pd.DataFrame({'infeasible_rxns_from_seeding': list(incon_rxns_cs)})
    infeasible_df = pd.concat([incon_df, infeas_exp_df, infeas_force_df, incon_df_cs], ignore_index=True, axis=1)

    return model_seeded_final, exp_idx_list, infeasible_df


def split_gene_expression_data(expression_data):
    """
    Splits genes that have mapped to multiple Entrez IDs are formated as "gene12//gene2//gene3"
    """
    expression_data.rename(columns={'ENTREZ_GENE_ID': 'Gene', 'Active': 'Data'}, inplace=True)
    expression_data = expression_data.loc[:, ['Gene', 'Data']]
    expression_data['Gene'] = expression_data['Gene'].astype(str)
    single_gene_names = expression_data[~expression_data.Gene.str.contains('//')].reset_index(drop=True)
    multiple_gene_names = expression_data[expression_data.Gene.str.contains('//')].reset_index(drop=True)
    breaks_gene_names = pd.DataFrame(columns=['Gene', 'Data'])
    # print(single_gene_names.shape)
    # print(multiple_gene_names.shape)
    for index, row in multiple_gene_names.iterrows():
        for genename in row['Gene'].split('///'):
            #breaks_gene_names = breaks_gene_names.append({'Gene': genename, 'Data': row['Data']}, ignore_index=True)
            breaks_gene_names = pd.concat([breaks_gene_names, pd.DataFrame([{'Gene': genename, 'Data': row['Data']}])],
                                          axis=0, ignore_index=True)
    # print(breaks_gene_names.shape)
    gene_expressions = pd.concat([single_gene_names, breaks_gene_names], axis=0, ignore_index=True)
    # print(gene_expressions.shape)
    gene_expressions.set_index('Gene', inplace=True)

    return gene_expressions


def map_expression_to_rxn(model_cobra, gene_expression_file):
    """
    Map gene ids to a reaction based on GPR (gene to protein to reaction) association rules
    which are defined in general genome-scale metabolic model
    """
    expression_data = pd.read_csv(gene_expression_file)
    gene_expressions = split_gene_expression_data(expression_data)
    expression_rxns = collections.OrderedDict()
    # expr_vector = []
    cnt = 0
    for rxn in model_cobra.reactions:
        gene_reaction_rule = correct_bracket(rxn.gene_reaction_rule, rxn.gene_name_reaction_rule)
        gene_ids = re.findall(r'\d+', gene_reaction_rule)
        expression_rxns[rxn.id] = -1.0
        if gene_reaction_rule.strip() == '':
            continue
        for gid in gene_ids:
            boolval = '-1'
            if gid in gene_expressions.index:
                boolval = '{}'.format(gene_expressions.at[gid, 'Data'])
            gene_reaction_rule = gene_reaction_rule.replace('{}'.format(gid), boolval, 1)

        try:
            gene_reaction_by_rule = gene_rule_evaluable(gene_reaction_rule)
            gene_reaction_by_rule = gene_reaction_by_rule.strip()
            expression_rxns[rxn.id] = eval(gene_reaction_by_rule)
            # expr_vector.append(expression_rxns[rxn.id])
        except BaseException:
            # print(gene_reaction_by_rule)
            cnt += 1
    print('Map gene expression to reactions, {} errors.'.format(cnt))
    expr_vector = np.array(list(expression_rxns.values()), dtype=float)
    return expression_rxns, expr_vector


def main(argv):
    """
    Seed a context-specific model from a list of expressed genes, a reference
    """
    parser = argparse.ArgumentParser(
        prog="create_context_specific_model.py",
        description="Seed a context-specific model from a list of expressed genes, a reference",
        epilog="For additional help, please post questions/issues in the MADRID GitHub repo at "
               "https://github.com/HelikarLab/MADRID or email babessell@gmail.com",
    )

    parser.add_argument("-n", "--context-name",
                        type=str,
                        required=True,
                        dest="context_name",
                        help="Name of context or context used consistent with outputs of merge_xomics.py."
                        )

    parser.add_argument("-m", "--reference-model-filepath",
                        type=str,
                        required=True,
                        dest="modelfile",
                        help="Name of Genome-scale metabolic model to seed the context model to. For example, the "
                             "GeneralModel.mat, is a modified Recon3D model. We also provide iMM_madrid for mouse."
                             "OT can be .mat, .xml, or .json."
                        )

    parser.add_argument("-g", "--active-genes-filepath",
                        type=str,
                        required=True,
                        dest="genefile",
                        help="Path to logical table of active genes output from merge_xomics.py called "
                             "ActiveGenes_contextName_Merged.csv. Should be in the corresponding context/context folder "
                             "inside /work/data/results/contextName/. The json file output from the function using "
                             "the context of interest as the key can be used here."
                        )

    parser.add_argument("-o", "--objective",
                        type=str,
                        required=False,
                        default="biomass_reaction",
                        dest="objective",
                        help="Reaction ID of the objective function in the model. Generally a biomass function. "
                        )

    parser.add_argument("-b", "--boundary-reactions-filepath",
                        type=str,
                        required=False,
                        default=None,
                        dest="bound_rxns_file",
                        help="Path to file contains the exchange (media), sink, and demand reactions which "
                             "the model should use to fulfill the reactions governed by transcriptomic and proteomics "
                             "data inputs. It must be a csv with three columns: Rxn, Lowerbound, Upperbound. If not "
                             "specified, MADRID will allow ALL BOUNDARY REACTIONS THAT ARE OPEN IN THE REFERENCE MODEL "
                             "TO BE USED!"
                        )

    parser.add_argument("-x", "--exclude-reactions-filepath",
                        type=str,
                        required=False,
                        default=None,
                        dest="exclude_rxns_file",
                        help="Filepath to file that contains reactions which will be removed from active reactions "
                             "the model to use when seeding, even if considered active from transcriptomic and "
                             "proteomics data inputs. It must be a csv with one column of reaction IDs consistent with "
                             "the reference model"
                        )

    parser.add_argument("-f", "--force-reactions-filepath",
                        type=str,
                        required=False,
                        default=None,
                        dest="force_rxns_file",
                        help="Filepath to file that contains reactions which will be added to active reactions for "
                             "the model to use when seeding (unless it causes the model to be unsolvable), regardless "
                             "of results of transcriptomic and proteomics data inputs. It must be a csv with one "
                             "column of reaction IDs consistent with the reference model"
                        )

    parser.add_argument("-a", "--algorithm",
                        type=str,
                        required=False,
                        default="GIMME",
                        dest="algorithm",
                        help="Algorithm used to seed context specific model to the Genome-scale model. Can be either "
                             "GIMME or FastCORE."
                        )

    parser.add_argument("-s", "--solver",
                        type=str,
                        required=False,
                        default="glpk",
                        dest="solver",
                        help="Solver used to seed model and attempt to solve objective. Default is GLPK, also takes "
                             "GUROBI but you must mount a container license to the Docker to use. An academic license "
                             "can be obtained for free. See the README on the Github or Dockerhub for information on "
                             "mounting this license."
                        )
    parser.add_argument("-t", "--output-filetypes",
                        type=str,
                        required=False,
                        default="mat",
                        dest="output_filetypes",
                        help="Filetypes to save seeded model type. Can be either a string with one filetype such as "
                             "'xml' or multiple in the format \"['extension1', 'extention2', ... etc]\". If you want "
                             "to output in all 3 accepted formats,  would be: \"['mat', 'xml', 'json']\" "
                             "Note the outer quotes required to be interpreted by cmd. This a string, not a python list"
                        )

    args = parser.parse_args(argv)

    context_name = args.context_name
    reference_model = args.modelfile
    genefile = args.genefile
    objective = args.objective
    bound_rxns_file = args.bound_rxns_file
    exclude_rxns_file = args.exclude_rxns_file
    force_rxns_file = args.exclude_rxns_file
    recon_alg = args.algorithm.upper()
    solver = args.solver.upper()
    output_filetypes = args.output_filetypes

    if not os.path.exists(reference_model):
        raise FileNotFoundError(f"Reference model not found at {reference_model}")

    if not os.path.exists(genefile):
        raise FileNotFoundError(f"Active genes file not found at {genefile}")


    if bound_rxns_file:
        try:
            print(f"Reading {bound_rxns_file} for boundary reactions")
            df = pd.read_csv(bound_rxns_file, header=0, sep=",")
            bound_rxns = df['Rxn'].tolist()
            bound_ub = df['Upperbound'].tolist()
            bound_lb = df['Lowerbound'].tolist()
            del df
        except FileNotFoundError:
            print("--boundary-reactions-file must be a path to a csv file with three columns: "
                  "Rxn, Lowerbound, Upperbound")
        except BaseException:
            print("Boundary reactions file must be a csv with three columns: Rxn, Lowerbound, Upperbound")
    else:
        bound_rxns = []
        bound_ub = []
        bound_lb = []

    if exclude_rxns_file:
        try:
            print(f"Reading {exclude_rxns_file} for boundary reactions")
            df = pd.read_csv(exclude_rxns_file, header=None).values.tolist()
            exclude_rxns = df.values.tolist()
        except FileNotFoundError:
            print("--exclude-reactions-file must be a path to a csv file with one column.")
            print("This column should be populated with reaction ids in the reference model which will not be " +
                  "included in the context specific model, regardless of omics expression")
        except BaseException:
            print("exclude reactions file must be a csv with one column.")
            print("This column should be populated with reaction ids in the reference model which will not be " +
                  "included in the context specific model, regardless of omics expression")
    else:
        exclude_rxns = []

    if force_rxns_file:
        try:
            print(f"Reading {force_rxns_file} for boundary reactions")
            df = pd.read_csv(force_rxns_file, header=None).values.tolist()
            force_rxns = df.values.tolist()
        except FileNotFoundError:
            print("--force-reactions-file must be a path to a csv file with one column.")
            print("This column should be populated with reaction ids in the reference model which will be force " +
                  "included in the context specific model, regardless of omics expression (unless the reaction " +
                  "causes the model to be infeasible).")
        except BaseException:
            print("force reactions file must be a csv with one column.")
            print("This column should be populated with reaction ids in the reference model which will be force " +
                  "included in the context specific model, regardless of omics expression (unless the reaction " +
                  "causes the model to be infeasible).")
    else:
        force_rxns = []

    def print_filetype_help():
        print("Unsupported model format. Supports 'xml', 'mat', and 'json''.")
        print("Or use multiple with: \"['.extension1', 'extention2', ... etc]\"")
        print("For example, if you want to output in all 3 accepted formats:")
        print("\"['mat', 'xml', 'json']\"")
        print("Note the outer quotes required to be interpreted by cmd. This a string, not a python list")

    if output_filetypes not in ["xml", "mat", "json"]:
        try:
            output_filetypes = output_filetypes.strip("[").strip("]").replace("'", "").replace(" ", "").split(",")
            if any(form not in ["xml", "json", "mat"] for form in output_filetypes):
                print_filetype_help()
                sys.exit()
        except:
            print_filetype_help()
    else:
        output_filetypes = [output_filetypes]

    if recon_alg not in ["FASTCORE", "GIMME"]:
        print(f"Algorithim {recon_alg} not supported. Please use 'GIMME' or 'FASTCORE'")
        sys.exit(2)

    if solver not in ["GUROBI", "GLPK"]:
        print(f"Solver {solver} not supported. Please use 'GLPK' or 'GUROBI'")
        sys.exit(2)

    print(f"Constructing model of {context_name} with {recon_alg} reconstruction algorithm using {solver} solver")
    outfile_basename = context_name + "_SpecificModel"

    context_model, core_list, infeas_df = create_context_specific_model(reference_model, genefile, recon_alg,
                                                                        objective, exclude_rxns, force_rxns,
                                                                        bound_rxns, bound_lb, bound_ub, solver)

    infeas_df.to_csv(os.path.join(configs.rootdir, "data", "results", context_name,
                                  context_name + "_infeasible_rxns.csv"),
                     index=False)

    if recon_alg == "FASTCORE":
        pd.DataFrame(core_list).to_csv(os.path.join(configs.rootdir, "data", "results", context_name,
                                                    context_name + "_core_rxns.csv"),
                                       index=False)
    if "mat" in output_filetypes:
        outputfile = outfile_basename + ".mat"
        print(f"Output file is '{outputfile}'")
        # cobra.io.mat.save_matlab_model(context_model, os.path.join(configs.rootdir, 'data', outputfile))
        cobra.io.save_matlab_model(context_model, os.path.join(configs.datadir, "results", context_name, outputfile))
    if "xml" in output_filetypes:
        outputfile = outfile_basename + ".xml"
        print(f"Output file is '{outputfile}'")
        print("**Note cobrapy only supports level 2 SBML model, while this model is level 3")
        cobra.io.write_sbml_model(context_model, os.path.join(configs.datadir, "results", context_name, outputfile))
    if "json" in output_filetypes:
        outputfile = outfile_basename + ".json"
        print(f"Output file is '{outputfile}'")
        cobra.io.save_json_model(context_model, os.path.join(configs.datadir, "results", context_name, outputfile))

    print("Number of Genes: " + str(len(context_model.genes)))
    print("Number of Metabolites: " + str(len(context_model.metabolites)))
    print("Number of Reactions: " + str(len(context_model.reactions)))
    print(context_model.objective._get_expression())
    print(context_model.optimize())
    return None


if __name__ == "__main__":
    main(sys.argv[1:])
