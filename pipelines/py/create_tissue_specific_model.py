#!/usr/bin/python3
import sys
import getopt
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
    Make expression rule evaluatable
    """
    gene_reaction_by_rule = gene_rule_logical(expression_in)
    gene_reaction_by_rule = gene_reaction_by_rule.replace('{', '(')
    gene_reaction_by_rule = gene_reaction_by_rule.replace('}', ')')

    return gene_reaction_by_rule


def create_tissue_specific_model(general_model_file, gene_expression_file, recon_algorithm, objective,
                                 exclude_rxns, force_rxns, bound_rxns, bound_lb, bound_ub, solver):
    """
    Seed a context specific model. Core reactions are determined from GPR associations with gene expression logicals.
    Core reactions that do not necessarily meet GPR association requirements can be forced if in the force reaction
    file. Metabolite exchange (media), sinks, and demands are determined from exchanges file. Reactions can also be
    force excluded even if they meet GPR association requirements using the force exclude file.
    """
    model_cobra = cobra.io.load_matlab_model(general_model_file)  # load reference model
    model_cobra.objective = {getattr(model_cobra.reactions, objective): 1}  # set objective
    all_rxns = model_cobra.reactions  # get all reactions

    # get boundary reactions
    exchange_rxns = [rxn.id for rxn in all_rxns if re.search("EX_", rxn.id)]
    sink_rxns = [rxn.id for rxn in all_rxns if re.search("SK_", rxn.id)]
    demand_rxns = [rxn.id for rxn in all_rxns if re.search("DM_", rxn.id)]

    # close sinks and demands
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
    del model_cobra_rm
    print(f''"Under given boundary assumptions, there are {incon_rxns_cnt} infeasible reactions in the\n \
            reference model.\n\n \
            These reactions will not be considered active in context specific model contstruction. If\n \
            any infeasible reactions are found to be active according to expression data, or, are\n \
            found in the force reactions list, they will be specified in 'InfeasibleRxns.csv'\n\n \
            \n\nNote that it is normal for this value to be very large, however, if many of these\n \
            reactions are active according to your expression data, it is likely that you are\n \
            missing some critical exchange (media) reactions.")

    # CoBAMP methods
    cobamp_model = COBRAModelObjectReader(model_cobra)  # load model in readable format for CoBAMP
    cobamp_model.get_irreversibilities(True)
    s_matrix = cobamp_model.get_stoichiometric_matrix()
    lb, ub = cobamp_model.get_model_bounds(False, True)
    rx_names = cobamp_model.get_reaction_and_metabolite_ids()[0]
    # pd.DataFrame(rx_names).to_csv(os.path.join(configs.rootdir, 'data', 'results', 'liver_control', 'rx_names'))

    # get expressed reactions
    expression_rxns, expr_vector = map_expression_to_rxn(model_cobra, gene_expression_file)
    # expressed_rxns = list({k: v for (k, v) in expression_rxns.items() if v > 0}.keys())

    # find reactions in the force reactions file that are not in general model
    for rxn in force_rxns:
        if rxn not in rx_names:
            warn(f"{rxn} from force reactions not in general model, check BiGG \
            (or database for reactions used in your general model) for synonyms")

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
        tissue_rxns = algorithm.fastcore()
        model_seeded_final = model_cobra.copy()
        r_ids = [r.id for r in model_seeded_final.reactions]
        remove_rxns = [r_ids[int(i)] for i in range(s_matrix.shape[1]) if i not in tissue_rxns]
        model_seeded_final.remove_reactions(remove_rxns, True)

    else:
        print("Invalid reconstruction algorithm")
        return None

    # check number of unsolvable reactions for seeded model under media assumptions
    model_cobra_rm = cobra.flux_analysis.fastcc(model_seeded_final)  # create flux consistant model
    incon_rxns_cs = set(model_seeded_final.reactions.list_attr("id")) - set(model_cobra_rm.reactions.list_attr("id"))
    incon_rxns_cnt_cs = len(incon_rxns_cs)
    del model_cobra_rm
    print(f"""
    Under given boundary assumptions, with infeasible reactions from the general model not considered
    there are {incon_rxns_cnt_cs} new infeasible reactions in the tissue-specific model\n\n
    These reactions will be removed from the output model to ensure they can be simulated.\n\n
    Note that this value should be very low.""")

    incon_df = pd.DataFrame({'general_infeasible_rxns': list(incon_rxns)})
    infeas_exp_df = pd.DataFrame({'expressed_infeasible_rxns': infeas_exp_rxns})
    infeas_force_df = pd.DataFrame({'infeasible_rxns_in_force_list': infeas_exp_rxns})
    incon_df_cs = pd.DataFrame({'infeasible_rxns_from_seeding': list(incon_rxns_cs)})
    infeasible_df = pd.concat([incon_df, infeas_exp_df, infeas_force_df, incon_df_cs], ignore_index=True, axis=1)
    infeasible_df.head()

    return model_seeded_final, exp_idx_list, infeasible_df


def split_gene_expression_data(expression_data):
    """
    Splits genes that have mapped to multiple Entrez IDs are formated as "gene12//gene2//gene3"
    """
    expression_data.rename(columns={'ENTREZ_GENE_ID': 'Gene', 'Express': 'Data'}, inplace=True)
    expression_data = expression_data.loc[:, ['Gene', 'Data']]
    expression_data['Gene'] = expression_data['Gene'].astype(str)
    single_gene_names = expression_data[~expression_data.Gene.str.contains('//')].reset_index(drop=True)
    multiple_gene_names = expression_data[expression_data.Gene.str.contains('//')].reset_index(drop=True)
    breaks_gene_names = pd.DataFrame(columns=['Gene', 'Data'])
    # print(single_gene_names.shape)
    # print(multiple_gene_names.shape)
    for index, row in multiple_gene_names.iterrows():
        for genename in row['Gene'].split('///'):
            breaks_gene_names = breaks_gene_names.append({'Gene': genename, 'Data': row['Data']}, ignore_index=True)
    # print(breaks_gene_names.shape)
    gene_expressions = single_gene_names.append(breaks_gene_names, ignore_index=True)
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
    # default arguments, all other args are required for script to run
    objective = 'biomass_reaction'
    exclude_rxns = []
    force_rxns = []
    bound_rxns = []
    bound_ub = []
    bound_lb = []
    recon_alg = "GIMME"
    solver = "GLPK"
    out_formats = [".mat"]

    try:
        opts, args = getopt.getopt(argv, "hn:m:g:o:b:x:f:a:s:t:",
                                   ["tissue_name=", "reference_model_file=", "gene_expression_file=",
                                    "objective=", "boundary_reactions_file", "exclude_reactions=",
                                    "forceRxns=", "algorithm=", "solver=", "output_model_file_types"])
    except getopt.GetoptError:
        print('python3 create_tissue_specific_model.py -n <tissue_name*> -m <reference_model_file*> \
               -g <gene_expression_file*> -o <objective*> -b <boundary_reactions_file> \
               -x <exclude_reactions_file> -f <force_reactions_file> -a <algorithm> -s <solver> \
               -t <output_model_file_types')
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print('python3 create_tissue_specific_model.py -n <tissue_name*> -m <reference_model_file*> \
                   -g <gene_expression_file*> -o <objective*> -b <boundary_reactions_file> \
                   -x <exclude_reactions_file> -f <force_reactions_file> -a <algorithm> -s <solver> \
                   -t <output_model_file_types')
            sys.exit()
        elif opt in ("-n", "--tissue_name"):
            tissue_name = arg
        elif opt in ("-m", "--reference_model_file"):
            modelfile = arg
            # output model format defaults to input model format unless otherwise specified
            if modelfile[-4:] == ".mat":
                out_formats = [".mat"]
            elif modelfile[-4] == ".xml":
                out_formats = [".xml"]
            elif modelfile[-5] == ".json":
                out_formats = [".json"]
        elif opt in ("-g", "--gene_expression_file"):
            genefile = arg
            # assign output model file
        elif opt in ("-o", "--objective"):
            objective = arg
        elif opt in ("-b", "--boundary_reactions_file"):
            try:
                df = pd.read_csv(arg, header=0, sep=",")
                bound_rxns = df['Rxn'].to_list()
                bound_ub = df['Upperbound'].to_list()
                bound_lb = df['Lowerbound'].to_list()
                del df
            except BaseException:
                raise NameError("Force reactions must have three columns: Rxn Name, Lowerbound, Upperbound")
        elif opt in ("-x", "--exclude_reactions_file"):
            exclude_rxns = pd.read_csv(arg, header=0)['Rxn'].to_list()
        elif opt in ("-f", "--force_reactions_file"):
            try:
                df = pd.read_csv(arg, header=0, sep=",")
                force_rxns = df['Rxn'].to_list()
                del df
            except:
                raise NameError("Force reactions must have three columns: Rxn Name, Lowerbound, Upperbound")
        elif opt in ("-a", "--algorithm"):
            recon_alg = arg.upper()
        elif opt in ("-t", "--output_model_file_types"):
            out_formats = arg.strip("[").strip("]").replace("'", "").replace(" ", "").split(",")
    try:
        print('Tissue Name is "{}"'.format(tissue_name))
    except:
        raise NameError("Missing tissue name argument, use -h for list of inputs")
    try:
        print('General Model file is "{}"'.format(modelfile))
    except:
        raise NameError("Missing general model file, use -h for list of inputs")
    try:
        print('Gene Expression file is "{}"'.format(genefile))
    except:
        raise NameError("Missing gene expression file, use -h for ")

    print('Constructing model with "{}" reconstruction algorithm using "{}" solver'.format(recon_alg, solver))
    outfile_basename = tissue_name + "_SpecificModel"

    # check if any unsupported filetypes are specified before creating CS models
    if any(form not in [".xml", ".json", ".mat"] for form in out_formats):
        raise NameError('Unsupported model format. Supports ".xml", ".mat", and ".json".')

    # check if any unsupported algorithms are specified before creating CS models
    if recon_alg not in ["GIMME", "FASTCORE"]:
        raise NameError('Unsupported algorithm. Supports "GIMME", "FastCORE"')

    general_model_file = os.path.join(configs.rootdir, 'data', modelfile)
    gene_expression_file = os.path.join(configs.rootdir, 'data', 'results', tissue_name, genefile)
    tissue_model, core_list, infeas_df = create_tissue_specific_model(general_model_file, gene_expression_file,
                                                                      recon_alg, objective, exclude_rxns,
                                                                      force_rxns, bound_rxns, bound_lb, bound_ub,
                                                                      solver)

    infeas_df.to_csv(os.path.join(configs.rootdir, "data", "results", tissue_name,
                                  tissue_name + "_infeasible_rxns.csv"),
                     index=False)

    if recon_alg == "FASTCORE":
        pd.DataFrame(core_list).to_csv(os.path.join(configs.rootdir, "data", "results", tissue_name,
                                                    tissue_name + "_core_rxns.csv"),
                                       index=False)
    if ".mat" in out_formats:
        outputfile = outfile_basename + ".mat"
        print('Output file is "{}"'.format(outputfile))
        # cobra.io.mat.save_matlab_model(tissue_model, os.path.join(configs.rootdir, 'data', outputfile))
        cobra.io.save_matlab_model(tissue_model, os.path.join(configs.rootdir, "data",
                                                              "results", tissue_name, outputfile))
    if ".xml" in out_formats:
        outputfile = outfile_basename + ".xml"
        print('Output file is "{}"'.format(outputfile))
        print('**Note cobrapy only supports level 2 SBML model, while this model is level 3')
        cobra.io.write_sbml_model(tissue_model, os.path.join(configs.rootdir, 'data',
                                                             "results", tissue_name, outputfile))
    if ".json" in out_formats:
        outputfile = outfile_basename + ".json"
        print('Output file is "{}"'.format(outputfile))
        cobra.io.save_json_model(tissue_model, os.path.join(configs.rootdir, 'data',
                                                            "results", tissue_name, outputfile))

    print("Number of Genes: " + str(len(tissue_model.genes)))
    print("Number of Metabolites: " + str(len(tissue_model.metabolites)))
    print("Number of Reactions: " + str(len(tissue_model.reactions)))
    print(tissue_model.objective._get_expression())
    print(tissue_model.optimize())
    return None


if __name__ == "__main__":
    main(sys.argv[1:])
