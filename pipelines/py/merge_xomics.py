#!/usr/bin/python3

import os
import sys
import unidecode
import time
import pandas as pd
import numpy as np
import json
from project import configs
from microarray_gen import *
from proteomics_gen import *
from rnaseq_gen import *
from create_tissue_specific_model import split_gene_expression_data


def merge_xomics(microarray_file=None,
                 proteomics_file=None,
                 rnaseq_file=None,
                 sheet="Sheet1",
                 expression_requirement='default'):
    """
    Merges microarray, rnaseq, and/or proteomics active gene logicals from outputs of their respective "_gen.py"
    scripts.

    :param microarray_file: filename of microarray config file in /work/data/config_sheets/
    :param proteomics_file: filename of proteomics config file in /work/data/config_sheets/
    :param rnaseq_file: filename of RNA-seq config file in /work/data/config_sheets/
    :param sheet: sheet name to use, should be context, tissue, cell type, etc
    :param expression_requirement: integer, minimum number of provided sources with active gene for a it to be in model

    :return: dictionary where keys are contexts, (tissue name, etc) and values are expression tables
    """

    microarray_dict= load_microarray_tests(filename=microarray_file)
    proteomics_dict = load_proteomics_tests(filename=proteomics_file)
    rnaseq_dict = load_rnaseq_tests(filename=rnaseq_file, model_name=sheet)
    files_dict = dict()

    keys1 = proteomics_dict.keys()
    keys2 = microarray_dict.keys()
    keys3 = rnaseq_dict.keys()
    tests = set(keys1).union(set(keys2))
    tests = tests.union(set(keys3))
    tests.discard("dummy")

    merge_data = None
    for test in tests:
        if proteomics_file:
            prote_data = proteomics_dict[test].loc[:, ['expressed', 'top']]
            prote_data.rename(columns={'expressed': 'prote_exp', 'top': 'prote_top'}, inplace=True)
            test = unidecode.unidecode(test)

            try:
                if not merge_data:
                    merge_data = prote_data
            except BaseException:
                merge_data = None

        if microarray_file:
            trans_data = microarray_dict[test].loc[:, ['expressed', 'top']]
            trans_data.rename(columns={'expressed': 'trans_exp', 'top': 'trans_top'}, inplace=True)
            test = unidecode.unidecode(test)

            try:
                if not merge_data:
                    merge_data = trans_data
            except BaseException:
                merge_data = merge_data.join(trans_data, how='outer')

        if rnaseq_file:
            rnaseq_data = rnaseq_dict[test].loc[:, ['expressed', 'top']]
            rnaseq_data.rename(columns={'expressed': 'rnaseq_exp', 'top': 'rnaseq_top'}, inplace=True)
            test = unidecode.unidecode(test)

            try:
                if not merge_data:
                    merge_data = rnaseq_data
            except BaseException:
                merge_data = merge_data.join(rnaseq_data, how='outer')

        merge_data = mergeLogicalTable(merge_data)

        exp_list = []
        top_list = []
        if proteomics_file:
            exp_list.append('prote_exp')
            top_list.append('prote_top')

        if microarray_file:
            exp_list.append('trans_exp')
            top_list.append('trans_top')

        if rnaseq_file:
            exp_list.append('rnaseq_exp')
            top_list.append('rnaseq_top')
        
        num_sources = len(exp_list)
        if expression_requirement=='default':
            expression_requirement = num_sources
        
        merge_data['Express'] = 0
        merge_data['Required'] = 0
        merge_data.loc[:,'Required'] = merge_data[exp_list].apply(
            lambda x: expression_requirement-(num_sources-x.count())
            if (expression_requirement-(num_sources-x.count()) > 0)
            else 1,
            axis=1)

        merge_data.loc[merge_data[exp_list].sum(axis=1)>=merge_data['Required'], 'Express'] = 1
        merge_data.loc[merge_data[top_list].sum(axis=1) > 0, 'Express'] = 1

        filepath = os.path.join(configs.rootdir, 'data', 'results', test, 'merged_{}.csv'.format(test))
        merge_data.to_csv(filepath, index_label='ENTREZ_GENE_ID')

        filepath = os.path.join(configs.rootdir, 'data', 'results', test, 'GeneExpression_{}_Merged.csv'.format(test))
        merge_data.reset_index(drop=False, inplace=True)

        splitEntrez = split_gene_expression_data(merge_data)
        splitEntrez.rename(columns={'Gene':'ENTREZ_GENE_ID', 'Data':'Express'}, inplace=True)
        splitEntrez.to_csv(filepath, index_label='ENTREZ_GENE_ID')
        files_dict[test] = filepath

        print('{}: save to {}\n'.format(test, filepath))

    return files_dict


def handle_tissue_batch(microarray_file, rnaseq_file, proteomics_file, expression_requirement):
    """
    Handle merging of different data sources for each tissue type
    """

    if microarray_file != None:
        print(f"Microarray file is '{microarray_file}'")
        micro_config_filepath = os.path.join(configs.rootdir, "data", "config_sheets", microarray_file)
        xl = pd.ExcelFile(micro_config_filepath)
        sheet_names = xl.sheet_names

    if rnaseq_file != None:
        print(f"Bulk RNA-seq file is '{rnaseq_file}'")
        rnaseq_config_filepath = os.path.join(configs.rootdir, "data", "config_sheets", rnaseq_file)
        xl = pd.ExcelFile(rnaseq_config_filepath)
        if 'sheet_names' in locals():
            if sheet_names != xl.sheet_names:
                print("Sheet names in config files must be the same!")
                sys.exit()
        else:
            sheet_names = xl.sheet_names

    if proteomics_file != None:
        print(f"Proteomics file is '{proteomics_file}'")
        proteomics_config_filepath = os.path.join(configs.rootdir, "data", "config_sheets", proteomics_file)
        xl = pd.ExcelFile(proteomics_config_filepath)
        if 'sheet_names' in locals():
            if sheet_names != xl.sheet_names:
                print("Sheet names in config files must be the same!")
                sys.exit()
        else:
            sheet_names = xl.sheet_names

    dict_list = {}
    files_json = os.path.join(configs.rootdir, "data", "results", "step1_results_files.json")
    for tissue_name in sheet_names:
        files_dict = merge_xomics(microarray_file=microarray_file,
                                  proteomics_file=proteomics_file,
                                  rnaseq_file=rnaseq_file,
                                  expression_requirement=expression_requirement,
                                  sheet=tissue_name)

        dict_list.update(files_dict)

    with open(files_json, 'w') as fp:
        json.dump(dict_list, fp)

    return


def main(argv):
    """
    Merge expression tables of multiple sources, microarray, RNA-seq, and/or proteomics into one list

    User can specify the number of sources with an active gene in order for it to be considered active in the model.
    Otherwise, it defaults to the number of sources provided. High-confidence genes from any source will be considered
    active in the model, regardless of agreement with other sources.
    """

    parser = argparse.ArgumentParser(
        prog="merge_xomics.py",
        description="Merge expression tables of multiple sources, microarray, RNA-seq, and/or proteomics into one",
        epilog="For additional help, please post questions/issues in the MADRID GitHub repo at "
               "https://github.com/HelikarLab/MADRID or email babessell@gmail.com",
        usage="python3 $(prog)s [options]"
    )


    parser.add_argument("-m", "--microarray-config-file",
                        type=str,
                        required=False,
                        default=None,
                        dest="microarray_file",
                        help="Name of microarray config .xlsx file in the /work/data/config_files/."
                        )

    parser.add_argument("-r", "--rnaseq-config-file",
                        type=str,
                        required=False,
                        default=None,
                        dest="rnaseq_file",
                        help="Name of RNA-seq config .xlsx file in the /work/data/config_files/."
                        )

    parser.add_argument("-p", "--proteomics-config-file",
                        type=str,
                        required=False,
                        default=None,
                        dest="proteomics_file",
                        help="Name of proteomics config .xlsx file in the /work/data/config_files/."
                        )

    parser.add_argument("-e", "--expression-requirement",
                        type=str,
                        required=False,
                        default=None,
                        dest="expression_requirement",
                        help="Number of sources with active gene for it to be considered active even if it is not a "
                             "high confidence-gene"
                        )

    args = parser.parse_args(argv)

    microarray_file = args.microarray_file
    proteomics_file = args.proteomics_file
    rnaseq_file = args.rnaseq_file
    expression_requirement = args.expression_requirement

    if expression_requirement.lower() == "default":
       expression_requirement = sum([1 for test in [microarray_file, rnaseq_file, proteomics_file] if test==None])

    else:
        try:
            expression_requirement = int(expression_requirement)
        except ValueError:
            print("Expression requirement must be able to be converted to an integer!")
            sys.out()

    handle_tissue_batch(microarray_file, rnaseq_file, proteomics_file, expression_requirement)

    return


if __name__ == "__main__":
   main(sys.argv[1:])
