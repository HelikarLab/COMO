#!/usr/bin/python3

import os
import sys
import unidecode
import time
import pandas as pd
import numpy as np
import json
from collections import Counter
from project import configs
from microarray_gen import *
from proteomics_gen import *
from rnaseq_gen import *
from create_tissue_specific_model import split_gene_expression_data


def merge_xomics(microarray_file=None,
                 proteomics_file=None,
                 trnaseq_file=None,
                 mrnaseq_file=None,
                 scrnaseq_file=None,
                 sheet="Sheet1",
                 expression_requirement='default'):
    """
    Merges microarray, rnaseq, and/or proteomics active gene logicals from outputs of their respective "_gen.py"
    scripts.

    :param microarray_file: filename of microarray config file in /work/data/config_sheets/
    :param proteomics_file: filename of proteomics config file in /work/data/config_sheets/
    :param trnaseq_file: filename of Total RNA-seq config file in /work/data/config_sheets/
    :param mrnaseq_file: filename of mRNA-seq config file in /work/data/config_sheets/
    :param scrnaseq_file: filename of single-cell RNA-seq config file in /work/data/config_sheets/
    :param sheet: sheet name to use, should be context, tissue, cell type, etc
    :param expression_requirement: integer, minimum number of provided sources with active gene for a it to be in model

    :return: dictionary where keys are contexts, (tissue name, etc) and values are expression tables
    """

    microarray_dict = load_microarray_tests(filename=microarray_file)
    proteomics_dict = load_proteomics_tests(filename=proteomics_file)
    trnaseq_dict = load_rnaseq_tests(filename=trnaseq_file, model_name=sheet)
    mrnaseq_dict = load_rnaseq_tests(filename=mrnaseq_file, model_name=sheet)
    scrnaseq_dict = load_rnaseq_tests(filename=scrnaseq_file, model_name=sheet)
    files_dict = dict()

    keys1 = proteomics_dict.keys()
    keys2 = microarray_dict.keys()
    keys3 = trnaseq_dict.keys()
    keys4 = mrnaseq_dict.keys()
    keys5 = scrnaseq_dict.keys()

    tests = set(keys1).union(set(keys2))
    tests = tests.union(set(keys3))
    tests = tests.union(set(keys4))
    tests = tests.union(set(keys5))

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

        if trnaseq_file:
            trnaseq_data = trnaseq_dict[test].loc[:, ['expressed', 'top']]
            trnaseq_data.rename(columns={'expressed': 'trnaseq_exp', 'top': 'trnaseq_top'}, inplace=True)
            test = unidecode.unidecode(test)

            try:
                if not merge_data:
                    merge_data = trnaseq_data
            except BaseException:
                merge_data = merge_data.join(trnaseq_data, how='outer')

        if mrnaseq_file:
            mrnaseq_data = mrnaseq_dict[test].loc[:, ['expressed', 'top']]
            mrnaseq_data.rename(columns={'expressed': 'mrnaseq_exp', 'top': 'mrnaseq_top'}, inplace=True)
            test = unidecode.unidecode(test)

            try:
                if not merge_data:
                    merge_data = mrnaseq_data
            except BaseException:
                merge_data = merge_data.join(mrnaseq_data, how='outer')

        if scrnaseq_file:
            scrnaseq_data = scrnaseq_dict[test].loc[:, ['expressed', 'top']]
            scrnaseq_data.rename(columns={'expressed': 'scrnaseq_exp', 'top': 'scrnaseq_top'}, inplace=True)
            test = unidecode.unidecode(test)

            try:
                if not merge_data:
                    merge_data = scrnaseq_data
            except BaseException:
                merge_data = merge_data.join(scrnaseq_data, how='outer')

        merge_data = mergeLogicalTable(merge_data)

        exp_list = []
        top_list = []
        if proteomics_file:
            exp_list.append('prote_exp')
            top_list.append('prote_top')

        if microarray_file:
            exp_list.append('trans_exp')
            top_list.append('trans_top')

        if trnaseq_file:
            exp_list.append('trnaseq_exp')
            top_list.append('trnaseq_top')

        if mrnaseq_file:
            exp_list.append('mrnaseq_exp')
            top_list.append('mrnaseq_top')

        if scrnaseq_file:
            exp_list.append('scrnaseq_exp')
            top_list.append('scrnaseq_top')
        
        num_sources = len(exp_list)
        if expression_requirement == 'default':
            expression_requirement = num_sources
        
        merge_data['Express'] = 0
        merge_data['Required'] = 0
        merge_data.loc[:, 'Required'] = merge_data[exp_list].apply(
            lambda x: expression_requirement-(num_sources-x.count())
            if (expression_requirement-(num_sources-x.count()) > 0)
            else 1,
            axis=1)

        merge_data.loc[merge_data[exp_list].sum(axis=1) >= merge_data['Required'], 'Express'] = 1
        merge_data.loc[merge_data[top_list].sum(axis=1) > 0, 'Express'] = 1

        filepath = os.path.join(configs.rootdir, 'data', 'results', test, 'merged_{}.csv'.format(test))
        merge_data.to_csv(filepath, index_label='ENTREZ_GENE_ID')

        filepath = os.path.join(configs.rootdir, 'data', 'results', test, 'GeneExpression_{}_Merged.csv'.format(test))
        merge_data.reset_index(drop=False, inplace=True)

        split_entrez = split_gene_expression_data(merge_data)
        split_entrez.rename(columns={'Gene': 'ENTREZ_GENE_ID', 'Data':'Express'}, inplace=True)
        split_entrez.to_csv(filepath, index_label='ENTREZ_GENE_ID')
        files_dict[test] = filepath

        print('{}: save to {}\n'.format(test, filepath))

    return files_dict


def handle_tissue_batch(microarray_file, trnaseq_file, mrnaseq_file,
                        scrnaseq_file, proteomics_file, expression_requirement):
    """
    Handle merging of different data sources for each tissue type
    """

    sheet_names = []
    for file in [microarray_file, trnaseq_file, mrnaseq_file, scrnaseq_file, proteomics_file]:
        config_filepath = os.path.join(configs.rootdir, "data", "config_sheets", file)
        xl = pd.ExcelFile(config_filepath)
        sheet_names += xl.sheet_names

    counts = Counter(sheet_names)
    sheet_names = list(set(sheet_names))
    dict_list = {}
    max_inputs = max(counts.values())
    diff_inputs = max_inputs - expression_requirement # get difference between expressoin requirement and max input
    files_json = os.path.join(configs.rootdir, "data", "results", "step1_results_files.json")
    for tissue_name in sheet_names:
        exp_req = counts[tissue_name] - diff_inputs  # correct expr requirement based on how many fewer inputs from max
        files_dict = merge_xomics(microarray_file=microarray_file,
                                  proteomics_file=proteomics_file,
                                  trnaseq_file=trnaseq_file,
                                  mrnaseq_file=mrnaseq_file,
                                  scrnaseq_file=scrnaseq_file,
                                  expression_requirement=exp_req,
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


    parser.add_argument("-a", "--microarray-config-file",
                        type=str,
                        required=False,
                        default=None,
                        dest="microarray_file",
                        help="Name of microarray config .xlsx file in the /work/data/config_files/."
                        )

    parser.add_argument("-t", "--total-rnaseq-config-file",
                        type=str,
                        required=False,
                        default=None,
                        dest="trnaseq_file",
                        help="Name of total RNA-seq config .xlsx file in the /work/data/config_files/."
                        )

    parser.add_argument("-m", "--mrnaseq-config-file",
                        type=str,
                        required=False,
                        default=None,
                        dest="mrnaseq_file",
                        help="Name of mRNA-seq config .xlsx file in the /work/data/config_files/."
                        )

    parser.add_argument("-s", "--scrnaseq-config-file",
                        type=str,
                        required=False,
                        default=None,
                        dest="scrnaseq_file",
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
    trnaseq_file = args.trnaseq_file
    mrnaseq_file = args.mrnaseq_file
    scrnaseq_file = args.scrnaseq_file
    expression_requirement = args.expression_requirement

    def_exp_req = sum([
        1 for test in [microarray_file, trnaseq_file, mrnaseq_file, scrnaseq_file, proteomics_file] if test is None])

    if expression_requirement.lower() == "default":
        expression_requirement = def_exp_req

    else:
        try:
            expression_requirement = int(expression_requirement)
        except ValueError:
            print("Expression requirement must be able to be converted to an integer!")
            sys.out()

    handle_tissue_batch(microarray_file, trnaseq_file, mrnaseq_file, scrnaseq_file, expression_requirement)

    return


if __name__ == "__main__":
   main(sys.argv[1:])
