#!/usr/bin/python3
import os
import sys
import getopt
import unidecode
import time
import pandas as pd
import numpy as np
import json
# from bioservices import BioDBNet
from project import configs
from microarray_gen import *
from proteomics_gen import *
from bulk_gen import *
from create_tissue_specific_model import splitGeneExpressionData

# input parameters
def merge_xomics(transcript_file=None,
                 prote_file=None,
                 bulk_file=None,
                 exp_req='default'):

    microarray_dict= load_microarray_tests(filename=transcript_file)
    #Proteomics = load_prot_supplementary_data(prote_file)
    proteomics_dict = load_proteomics_tests(filename=prote_file)
    #Bulk = load_bulk_supplementary_data(bulk_file)
    bulk_dict = load_bulk_tests(filename=bulk_file)
    files_dict = dict()

    keys1 = proteomics_dict.keys()
    keys2 = microarray_dict.keys()
    keys3 = bulk_dict.keys()
    tests = set(keys1).union(set(keys2))
    tests = tests.union(set(keys3))
    tests.discard("dummy")

    merge_data = None
    for test in tests:
        if prote_file:
            prote_data = proteomics_dict[test].loc[:, ['expressed', 'top']]
            prote_data.rename(columns={'expressed': 'prote_exp',
                                       'top': 'prote_top'}, inplace=True)
            #prote_data = mergeLogicalTable(prote_data)
            test = unidecode.unidecode(test)
            try:
                if not merge_data:
                    merge_data = prote_data
            except:
                merge_data = None

        if transcript_file:
            trans_data = microarray_dict[test].loc[:, ['expressed', 'top']]
            trans_data.rename(columns={'expressed': 'trans_exp',
                                       'top': 'trans_top'}, inplace=True)

            test = unidecode.unidecode(test)
            try:
                if not merge_data:
                    merge_data = trans_data
            except:
                merge_data = merge_data.join(trans_data, how='outer')

        if bulk_file:
            bulk_data = bulk_dict[test].loc[:, ['expressed', 'top']]
            bulk_data.rename(columns={'expressed': 'bulk_exp',
                                      'top': 'bulk_top'}, inplace=True)

            test = unidecode.unidecode(test)
            try:
                if not merge_data:
                    merge_data = bulk_data
            except:
                merge_data = merge_data.join(bulk_data, how='outer')

        merge_data = mergeLogicalTable(merge_data)
        #test = unidecode.unidecode(test)
        exp_list = []
        top_list = []
        if prote_file:
            exp_list.append('prote_exp')
            top_list.append('prote_top')
        if transcript_file:
            exp_list.append('trans_exp')
            top_list.append('trans_top')
        if bulk_file:
            exp_list.append('bulk_exp')
            top_list.append('bulk_top')
        
        num_sources = len(exp_list)
        if exp_req=='default': exp_req = num_sources
        
        merge_data['Express'] = 0
        merge_data['Required'] = 0
        merge_data.loc[:,'Required'] = merge_data[exp_list].apply(
                                    lambda x: exp_req-(num_sources-x.count()) \
                                           if (exp_req-(num_sources-x.count()) > 0) \
                                           else 1, \
                                    axis=1)
                                    
        #merge_data.loc[:,'Required'] = exp_req - merge_data[exp_list].isnull().sum(axis=0)
        merge_data.loc[merge_data[exp_list].sum(axis=1)>=merge_data['Required'], 'Express'] = 1
        merge_data.loc[merge_data[top_list].sum(axis=1) > 0, 'Express'] = 1
        #merge_data = merge_data['Express'].astype(int)

        filepath = os.path.join(configs.rootdir, 'data', 'results', test, 'merged_{}.csv'.format(test))
        merge_data.to_csv(filepath, index_label='ENTREZ_GENE_ID')

        filepath = os.path.join(configs.rootdir, 'data', 'results', test, 'GeneExpression_{}_Merged.csv'.format(test))
        merge_data.reset_index(drop=False, inplace=True)

        splitEntrez = splitGeneExpressionData(merge_data)
        splitEntrez.rename(columns={'Gene':'ENTREZ_GENE_ID', 'Data':'Express'}, inplace=True)
        splitEntrez.to_csv(filepath, index_label='ENTREZ_GENE_ID')
        files_dict[test] = filepath

        print('{}: save to {}\n'.format(test, filepath))

    return files_dict

def main(argv):
    transfile = None
    protefile = None
    bulkfile = None
    expression_requirement = 'default'
    try:
        opts, args = getopt.getopt(argv, "ht:p:b:r:", ["transfile=", "protefile=", "bulkfile=", "expression_requirement="])
    except getopt.GetoptError:
        print('merge_xomics.py -t <transfile> -p <protefile> -r <expression_requirement>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('merge_xomics.py -t <transfile> -p <protefile> -r <expression_requirement>')
            sys.exit()
        elif opt in ("-t", "--transfile"):
            transfile = arg
        elif opt in ("-p", "--protefile"):
            protefile = arg
        elif opt in ('-b', "--bulkfile"):
            bulkfile = arg
        elif opt in ('-r', "--expression_requirement"):
            expression_requirement = int(arg)
            
    print('Microarray file is "{}"'.format(transfile))
    print('Proteomics file is "{}"'.format(protefile))
    print('Bulk RNA-seq file is "{}"'.format(bulkfile))
    files_dict = merge_xomics(transcript_file=transfile,
                              prote_file=protefile,
                              bulk_file=bulkfile,
                              exp_req=expression_requirement)
    
    files_json = os.path.join(configs.rootdir, 'data', 'results', 'step1_results_files.json')
    with open(files_json, 'w') as fp:
        json.dump(files_dict, fp)


if __name__ == "__main__":
   main(sys.argv[1:])
