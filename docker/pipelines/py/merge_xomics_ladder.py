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
from transcriptomic_gen import *
from proteomics_gen import *
from bulk_gen import *
from create_tissue_specific_model import splitGeneExpressionData

# input parameters
def merge_xomics(transcript_file=None,
                 prote_file=None,
                 bulk_file=None):

    transcriptomics_dict= load_transcriptomics_tests(filename = transcript_file)
    Proteomics = load_prot_supplementary_data(prote_file)
    proteomics_dict = load_proteomics_tests(Proteomics)
    print(proteomics_dict)
    Bulk = load_bulk_supplementary_data(bulk_file)
    bulk_dict = load_bulk_tests(Bulk)
    print(bulk_dict)
    files_dict = dict()

    #keys1 = proteomics_dict.keys()
    #keys2 = transcriptomics_dict.keys()
    #keys3 = bulk_dict.keys()

    #tests = set(keys2).intersection(set(keys3))

    if prote_file and transcript_file and bulk_file:
        keys1 = proteomics_dict.keys()
        keys2 = transcriptomics_dict.keys()
        keys3 = bulk_dict.keys()
        tests = set(keys1).intersection(set(keys2))
        tests = tests.intersection(set(keys3))
        tests2_1 = (set(keys1).union(set(keys2))).difference(tests)
        tests2_2 = (set(keys1).union(set(keys3))).difference(tests)
        tests2_3 = (set(keys2).union(set(keys3))).difference(tests)
        tests2 = test2_1.union(tests2_2)
        tests2 = tests2.union(tests2_3)
    elif prote_file and transcript_file:
        keys1 = proteomics_dict.keys()
        keys2 = transcriptomics_dict.keys()
        tests = set(keys1).intersection(set(keys2))
        tests2 = (set(keys1).union(set(keys2))).difference(tests)
    elif prote_file and bulk_file:
        keys1 = proteomics_dict.keys()
        print(keys1)
        keys3 = bulk_dict.keys()
        print(keys3)
        tests = set(keys1).intersection(set(keys3))
        print(tests)
        tests2 = (set(keys1).union(set(keys3))).difference(tests)
        print(tests2)
    elif prote_file:
        keys1 = proteomics_dict.keys()
        tests = set(keys1)
        tests2 = tests
    elif transcript_file and bulk_file:
        keys2 = transcriptomics_dict.keys()
        keys3 = bulk_dict.keys()
        tests = set(keys2).intersection(keys3)
        tests2 = (set(keys2).union(set(keys3))).difference(tests)
    elif transcript_file:
        keys2 = transcriptomics_dict.keys()
        tests = set(keys2)
        tests2 = tests
    elif bulk_file:
        keys3 = bulk_dict.keys()
        tests = set(keys3)
        tests2 = tests
    else:
        getopt.GetoptError
        print("No files given")
        sys.exit(2)
        
    #tests = set(keys1).intersection(set(keys2))
    print("keys1 generated")
    print(poop)
    merge_data = None
    for test in tests:
        print(test)
        if prote_file:
            prote_data = proteomics_dict[test].loc[:, ['expressed', 'top']]
            prote_data.rename(columns={'expressed': 'prote_exp',
                                       'top': 'prote_top'}, inplace=True)
            #prote_data = mergeLogicalTable(prote_data)
            print(prote_data)
            test = unidecode.unidecode(test)
            try:
                if not merge_data:
                    merge_data = prote_data
            except:
                merge_data = None

        if transcript_file:
            trans_data = transcriptomics_dict[test].loc[:, ['expressed', 'top']]
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
            #bulk_data = mergeLogicalTable(bulk_data)
            print(bulk_data)
            test = unidecode.unidecode(test)
            try:
                if not merge_data:
                    merge_data = bulk_data
            except:
                merge_data = merge_data.join(bulk_data, how='outer')

        print(merge_data)
        merge_data = mergeLogicalTable(merge_data)
        #test = unidecode.unidecode(test)
        print(merge_data)
        if prote_file and transcript_file and bulk_file:
            merge_data['Express'] = 0
            merge_data.loc[merge_data[['prote_exp', 'trans_exp', 'bulk_exp']].sum(axis=1) == 3, 'Express'] = 1
            merge_data.loc[merge_data[['prote_top', 'trans_top', 'bulk_top']].sum(axis=1) > 0, 'Express'] = 1

        elif prote_file and transcript_file:
            merge_data['Express'] = 0
            merge_data.loc[merge_data[['prote_exp', 'trans_exp']].sum(axis=1) == 2, 'Express'] = 1
            merge_data.loc[merge_data[['prote_top', 'trans_top']].sum(axis=1) > 0, 'Express'] = 1
        elif prote_file and bulk_file:
            merge_data['Express'] = 0
            merge_data.loc[merge_data[['prote_exp', 'bulk_exp']].sum(axis=1) == 2, 'Express'] = 1
            merge_data.loc[merge_data[['prote_top', 'bulk_top']].sum(axis=1) > 0, 'Express'] = 1
        elif prote_file:
            merge_data['Express'] = 0
            merge_data.loc[merge_data[['prote_exp']].sum(axis=1) == 1, 'Express'] = 1
            merge_data.loc[merge_data[['prote_top']].sum(axis=1) > 0, 'Express'] = 1
        elif trasnscript_file and bulk_file:
            merge_data['Express'] = 0
            merge_data.loc[merge_data[['trans_exp', 'bulk_exp']].sum(axis=1) == 2, 'Express'] = 1
            merge_data.loc[merge_data[['trans_top', 'bulk_top']].sum(axis=1) > 0, 'Express'] = 1
        elif transcript_file:
            merge_data['Express'] = 0
            merge_data.loc[merge_data[['trans_exp']].sum(axis=1) == 1, 'Express'] = 1
        elif bulk_file:
            merge_data['Express'] = 0
            merge_data.loc[merge_data[['bulk_exp']].sum(axis=1) == 1, 'Express'] = 1
            merge_data.loc[merge_data[['bulk_top']].sum(axis=1) > 0, 'Express'] = 1


        #merge_data = prote_data.join(trans_data, how='outer')
        #merge_data = mergeLogicalTable(merge_data)
        #merge_data['Express'] = 0
        #merge_data.loc[merge_data[['prote_0.5', 'trans_0.5']].sum(axis=1) == 2, 'Express'] = 1
        #merge_data.loc[merge_data[['prote_top', 'trans_top']].sum(axis=1) > 0, 'Express'] = 1

        filepath = os.path.join(configs.rootdir, 'data', 'merged_{}.csv'.format(test))
        merge_data.to_csv(filepath, index_label='ENTREZ_GENE_ID')

        filepath = os.path.join(configs.rootdir, 'data', 'GeneExpression_{}_Merged.csv'.format(test))
        merge_data.reset_index(drop=False, inplace=True)

        splitEntrez = splitGeneExpressionData(merge_data)
        splitEntrez.rename(columns={'Gene':'ENTREZ_GENE_ID', 'Data':'Express'}, inplace=True)
        splitEntrez.to_csv(filepath, index_label='ENTREZ_GENE_ID')
        files_dict[test] = filepath

        print('{}: save to {}\n'.format(test, filepath))

    # Those only with transcriptomics or proteomics
    #tests2 = (set(keys1).union(set(keys2))).difference(tests)
    #for test in tests2:
    while False:
        test = unidecode.unidecode(test)
        if prote_file:
            prote_data = proteomics_dict[test].loc[:, ['expressed', 'top']]
            prote_data.rename(columns={'expressed': 'prote_exp',
                                       'top': 'prote_top'}, inplace=True)
            print(prote_data.shape)
            test = unidecode.unidecode(test)
            if not merge_data:
                merge_data = prote_data

        if transcript_file:
            trans_data = transcriptomics_dict[test].loc[:, ['expressed', 'top']]
            trans_data.rename(columns={'expressed': 'trans_exp',
                                       'top': 'trans_top'}, inplace=True)
            print(trans_data.shape)
            test = unidecode.unidecode(test)
            if not merge_data:
                merge_data = trans_data
            else:
                merge_data = merge_data.join(trans_data, how='outer')

        if bulk_file:
            bulk_data = bulk_dict[test].loc[:, ['expressed', 'top']]
            bulk_data.rename(columns={'expressed': 'bulk_exp',
                                      'top': 'bulk_top'}, inplace=True)
            print(bulk_data.shape)
            test = unidecode.unidecode(test)
            if not merge_data:
                merge_data = bulk_data
            else:
                merge_data = merge_data.join(bulk_data, how='outer')

        #test = unidecode.unidecode(test)
        if prote_file and transcript_file and bulk_file:
            merge_data['Express'] = 0
            merge_data.loc[merge_data[['prote_exp', 'trans_exp', 'bulk_exp']].sum(axis=1) == 3, 'Express'] = 1
            merge_data.loc[merge_data[['prote_top', 'trans_top', 'bulk_top']].sum(axis=1) > 0, 'Express'] = 1
        elif prote_file and transcript_file:
            merge_data['Express'] = 0
            merge_data.loc[merge_data[['prote_exp', 'trans_exp']].sum(axis=1) == 2, 'Express'] = 1
            merge_data.loc[merge_data[['prote_top', 'trans_top']].sum(axis=1) > 0, 'Express'] = 1
        elif prote_file and bulk_file:
            merge_data['Express'] = 0
            merge_data.loc[merge_data[['prote_exp', 'bulk_exp']].sum(axis=1) == 2, 'Express'] = 1
            merge_data.loc[merge_data[['prote_top', 'bulk_top']].sum(axis=1) > 0, 'Express'] = 1
        elif prote_file:
            merge_data['Express'] = 0
            merge_data.loc[merge_data[['prote_exp']].sum(axis=1) == 1, 'Express'] = 1
            merge_data.loc[merge_data[['prote_top']].sum(axis=1) > 0, 'Express'] = 1
        elif trasnscript_file and bulk_file:
            merge_data['Express'] = 0
            merge_data.loc[merge_data[['trans_exp', 'bulk_exp']].sum(axis=1) == 2, 'Express'] = 1
            merge_data.loc[merge_data[['trans_top', 'bulk_top']].sum(axis=1) > 0, 'Express'] = 1
        elif transcript_file:
            merge_data['Express'] = 0
            merge_data.loc[merge_data[['trans_exp']].sum(axis=1) == 1, 'Express'] = 1
        elif bulk_file:
            merge_data['Express'] = 0
            merge_data.loc[merge_data[['bulk_exp']].sum(axis=1) == 1, 'Express'] = 1
            merge_data.loc[merge_data[['bulk_top']].sum(axis=1) > 0, 'Express'] = 1
        else:
            print("Unknown Error")
            return None

        merge_data = mergeLogicalTable(merge_data)


        old = '''
        if test in keys1:
            merge_data = proteomics_dict[test].loc[:, ['0.5', 'top 0.25']]
            merge_data.rename(columns={'0.5': '0.5', 'top 0.25': 'top'}, inplace=True)
        elif test in keys2:
            merge_data = transcriptomics_dict[test].loc[:, ['0.5', '0.9']]
            merge_data.rename(columns={'0.5': '0.5', '0.9': 'top'}, inplace=True)
        else:
            print('Unknown Error')
            return None

        test = unidecode.unidecode(test)

        merge_data = mergeLogicalTable(merge_data)
        merge_data['Express'] = 0
        merge_data.loc[merge_data[['0.5', 'top']].sum(axis=1) > 0, 'Express'] = 1
        filepath = os.path.join(configs.rootdir, 'data', 'merged_{}.csv'.format(test))
        merge_data.to_csv(filepath, index_label='ENTREZ_GENE_ID')
        '''

        filepath = os.path.join(configs.rootdir, 'data', 'GeneExpression_{}_Merged.csv'.format(test))
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
    try:
        opts, args = getopt.getopt(argv, "ht:p:b:", ["transfile=", "protefile=", "bulkfile="])
    except getopt.GetoptError:
        print('merge_xomics.py -t <transfile> -p <protefile>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('merge_xomics.py -t <transfile> -p <protefile>')
            sys.exit()
        elif opt in ("-t", "--transfile"):
            transfile = arg
        elif opt in ("-p", "--protefile"):
            protefile = arg
        elif opt in ('-b', "--bulkfile"):
            bulkfile = arg
    print('Transcriptomics file is "{}"'.format(transfile))
    print('Proteomics file is "{}"'.format(protefile))
    print('Bulk RNA-seq file is "{}"'.format(bulkfile))
    files_dict = merge_xomics(transcript_file=transfile,
                              prote_file=protefile,
                              bulk_file=bulkfile)
    print(files_dict)
    files_json = os.path.join(configs.rootdir, 'data', 'step1_results_files.json')
    with open(files_json, 'w') as fp:
        json.dump(files_dict, fp)


if __name__ == "__main__":
   main(sys.argv[1:])
