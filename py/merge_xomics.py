#!/usr/bin/python3
import os
import sys
import time
import pandas as pd
import numpy as np
# from bioservices import BioDBNet
from project import configs
from transcriptomic_gen import *
from proteomics_gen import *

# input parameters
def merge_xomics(transcript_file='GeneExpressionDataUsed.xlsx', prote_file='Supplementary Data 1.xlsx'):
    transcriptomics_dict= load_transcriptomics_tests(filename = transcript_file)
    Proteomics = load_supplementary_data(prote_file)
    proteomics_dict = load_proteomics_tests(Proteomics)

    files_dict = dict()
    keys1 = proteomics_dict.keys()
    keys2 = transcriptomics_dict.keys()
    tests = set(keys1).intersection(set(keys2))
    for test in tests:
        prote_data = proteomics_dict[test].loc[:, ['0.5', 'top 0.25']]
        prote_data.rename(columns={'0.5': 'prote_0.5', 'top 0.25': 'prote_top'}, inplace=True)

        trans_data = transcriptomics_dict[test].loc[:, ['0.5', '0.9']]
        trans_data.rename(columns={'0.5': 'trans_0.5', '0.9': 'trans_top'}, inplace=True)

        print(trans_data.shape)
        print(prote_data.shape)

        merge_data = prote_data.join(trans_data, how='outer')
        merge_data = mergeLogicalTable(merge_data)
        merge_data['Express'] = 0
        merge_data.loc[merge_data[['prote_0.5', 'trans_0.5']].sum(axis=1) == 2, 'Express'] = 1
        merge_data.loc[merge_data[['prote_top', 'trans_top']].sum(axis=1) > 0, 'Express'] = 1
        filepath = os.path.join(configs.rootdir, 'data', 'merged_{}.csv'.format(test))
        merge_data.to_csv(filepath, index_label='ENTREZ_GENE_ID')
        files_dict[test] = filepath

        print('{}: save to {}\n'.format(test, filepath))

    # Those only with transcriptomics or proteomics
    tests2 = (set(keys1).union(set(keys2))).difference(tests)
    for test in tests2:
        if test in keys1:
            merge_data = proteomics_dict[test].loc[:, ['0.5', 'top 0.25']]
            merge_data.rename(columns={'0.5': '0.5', 'top 0.25': 'top'}, inplace=True)
        elif test in keys2:
            merge_data = transcriptomics_dict[test].loc[:, ['0.5', '0.9']]
            merge_data.rename(columns={'0.5': '0.5', '0.9': 'top'}, inplace=True)
        else:
            print('Unknown Error')
            return None

        merge_data = mergeLogicalTable(merge_data)
        merge_data['Express'] = 0
        merge_data.loc[merge_data[['0.5', 'top']].sum(axis=1) > 0, 'Express'] = 1
        filepath = os.path.join(configs.rootdir, 'data', 'merged_{}.csv'.format(test))
        merge_data.to_csv(filepath, index_label='ENTREZ_GENE_ID')
        files_dict[test] = filepath

        print('{}: save to {}\n'.format(test, filepath))

    return files_dict

def main(argv):
    files_dict = merge_xomics(transcript_file='GeneExpressionDataUsed.xlsx',
                              prote_file='Supplementary Data 1.xlsx')
    print(files_dict)

if __name__ == "__main__":
   main(sys.argv[1:])