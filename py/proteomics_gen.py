#!/usr/bin/python3

import os
import sys
import time
import pandas as pd
import numpy as np
# from bioservices import BioDBNet
from transcriptomic_gen import *

# Load Proteomics
def load_proteomics_data(datafilename):
    data_sheet = list(range(1)) # first 5 sheets
    dataFullPath = os.path.join(projectdir, 'data', datafilename)
    if os.path.isfile(dataFullPath):
        fulldata = pd.read_excel(dataFullPath, sheet_name=data_sheet, header=0)
    else:
        print("Error: file not found: {}".format(dataFullPath))
        return None
    proteomics_data = fulldata[0]

    # Preprocess data, drop na, duplicate ';' in Gene names,
    proteomics_data['Gene names'] = proteomics_data['Gene names'].astype(str)
    proteomics_data.dropna(subset=['Gene names'],inplace=True)
    pluralnames = proteomics_data[proteomics_data['Gene names'].str.contains(';')==True]
    for idx, row in pluralnames.iterrows():
        # print('{}: {}'.format(idx, row['Gene names']))
        names = row['Gene names'].split(';')
        rows = []
        for name in names:
            rowcopy = row.copy()
            rowcopy['Gene names'] = name
            rows.append(rowcopy)
        proteomics_data.drop(index=idx,inplace=True)
        proteomics_data.append(rows, ignore_index=True)

    proteomics_data.rename(columns={'Gene names':'Gene Symbol'}, inplace=True)
    return proteomics_data


def load_supplementary_data(suppfilename):
    suppFullPath = os.path.join(projectdir, 'data', suppfilename)
    supp_sheet = ['Proteomics']
    supplements = pd.read_excel(suppFullPath, sheet_name=supp_sheet, header=0)
    Proteomics = supplements['Proteomics']
    return Proteomics


def load_gene_symbol_map(gene_symbols, filename = 'proteomics_entrez_map.csv'):
    filepath = os.path.join(projectdir, 'data', 'proteomics_entrez_map.csv')
    if os.path.isfile(filepath):
        sym2id = pd.read_csv(filepath, index_col='Gene Symbol')
    else:
        sym2id = fetch_entrez_gene_id(gene_symbols, input_db='Gene Symbol')
        sym2id.loc[sym2id['Gene ID'] == '-', ['Gene ID']] = np.nan
        sym2id.to_csv(filepath, index_label='Gene Symbol')
    return sym2id


def save_proteomics_tests(Proteomics, proteomics_data):
    tests = []
    files = []
    datas = []
    for test in list(Proteomics):
        if test == 'Statistics':
            break
        cols = Proteomics[test].dropna().tolist()
        testdata = proteomics_data.loc[:, cols]
        # Logical Calculation
        thresholds = testdata.quantile(0.75, axis=0)
        testbool = pd.DataFrame(0, columns=list(testdata), index=testdata.index)
        for col in list(testdata):
            testbool.loc[testdata[col] > thresholds[col], [col]] = 1

        testdata['pos'] = (testdata > 0).sum(axis=1) / testdata.count(axis=1)
        testdata['0.5'] = 0
        testdata.loc[(testdata['pos'] > 0.5), ['0.5']] = 1
        testdata['top 0.25'] = 0
        testdata.loc[testbool.any(axis=1), 'top 0.25'] = 1

        output_path = os.path.join(projectdir, 'data', 'Proteomics_{}.csv'.format(test.strip()))
        testdata.to_csv(output_path, index_label='ENTREZ_GENE_ID')
        print('Test Data Saved to {}'.format(output_path))
        tests.append(test.strip())
        files.append(output_path)
        datas.append(testdata)

    proteomics_dict = dict(zip(tests, files))
    testdata_dict = dict(zip(tests, datas))
    return proteomics_dict, testdata_dict


# Read data from csv files
def read_proteomics_tests(Proteomics):
    tests = []
    datas = []
    for test in list(Proteomics):
        if test == 'Statistics':
            break
        output_path = os.path.join(projectdir, 'data', 'Proteomics_{}.csv'.format(test.strip()))
        testdata = pd.read_csv(output_path, index_col='ENTREZ_GENE_ID')
        print('Test Data Load From {}'.format(output_path))
        tests.append(test.strip().replace('ï', 'i'))
        datas.append(testdata)

    proteomics_dict = dict(zip(tests, datas))
    return proteomics_dict


def main(argv):
    datafilename = 'ni.3693-S5.xlsx'
    suppfilename = 'Supplementary Data 1.xlsx'
    proteomics_data = load_proteomics_data(datafilename)
    Proteomics = load_supplementary_data(suppfilename)
    sym2id = load_gene_symbol_map(gene_symbols=proteomics_data['Gene Symbol'].tolist(),
                                  filename='proteomics_entrez_map.csv')

    # map gene symbol to ENTREZ_GENE_ID
    proteomics_data.set_index('Gene Symbol', inplace=True)
    proteomics_data['ENTREZ_GENE_ID'] = sym2id['Gene ID']
    proteomics_data.dropna(subset=['ENTREZ_GENE_ID'], inplace=True)
    proteomics_data.set_index('ENTREZ_GENE_ID', inplace=True)
    proteomics_data.drop(columns=['Majority protein IDs'], inplace=True)
    prote_data_filepath = os.path.join(projectdir, 'data', 'proteomics_data.csv')
    if not os.path.isfile(prote_data_filepath):
        proteomics_data.to_csv(prote_data_filepath, index_label='ENTREZ_GENE_ID')

    # save proteomics data by test
    proteomics_dict, testdata_dict = save_proteomics_tests(Proteomics, proteomics_data)
    return True


if __name__ == "__main__":
   main(sys.argv[1:])