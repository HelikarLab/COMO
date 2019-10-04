#!/usr/bin/python3
import re
import os
import pandas as pd
import numpy as np
from project import configs
from GSEpipelineFast import *


def get_entrez_id(regulated, outputFullPath):
    RA_UP = fetch_entrez_gene_id(list(regulated.index.values), input_db='Affy ID')
    RA_UP.replace(to_replace='-', value=np.nan, inplace=True)
    RA_UP.dropna(how='any', subset=['Gene ID'], inplace=True)
    RA_UP['Gene ID'].to_csv(outputFullPath, index=False)
    return RA_UP


def pharse_configs(inqueryFullPath):
    xl = pd.ExcelFile(inqueryFullPath)
    sheet_name = xl.sheet_names
    inqueries = pd.read_excel(inqueryFullPath, sheet_name=sheet_name, header=0)
    inqueries['Sheet1'].fillna(method='ffill',inplace=True)
    df = inqueries['Sheet1'].loc[:,['GSE ID','Samples','GPL ID','Instrument']]
    df_target = inqueries['Sheet1'].loc[:,['Samples','Experiment']]
    df_target.rename(columns={'Samples':'FileName','Experiment':'Condition'},inplace=True)
    df_target['FileName'] = df_target['FileName'].astype(str) + '.txt.gz'
    df_target['SampleNumber']= 1 + df_target.index.values
    df_target = df_target[['SampleNumber','FileName','Condition']]
    return df, df_target



def main(argv):
    print(configs.rootdir)
    filename = 'Disease_Gene_Analyzed.xlsx'
    inqueryFullPath = os.path.join(configs.rootdir, 'data', filename)
    querytable, df_target = pharse_configs(inqueryFullPath)
    targetdir = os.path.join(configs.datadir,'targets.txt')
    df_target.to_csv(targetdir, index=False, sep='\t')
    sr = querytable['GSE ID']
    gse_ids = sr[sr.str.match('GSE')].unique()
    GSE_ID = gse_ids[0]
    gseXXX = GSEproject(GSE_ID,querytable,configs.rootdir)
    for key,val in gseXXX.platforms.items():
        rawdir = os.path.join(gseXXX.genedir,key)
        print('{}:{}, {}'.format(key, val, rawdir))
        data2 = affyio.fitaffydir(rawdir, targetdir)

    data2['abs_logFC'] = data2['logFC'].abs()
    data2.sort_values(by='abs_logFC', ascending=False, inplace=True)
    regulated = data2[data2['abs_logFC']>=1.0]
    down_regulated = regulated[regulated['logFC']<0]
    up_regulated = regulated[regulated['logFC']>0]

    up_file = os.path.join(configs.datadir,'RA_UP_{}.txt'.format(GSE_ID))
    down_file = os.path.join(configs.datadir,'RA_DOWN_{}.txt'.format(GSE_ID))
    RA_UP = get_entrez_id(up_regulated, up_file)
    RA_DOWN = get_entrez_id(down_regulated, down_file)
    print(RA_DOWN)
    print(RA_UP)

    all_file = os.path.join(configs.datadir, 'RA_FULL_{}.txt'.format(GSE_ID))
    RA_FULL = get_entrez_id(data2, all_file)
    data2.index.name = 'Affy ID'
    data2['ENTREZ_GENE_ID'] = RA_FULL['Gene ID']
    data2.dropna(how='any', subset=['ENTREZ_GENE_ID'], inplace=True)
    raw_file = os.path.join(configs.datadir, 'Raw_Fit_{}.csv'.format(GSE_ID))
    print('Raw Data saved to\n{}'.format(raw_file))
    data2.to_csv(raw_file, index=False)

if __name__ == "__main__":
   main(sys.argv[1:])
