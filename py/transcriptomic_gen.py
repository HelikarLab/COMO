#!/usr/bin/python3
import re
import os
import sys

import pandas as pd
import numpy as np
from scipy import stats
from sqlalchemy import Column, ForeignKey, Integer, String, Float
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import relationship, sessionmaker, load_only
from sqlalchemy import create_engine

from GSEpipelineFast import *

#create declarative_base instance
Base = declarative_base()

#creates a create_engine instance at the bottom of the file
engine = create_engine('sqlite:///transcriptomics.db')

#we'll add classes here
class Sample(Base):
    __tablename__ = 'sample'

    id = Column(Integer, primary_key=True)
    ENTREZ_GENE_ID = Column(String(250), nullable=False)
    VALUE = Column(Float, nullable=False)
    P_VALUE = Column(Float)
    ABS_CALL = Column(String(2))
    Sample = Column(String(64))

class IDMapps(Base):
    __tablename__ = 'id_entrez_map'

    idx = Column(Integer, primary_key=True)
    ID = Column(String(250))
    ENTREZ = Column(String(250))


class GSEinfo(Base):
    __tablename__ = 'gseinfo'

    Sample = Column(String(64), primary_key=True)
    GSE = Column(String(64))
    Platform = Column(String(64))

Base.metadata.create_all(engine)

DBSession = sessionmaker(bind=engine)
session = DBSession()

# find project root dir
currentdir = os.getcwd()
dirlist = currentdir.split('/')
projectdir = '/'.join(dirlist[0:-1])


def lookupTranscriptomicsDB(gseXXX):
    '''
    check if gse already in database
    :param gseXXX:
    :return:
    '''
    df = pd.read_sql(session.query(GSEinfo).filter_by(GSE=gseXXX.gsename).options(load_only("Sample", "GSE")).statement,
                     session.bind)
    gsm_list = gseXXX.gsm_platform.keys()
    gsm_db = df.Sample.tolist()
    if set(gsm_list).issubset(set(gsm_db)):
        return True
    else:
        return False
    # df = pd.read_sql(
    #     session.query(Sample).filter_by(Sample='GSM60728').options(load_only("ABS_CALL", "ENTREZ_GENE_ID")).statement,
    #     session.bind,
    #     index_col='ENTREZ_GENE_ID'
    #     )
    # df.drop('id', axis=1, inplace=True)


def updateTranscriptomicsDB(gseXXX):
    '''
    Update GSE info to transcriptomics.db
    :param gseXXX: gse object
    :return:
    '''
    # df_clean = gseXXX.get_entrez_table_pipeline()
    df_clean = gseXXX.get_entrez_table_pipeline()
    if df_clean.empty:
        return None

    df_clean.sort_index(inplace=True)
    df_clean.reset_index(inplace=True)

    # write to database table sample
    df_samples = pd.DataFrame([],columns=['ENTREZ_GENE_ID','VALUE','P_VALUE','ABS_CALL','Sample'])
    df_samples.index.name = 'id'

    cols_clean = list(df_clean)
    for key,val in gseXXX.gsm_platform.items():
        col_val = '{}.cel.gz'.format(key.lower())
        col_abs = '{}.cel.gz.1'.format(key.lower())
        col_p = '{}.cel.gz.2'.format(key.lower())
        if not col_val in cols_clean:
            continue
        df_s = pd.DataFrame([])
        df_s['ENTREZ_GENE_ID'] = df_clean['ENTREZ_GENE_ID']
        df_s['VALUE'] = df_clean[col_val]
        df_s['ABS_CALL'] = df_clean[col_abs]
        df_s['P_VALUE'] = df_clean[col_p]
        df_s['Sample'] = key.upper()
        df_s.index.name='id'
        df_samples = pd.concat([df_samples, df_s.dropna(how='any')], ignore_index=True, sort=True)

    df_samples.index.name = 'id'
    df_samples.to_sql(con=engine,name='sample',if_exists='append',index=False)

    # write to database table GSEinfo
    df_gseinfo = pd.DataFrame([],columns=['Sample','GSE','Platform'])
    df_gseinfo['Sample'] = pd.Series(list(gseXXX.gsm_platform.keys()))
    df_gseinfo['Platform'] = pd.Series(list(gseXXX.gsm_platform.values()))
    df_gseinfo['GSE'] = gseXXX.gsename

    df_gseinfo.to_sql(con=engine,name='gseinfo',if_exists='append',index=False)


# function to complete the inquery of a sheet
def queryTest(df):
    sr = df['GSE ID']
    gse_ids = sr[sr.str.match('GSE')].unique()
    sr = df['Samples'].dropna()
    gsm_ids = sr.unique()
    print('---\nStart Collecting Data for:')
    print(gse_ids)
    print(gsm_ids)
    print('---\n')
    # fetch data of each gse if it is not in the database, update database
    for gse_id in gse_ids:
        querytable = df[df['GSE ID']==gse_id]
        gseXXX = GSEproject(gse_id,querytable,projectdir)
        if lookupTranscriptomicsDB(gseXXX):
            print("{} already in database, skip over.".format(gseXXX.gsename))
            continue
        updateTranscriptomicsDB(gseXXX)

    df_results = fetchLogicalTable(gsm_ids)
    df_output = mergeLogicalTable(df_results)
    return df_output


def fetchLogicalTable(gsm_ids):
    '''
    Fetch the Logical Table Based on ABS_CALL of Samples
    :param gsm_ids: list of sample names to fetch
    :return: pandas dataframe
    '''
    df_results = pd.DataFrame([], columns=['ENTREZ_GENE_ID'])
    df_results.set_index('ENTREZ_GENE_ID', inplace=True)
    for gsm in gsm_ids:
        # print(gsm)
        df = pd.read_sql(
            session.query(Sample).filter_by(Sample=gsm).options(load_only("ABS_CALL", "ENTREZ_GENE_ID")).statement,
            session.bind,
            index_col='ENTREZ_GENE_ID')
        df.drop('id', axis=1, inplace=True)
        df.rename(columns={'ABS_CALL': gsm}, inplace=True)
        df.loc[df[gsm] == 'A', gsm] = 0
        df.loc[df[gsm] == 'P', gsm] = 1
        df.loc[df[gsm] == 'M', gsm] = 1
        df.sort_values(by=['ENTREZ_GENE_ID',gsm],inplace=True)
        df = df[~df.index.duplicated(keep='last')]

        df_results = pd.concat([df_results, df], axis=1, sort=False)

    # Need to set index name after merge
    df_results.index.name = 'ENTREZ_GENE_ID'
    return df_results


## Merge Output
def mergeLogicalTable(df_results):
    '''
    Merge the Rows of Logical Table belongs to the same ENTREZ_GENE_ID
    :param df_results:
    :return: pandas dataframe of merged table
    '''
    # step 1: get all plural ENTREZ_GENE_IDs in the input table, extract unique IDs
    df_results.reset_index(drop=False, inplace=True)
    df_results['ENTREZ_GENE_ID'] = df_results['ENTREZ_GENE_ID'].str.replace(" /// ", "//")
    id_list = []
    entrez_single_id_list = df_results[~df_results['ENTREZ_GENE_ID'].str.contains("//")]['ENTREZ_GENE_ID'].tolist()
    entrez_id_list = df_results[df_results['ENTREZ_GENE_ID'].str.contains("//")]['ENTREZ_GENE_ID'].tolist()
    for entrez_id in entrez_id_list:
        entrez_ids = entrez_id.split("//")
        id_list.extend(entrez_ids)
        df_dups = pd.DataFrame([], columns=list(df_results), index=list(range(len(entrez_ids))))
        dup_rows = []
        for eid in entrez_ids:
            rows = df_results.loc[df_results['ENTREZ_GENE_ID']==entrez_id].copy()
            rows['ENTREZ_GENE_ID'] = eid
            dup_rows.append(rows)
        df_results = df_results.append(dup_rows,ignore_index=True)
        df_results.drop(df_results[df_results['ENTREZ_GENE_ID']==entrez_id].index, inplace=True)

    # step 2: print out information about merge
    common_elements = list(set(entrez_single_id_list).intersection(set(id_list)))
    # information of merge
    print('{} single ENTREZ_GENE_IDs to merge'.format(len(common_elements)))
    print('id_list: {}, set: {}'.format(len(id_list),len(set(id_list))))
    print('entrez_single_id_list: {}, set: {}'.format(len(entrez_single_id_list),len(set(entrez_single_id_list))))
    print('entrez_id_list: {}, set: {}'.format(len(entrez_id_list),len(set(entrez_id_list))))

    dups = [x for x in id_list if id_list.count(x) > 1]
    # dups = list(set(id_list))
    print('dups: {}, set: {}'.format(len(dups),len(set(dups))))

    full_entre_id_sets = []
    cnt = 0
    entrez_dups_list = []
    idx_list = list(range(len(entrez_id_list)))
    for idx1 in range(len(entrez_id_list)):
        if not idx1 in idx_list:
            continue
        set1 = set(entrez_id_list[idx1].split("//"))
        idx_list.remove(idx1)
        toremove = []
        for idx2 in idx_list:
            set2 = set(entrez_id_list[idx2].split("//"))
            intersect = set1.intersection(set2)
            if bool(intersect):
                set1 = set1.union(set2)
                toremove.append(idx2)
        for idx3 in toremove:
            idx_list.remove(idx3)
        sortlist = list(set1)
        sortlist.sort(key=int)
        new_entrez_id = " /// ".join(sortlist)
        full_entre_id_sets.append(new_entrez_id)
    full_entre_id_sets = list(set(full_entre_id_sets))

    for full_entrez_id in full_entre_id_sets:
        singles = full_entrez_id.split(" /// ")
        entrez_dups_list.append(singles)
        cnt += 1

    print('{} id merged'.format(cnt))
    entrez_dups_dict = dict(zip(full_entre_id_sets,entrez_dups_list))
    # full_entre_id_sets = list(set(full_entre_id_sets))

    # df_results.reset_index(inplace=True)
    for merged_entrez_id, entrez_dups_list in entrez_dups_dict.items():
        df_results['ENTREZ_GENE_ID'].replace(to_replace=entrez_dups_list,
                                             value=merged_entrez_id,
                                             inplace=True)
    # df_results.drop_duplicates(subset=['ENTREZ_GENE_ID'], keep='first', inplace=True)
    df_results.set_index('ENTREZ_GENE_ID', inplace=True)

    df_output = df_results.fillna(-1).groupby(level=0).max()
    df_output.replace(-1, np.nan, inplace=True)
    posratio = df_output.sum(axis=1,skipna=True)/df_output.count(axis=1)
    df_output['Pos'] = posratio
    df_output['0.5'] = np.where(posratio >= 0.5, 1, 0)
    df_output['0.9'] = np.where(posratio >= 0.9, 1, 0)
    return df_output


def main(args):
    # input from user
    filename = 'GeneExpressionDataUsed.xlsx'
    sheet_name = list(range(5)) # first 5 sheets

    # gse2770 = GSEproject('GSE2770',projectdir)

    inqueryFullPath = os.path.join(projectdir, 'data', filename)
    inqueries = pd.read_excel(inqueryFullPath, sheet_name=sheet_name, header=0)
    
    for i in sheet_name:
        # print(list(inqueries[i]))
        inqueries[i].fillna(method='ffill',inplace=True)
        df = inqueries[i].loc[:,['GSE ID','Samples','GPL ID','Instrument']]
        df_output = queryTest(df)
        filename = 'logicaltable_sheet_{}.csv'.format(i+1)
        fullsavepath = os.path.join(projectdir,'data',filename)
        df_output.to_csv(fullsavepath)
        print('Save to {}'.format(fullsavepath))

if __name__ == "__main__":
   main(sys.argv[1:])