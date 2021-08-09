#!/usr/bin/python3

import re
import os
import pandas as pd
import numpy as np
#import GEOparse
import urllib.request
import tarfile
from instruments import *
from GSEpipeline import load_gse_soft

from rpy2.robjects import pandas2ri
pandas2ri.activate()
# Input: Extract Gene Info from GEO DataSets

# gse = load_gse_soft(gsename)

# Extract Platform Information

def download_gsm_id_maps(datadir, gse, gpls = ['GPL96','GPL97','GPL8300'],vendor='affy'):
    '''
    download ID to ENTREZ_GENE_ID maps, create a csv file for each platform, and return dictionary
    :param gse: gse object
    :param datadir: path to save the csv files
    :param platforms: list of platforms
    :return: dictionary of maps indexed by platform
    '''
    # maps_list = []
    #gene_maps = pd.DataFrame([],columns=['GPL96','GPL97','GPL8300','ENTREZ_GENE_ID'])
    #gene_maps.set_index('ENTREZ_GENE_ID',inplace=True)
    for gpl in gpls:
        table =gse.gpls[gpl].table.copy()
        if vendor.lower() == 'affy':
            temp = table[['ID','ENTREZ_GENE_ID']]
        elif vendor.lower() == 'agilent':
            input_values = table.loc[table['CONTROL_TYPE']=='FALSE','SPOT_ID'].tolist()
            temp = fetch_entrez_gene_id(input_values)
            temp.drop(columns=['Ensembl Gene ID'],inplace=True)
            temp.reset_index(inplace=True)
            temp.rename(columns={"Agilent ID":"ID", "Gene ID":"ENTREZ_GENE_ID"}, inplace=True)
            temp.replace(to_replace='-', value=np.nan, inplace=True)
        else:
            print('Unsupported Platform: {}'.format(gpl))
            continue
        # Save to file
        filefullpath = os.path.join(datadir,'{}entrez.csv'.format(gpl.lower()))
        print(filefullpath)
        temp.to_csv(filefullpath, index=False)
        # Single Table
        # temp.dropna(axis=0,inplace=True)
        # temp.set_index('ENTREZ_GENE_ID',inplace=True)
        # maps_list.append(temp)

    # maps_dict = dict(zip(platforms, maps_list))
    # return maps_dict





# df_outer = create_entrez_table_default(gse)


# initialize with platform
# gsm_maps = get_gsm_tables(gse)

class GSEproject:
    def __init__(self,gsename, querytable, rootdir='../'):
        self.gsename = gsename
        ## Setup paths
        self.querytable = querytable
        self.rootdir = rootdir
        self.datadir = os.path.join(self.rootdir,'data')
        self.outputdir = os.path.join(self.rootdir,'output')
        self.genedir = os.path.join(self.datadir,self.gsename + '_RAW')
        print('Initialize project ({}):\nRoot: {}\nRaw data: {}'.format(self.gsename, self.rootdir, self.genedir))
        self.gsm_platform = self.get_gsm_platform()
        pairs = querytable.loc[:,['GPL ID','Instrument']].drop_duplicates()
        gpls = pairs['GPL ID'].tolist()
        vendors = pairs['Instrument'].tolist()
        self.platforms = dict(zip(gpls, vendors))
        self.download_samples()

    def organize_gse_raw_data(self):
        """
        Organize raw data at local folder
        :return:
        """
        # create a folder for each platform
        for key in self.platforms.keys():
            platformdir = os.path.join(self.genedir, key)
            if not os.path.exists(platformdir):
                os.makedirs(platformdir)
                print('Path created: {}'.format(platformdir))
            else:
                print('Path already exist: {}'.format(platformdir))

        # Move Corresponding Cel files to Folders
        onlyfiles = [f for f in os.listdir(self.genedir) if f.endswith('.gz')]
        cnt = 0
        for file in onlyfiles:
            filelist = file.split('.')
            prefix = filelist[0]
            if prefix in self.gsm_platform:
                platform = self.gsm_platform[prefix]
                platformdir = os.path.join(self.genedir, platform)
                src_path = os.path.join(self.genedir, file)
                dst_path = os.path.join(platformdir, file)
                os.rename(src_path, dst_path)
                print('Move {} to {}'.format(src_path, dst_path))
                cnt += 1
        print('{} raw data files moved.'.format(cnt))

    def get_gsm_tables(self):
        '''
        get gsm maps in table
        :param gse:
        :return:
        '''
        gsm_tables = {}
        for gpl,vendor in self.platforms.items():
            filename = '{}entrez.csv'.format(gpl.lower())
            filepath = os.path.join(self.datadir, filename)
            if not os.path.isfile(filepath):
                # Could improve to automatic download new tables based on platform
                gse = load_gse_soft(self.gsename)
                download_gsm_id_maps(self.datadir, gse, gpls=[gpl], vendor=vendor)
                print('Skip Unsupported Platform: {}, {}'.format(gpl, vendor))
                # continue
            temp = pd.read_csv(filepath)
            df = temp[['ID', 'ENTREZ_GENE_ID']]
            df.set_index('ID', inplace=True)
            # df.drop_duplicates(keep='last', inplace=True)
            gsm_tables[gpl] = df
        return gsm_tables

    def get_gsm_platform(self):
        '''
        Classify the Samples by Platform
        get the platform of each gsm
        :param gse: gse object
        :return: dictionary key: gsm, value: platform such as 'GEL96', 'GEL97', 'GEL8300'
        '''
        keys = self.querytable['Samples'].str.upper().tolist()
        values = self.querytable['GPL ID'].str.upper().tolist()
        gsm_platform = dict(zip(keys, values))
        return gsm_platform

    def gsms_included_by(self, df):
        for gsm in self.gsm_platform.keys():
            included = False
            gsml = gsm.lower()
            for key in list(df):
                if gsml == key.lower().split('.')[0]:
                    included = True
                    break
            if not included:
                return False
        return True

    def get_entrez_table_pipeline(self, fromcsv=True):
        '''
        create ENTREZ ID based table from gse
        :param gse: gse object
        :return: pandas dataframe for table of GSE
        '''
        print("ENTREZ_TABLE PIPE")
        filefullpath = os.path.join(self.genedir, '{}_sc500_full_table.csv'.format(self.gsename))
        print(filefullpath)
        if fromcsv and os.path.isfile(filefullpath):
            try:
                df_clean_sc500 = pd.read_csv(filefullpath)
                df_clean_sc500.dropna(axis='columns', how='all', inplace=True)
                df_clean_sc500.dropna(how='all', inplace=True)
                df_clean_sc500 = df_clean_sc500[~df_clean_sc500.index.duplicated(keep='last')]
                print(df_clean_sc500)
                if self.gsms_included_by(df_clean_sc500) and not df_clean_sc500.empty:
                    return df_clean_sc500
                else:
                    print("Need Append GSMs")
            except:
                print("Unable to read {}")
                # self.download_raw(overwrite=True)

        print('Create new table: {}'.format(filefullpath))
        gsm_maps = self.get_gsm_tables()
        if not any(gsm_maps):
            print("Not available, return empty dataframe")
            return pd.DataFrame([])

        # step 1: Ready Affy files from folders
        gsm_tables_sc500 = {}
        for key, vendor in self.platforms.items():
            platformdir = os.path.join(self.genedir, key)
            print('{} Read Path: {}'.format(vendor, platformdir))
            if os.path.exists(platformdir):
                if vendor.lower() == 'affy':
                    outputdf = affyio.readaffydir(platformdir)
                elif vendor.lower() == 'agilent':
                    outputdf = readagilent(platformdir, list(self.gsm_platform.keys()))
                else:
                    print("Unsupported Platform {} and Vendor {}".format(key,vendor))
                    continue
            else:
                print('Path not exist: {}'.format(platformdir))
                continue
            outputdf = ro.conversion.rpy2py(outputdf)
            outputdf['ENTREZ_GENE_ID'] = gsm_maps[key]['ENTREZ_GENE_ID']
            print(outputdf)
            gsm_tables_sc500[key] = outputdf
            print(gsm_tables_sc500[key])
            print(gsm_tables_sc500)

        # step 2: Drop rows without ENTREZ GENE ID, set index to ENTREZ
        for key in self.platforms.keys():
            gsm_tables_sc500[key].dropna(subset=['ENTREZ_GENE_ID'], inplace=True)
            gsm_tables_sc500[key].set_index('ENTREZ_GENE_ID', inplace=True)

        # step 3: Merge tables of platforms
        df_outer_sc500 = None
        for key in self.platforms.keys():
            print('{}: {}'.format(key, gsm_tables_sc500[key].shape))
            if df_outer_sc500 is None:
                df_outer_sc500 = gsm_tables_sc500[key]
            else:
                df_outer_sc500 = pd.merge(df_outer_sc500, gsm_tables_sc500[key], on='ENTREZ_GENE_ID', how='outer')

        df_outer_sc500.dropna(how='all', inplace=True)
        print('Full: {}'.format(df_outer_sc500.shape))
        df_outer_sc500.rename(str.lower, axis='columns', inplace=True)
        keys = []
        vals = []
        gsms_loaded = []
        for col in list(df_outer_sc500):
            if '.cel.gz' in col:
                strs = col.split('.cel.gz')
                gsm = strs[0].split('_')[0]
                newcol = '{}.cel.gz{}'.format(gsm, strs[-1])
                vals.append(newcol)
                keys.append(col)
                gsms_loaded.append(gsm)
        df_outer_sc500.rename(columns=dict(zip(keys, vals)), inplace=True)
        gsms_loaded = list(set(gsms_loaded).union(set(self.gsm_platform.keys())))

        # step 4: Remove duplicated items, keep largest VALUE for each GSM
        if not 'df_clean_sc500' in locals():
            df_clean_sc500 = pd.DataFrame([], index=df_outer_sc500.index)
            df_clean_sc500 = df_clean_sc500[~df_clean_sc500.index.duplicated(keep='first')]
        elif df_clean_sc500.empty:
            df_clean_sc500 = pd.DataFrame([], index=df_outer_sc500.index)
            df_clean_sc500 = df_clean_sc500[~df_clean_sc500.index.duplicated(keep='first')]
        else:
            print(df_clean_sc500)
            df_clean_sc500.set_index('ENTREZ_GENE_ID', inplace=True)
            placeholder = pd.DataFrame([], columns=['placeholder'],index=df_outer_sc500.index)
            placeholder['placeholder'] = 0
            placeholder.index.name = 'ENTREZ_GENE_ID'
            df_clean_sc500 = pd.merge(df_clean_sc500, placeholder, on='ENTREZ_GENE_ID', how='outer')
            df_clean_sc500 = df_clean_sc500[~df_clean_sc500.index.duplicated(keep='last')]

        for key in gsms_loaded:
            key_low = key.lower()
            col1, col2, col3 = '{}.cel.gz'.format(key_low), '{}.cel.gz.1'.format(key_low), '{}.cel.gz.2'.format(key_low)
            try:
                temp = df_outer_sc500.loc[:, [col1, col2, col3]]
            except:
                if key in list(self.gsm_platform.keys()):
                    print('{} not in df_outer_sc500'.format(key))
                continue
            temp.sort_values(by=['ENTREZ_GENE_ID', col1], inplace=True)
            temp = temp[~temp.index.duplicated(keep='last')]
            df_clean_sc500[col1] = temp[col1]
            df_clean_sc500[col2] = temp[col2]
            df_clean_sc500[col3] = temp[col3]

        # if 'express' in list(df_outer_sc500):
        #     df_clean_sc500 = df_outer_sc500.sort_values(by=['ENTREZ_GENE_ID', 'express'])
        #     df_clean_sc500 = df_clean_sc500[~df_clean_sc500.index.duplicated(keep='last')]

        # step 5: save to csv file
        try:
            df_clean_sc500.set_index('ENTREZ_GENE_ID',inplace=True)
        except:
            pass
        df_clean_sc500.dropna(axis='columns', how='all', inplace=True)
        df_clean_sc500.dropna(how='all', inplace=True)
        try:
            df_clean_sc500.drop(columns=['placeholder'], inplace=True)
        except:
            pass
        df_clean_sc500.sort_index(inplace=True)
        df_clean_sc500.to_csv(filefullpath)
        print('Full table saved to:\n{}'.format(filefullpath))
        return df_clean_sc500

    def download_raw(self, overwrite=False):
        # check if path created
        if (not os.path.isdir(self.genedir)) or overwrite:
            os.makedirs(self.genedir,exist_ok=True)
            url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc={}&format=file".format(self.gsename)
            filefullpath = os.path.join(self.genedir,'{}_RAW.tar'.format(self.gsename))
            if not os.path.isfile(filefullpath):
                print('Download Raw File: {}'.format(filefullpath))
                urllib.request.urlretrieve(url, filefullpath)
            else:
                print("File already exist: {}".format(filefullpath))
            tfile = tarfile.open(filefullpath)
            tfile.extractall(path=self.genedir)
            os.remove(filefullpath)
            print('Remove Raw File: {}'.format(filefullpath))
            self.organize_gse_raw_data()
        else:
            pass

    def download_samples(self, overwrite=False):
        os.makedirs(self.genedir,exist_ok=True)
        for gsm,gpl in self.gsm_platform.items():
            platformdir = os.path.join(self.genedir,gpl)
            os.makedirs(platformdir,exist_ok=True)
            sample_url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc={}&format=file".format(gsm)
            filefullpath = os.path.join(self.genedir,'{}.tar'.format(gsm))
            if (not os.path.isfile(filefullpath)) or overwrite:
                urllib.request.urlretrieve(sample_url, filefullpath)
                print('Retrieve Sample: {}'.format(filefullpath))
            else:
                print("Sample exist: {}".format(filefullpath))
                continue
            tfile = tarfile.open(filefullpath)
            tfile.extractall(path=platformdir)
            # os.remove(filefullpath) # keep to avoid re-download
            print('Extract to: {}'.format(platformdir))
        print('Retrieve Samples Completed.')

    def calculate_z_score(self,df,to_csv=False):
        cols = list(df)
        result = pd.DataFrame([],index=df.index)
        for col in cols:
            if not '.gz.' in col:
                score = np.log2(df[col])
                result[col] = (score-score.mean())/score.std(ddof=0)
        if to_csv:
            filefullpath = os.path.join(self.genedir,'{}_data_z.csv'.format(self.gsename))
            result.to_csv(filefullpath)
        return result