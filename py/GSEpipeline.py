#!/usr/bin/python3

import re
import os
import pandas as pd
import numpy as np
import GEOparse
from rpy2.robjects.packages import importr
from rpy2.robjects import r, pandas2ri
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage

# automatically convert ryp2 dataframe to Pandas dataframe
pandas2ri.activate()
# Initialize Rpy2 for Affy package
affy = importr("affy")
# Define R function for read affy
string = """
readaffydir <- function(addr){
    crd <- getwd()
    setwd(addr)
    mydata = ReadAffy()
    setwd(crd)
    eset = mas5(mydata)
    eset_PMA <- mas5calls(mydata)
    y <- data.frame(exprs(eset), exprs(eset_PMA), assayDataElement(eset_PMA, "se.exprs"))
    y <- y[,sort(names(y))]
    return(y)
}
"""
affyio = SignatureTranslatedAnonymousPackage(string, "affyio")


# Input: Extract Gene Info from GEO DataSets
def load_gse_soft(name='GSE2770'):
    '''
    Read GSE information from local soft file, otherwise read online.
    :param name: name of gse
    :return: gse object by GEOparse
    '''
    softfile = "./{}_family.soft.gz".format(name)

    if os.path.isfile(softfile):
        gse = GEOparse.get_GEO(filepath=softfile)
    else:
        gse = GEOparse.get_GEO(geo=name, destdir="./")

    return gse

# gse = load_gse_soft(gsename)

# Extract Platform Information
def get_platform_probe(gse):
    '''
    extract platform information
    :param gse: gse object
    :return: dictionary of platform name and probe
    '''
    keys = []
    values = []
    for gpl in gse.gpls:
        keys.append(gse.gpls[gpl].name)
        r1 = re.findall(r"\[.*?\]", gse.gpls[gpl].metadata['title'][0])
        values.append(r1[0][4:-1])
        print(gse.gpls[gpl].name)
        print(gse.gpls[gpl].metadata['title'][0])

    celformat = dict(zip(keys, values))
    return celformat



def download_gsm_id_maps(datadir, gse, platforms = ['GPL96','GPL97','GPL8300']):
    '''
    download ID to ENTREZ_GENE_ID maps, create a csv file for each platform, and return dictionary
    :param gse: gse object
    :param datadir: path to save the csv files
    :param platforms: list of platforms
    :return: dictionary of maps indexed by platform
    '''
    maps_list = []
    #gene_maps = pd.DataFrame([],columns=['GPL96','GPL97','GPL8300','ENTREZ_GENE_ID'])
    #gene_maps.set_index('ENTREZ_GENE_ID',inplace=True)
    for platform in platforms:
        table =gse.gpls[platform].table.copy()
        temp = table[['ID','ENTREZ_GENE_ID']]
        # Save to file
        filefullpath = os.path.join(datadir,'{}entrez.csv'.format(platform))
        print(filefullpath)
        temp.to_csv(filefullpath, index=False)
        # Single Table
        temp.dropna(axis=0,inplace=True)
        temp.set_index('ENTREZ_GENE_ID',inplace=True)
        maps_list.append(temp)

    maps_dict = dict(zip(platforms, maps_list))
    return maps_dict



def get_gsm_tables(gse):
    '''
    get gsm maps in table
    :param gse:
    :return:
    '''
    celformat = get_platform_probe(gse)
    gsm_tables = {}
    for key, val in celformat.items():
        temp = gse.gpls[key].table.copy()
        df = temp[['ID', 'ENTREZ_GENE_ID']]
        df.set_index('ID', inplace=True)
        # df.drop_duplicates(keep='last', inplace=True)
        gsm_tables[key] = df
    return gsm_tables


# df_outer = create_entrez_table_default(gse)


# initialize with platform
# gsm_maps = get_gsm_tables(gse)

class GSEproject:
    def __init__(self,gsename="GSE2770",rootdir='../'):
        self.gsename = gsename
        ## Setup paths
        self.rootdir = rootdir
        self.datadir = os.path.join(self.rootdir,'data')
        self.outputdir = os.path.join(self.rootdir,'output')
        self.genedir = os.path.join(self.datadir,self.gsename + '_RAW')
        print('Initialize project ({}):\nRoot: {}\nRaw data: {}'.format(self.gsename, self.rootdir, self.genedir))
        self.gse = load_gse_soft(self.gsename)
        self.celformat = get_platform_probe(self.gse)
        self.gsm_platform = self.get_gsm_platform()

    def organize_gse_raw_data(self):
        """
        Organize raw data at local folder
        :return:
        """
        # create a folder for each platform
        for key in self.celformat.keys():
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

    def get_gsm_platform(self):
        '''
        Classify the Samples by Platform
        get the platform of each gsm
        :param gse: gse object
        :return: dictionary key: gsm, value: platform such as 'GEL96', 'GEL97', 'GEL8300'
        '''
        keys = []
        values = []
        for gsm in self.gse.gsms:
            keys.append(self.gse.gsms[gsm].name)
            for key, val in self.celformat.items():
                r1 = re.findall(r"{}".format(val), self.gse.gsms[gsm].columns.loc['ID_REF', 'description'])
                if not r1:
                    pass
                else:
                    values.append(key)

        gsm_platform = dict(zip(keys, values))
        return gsm_platform

    # create 3 tables by platform
    def get_entrez_table_default(self, fromcsv=True):
        '''
        Create entrez table from online data
        :param gse: gse object
        :return:
        '''
        filefullpath = os.path.join(self.genedir, '{}_full_table.csv'.format(self.gsename))
        if fromcsv:
            df_outer = pd.read_csv(filefullpath)
            return df_outer
        else:
            print('Create new table: {}'.format(filefullpath))
        # step 1: initialize from ID maps
        # celformat = get_platform_probe(gse)
        gsm_tables = get_gsm_tables(self.gse)
        for key,val in self.celformat.items():
            gsm_tables[key].drop_duplicates(keep='last', inplace=True)

        # gsm_platform = get_gsm_platform(gse)
        # step 2: create tables by platform
        for key, val in self.gsm_platform.items():
            temp = self.gse.gsms[key].table.copy()
            temp.rename(columns={"ID_REF": "ID"}, inplace=True)
            temp.dropna(subset=['ID'], how='any', inplace=True)
            temp.set_index('ID', inplace=True)
            # col1,col2,col3 = '{}.VALUE'.format(key), '{}.ABS_CALL'.format(key), '{}.DETECTION P-VALUE'.format(key)
            col1, col2, col3 = '{}.CEL.gz'.format(key), '{}.CEL.gz.1'.format(key), '{}.CEL.gz.2'.format(key)
            # gsm_tables[val].loc[:,[col1,col2,col3]] = temp[['VALUE','ABS_CALL','DETECTION P-VALUE']]
            # Sort by VALUE and drop duplicated ID
            temp['ENTREZ_GENE_ID'] = gsm_tables[val]['ENTREZ_GENE_ID']
            temp.sort_values(by=['ENTREZ_GENE_ID', 'VALUE'], inplace=True)
            temp.drop_duplicates(subset=['ENTREZ_GENE_ID'], keep='last', inplace=True)
            gsm_tables[val][col1] = temp['VALUE']
            gsm_tables[val][col2] = temp['ABS_CALL']
            gsm_tables[val][col3] = temp['DETECTION P-VALUE']

        # step 3: drop NANs
        for key, val in self.celformat.items():
            # df = pd.DataFrame([],columns=['ID','ENTREZ_GENE_ID'])
            gsm_tables[key].dropna(subset=['ENTREZ_GENE_ID'], inplace=True)
            gsm_tables[key].set_index('ENTREZ_GENE_ID', inplace=True)
            # gsm_tables[key].dropna(how='all',inplace=True)
            # gsm_tables[key].drop_duplicates(keep='last',inplace=True)

        # step 4: merge tables
        df_outer = None
        for key, val in self.celformat.items():
            # df = pd.DataFrame([],columns=['ID','ENTREZ_GENE_ID'])
            print('{}: {}'.format(key, gsm_tables[key].shape))
            if df_outer is None:
                df_outer = gsm_tables[key]
            else:
                df_outer = pd.merge(df_outer, gsm_tables[key], on='ENTREZ_GENE_ID', how='outer')

        # step 5: save files
        df_outer.dropna(how='all', inplace=True)
        print('full : {}'.format(df_outer.shape))
        df_outer.sort_index(inplace=True)
        df_outer.to_csv(filefullpath)
        print('Full table saved to:\n{}'.format(filefullpath))
        return df_outer

    def get_entrez_table_pipeline(self, fromcsv=True):
        '''
        create ENTREZ ID based table from gse
        :param gse: gse object
        :return: pandas dataframe for table of GSE
        '''
        filefullpath = os.path.join(self.genedir, '{}_sc500_full_table.csv'.format(self.gsename))
        if fromcsv:
            df_clean_sc500 = pd.read_csv(filefullpath)
            return df_clean_sc500
        else:
            print('Create new table: {}'.format(filefullpath))
        gsm_maps = get_gsm_tables(self.gse)
        # step 1: Ready Affy files from folders
        gsm_tables_sc500 = {}
        for key in self.celformat.keys():
            platformdir = os.path.join(self.genedir, key)
            print('Affy Read Path: {}'.format(platformdir))
            if os.path.exists(platformdir):
                outputdf = affyio.readaffydir(platformdir)
            else:
                print('Path not exist: {}'.format(platformdir))
            outputdf['ENTREZ_GENE_ID'] = gsm_maps[key]['ENTREZ_GENE_ID']
            gsm_tables_sc500[key] = outputdf

        # step 2: Drop rows without ENTREZ GENE ID, set index to ENTREZ
        for key, val in self.celformat.items():
            gsm_tables_sc500[key].dropna(subset=['ENTREZ_GENE_ID'], inplace=True)
            gsm_tables_sc500[key].set_index('ENTREZ_GENE_ID', inplace=True)

        # step 3: Merge tables of platforms
        df_outer_sc500 = None
        for key, val in self.celformat.items():
            print('{}: {}'.format(key, gsm_tables_sc500[key].shape))
            if df_outer_sc500 is None:
                df_outer_sc500 = gsm_tables_sc500[key]
            else:
                df_outer_sc500 = pd.merge(df_outer_sc500, gsm_tables_sc500[key], on='ENTREZ_GENE_ID', how='outer')

        df_outer_sc500.dropna(how='all', inplace=True)
        print('Full: {}'.format(df_outer_sc500.shape))
        df_outer_sc500.rename(str.lower, axis='columns', inplace=True)

        # step 4: Remove duplicated items, keep largest VALUE for each GSM
        df_clean_sc500 = pd.DataFrame([], index=df_outer_sc500.index)
        df_clean_sc500 = df_clean_sc500[~df_clean_sc500.index.duplicated(keep='first')]
        for key, val in self.gsm_platform.items():
            key_low = key.lower()
            col1, col2, col3 = '{}.cel.gz'.format(key_low), '{}.cel.gz.1'.format(key_low), '{}.cel.gz.2'.format(key_low)
            try:
                temp = df_outer_sc500.loc[:, [col1, col2, col3]]
            except:
                print('{} not in df_outer_sc500'.format(key))
                continue
            temp.sort_values(by=['ENTREZ_GENE_ID', col1], inplace=True)
            temp = temp[~temp.index.duplicated(keep='last')]
            df_clean_sc500[col1] = temp[col1]
            df_clean_sc500[col2] = temp[col2]
            df_clean_sc500[col3] = temp[col3]

        # step 5: save to csv file
        df_clean_sc500.to_csv(filefullpath)
        df_clean_sc500.sort_index(inplace=True)
        print('Full table saved to:\n{}'.format(filefullpath))
        return df_clean_sc500

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