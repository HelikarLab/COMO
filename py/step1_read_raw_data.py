import GEOparse
import re
import os
import pandas as pd
from rpy2.robjects.packages import importr
from rpy2.robjects import r, pandas2ri
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage

pandas2ri.activate()
# Input: Extract Gene Info from GEO DataSets

gsename = "GSE2770"
softfile = "./{}_family.soft.gz".format(gsename)
if os.path.isfile(softfile):
    gse = GEOparse.get_GEO(filepath=softfile)
else:
    gse = GEOparse.get_GEO(geo=gsename, destdir="./")

# Extract Platform Information
keys = []
values = []
for gpl in gse.gpls:
    keys.append(gse.gpls[gpl].name)
    r1 = re.findall(r"\[.*?\]", gse.gpls[gpl].metadata['title'][0])
    values.append(r1[0][4:-1])
    # values.append(gse.gpls[gpl].metadata['title'][0])
    print(gse.gpls[gpl].name)
    print(gse.gpls[gpl].metadata['title'][0])

celformat = dict(zip(keys, values))

## Classify the Samples by Platform
keys = []
values = []
for gsm in gse.gsms:
    #print(gse.gsms[gsm].name)
    keys.append(gse.gsms[gsm].name)
    for key,val in celformat.items():
        r1 = re.findall(r"{}".format(val),gse.gsms[gsm].columns.loc['ID_REF','description'])
        if not r1:
            pass
        else:
            values.append(key)
    #r1 = re.findall(r"\(.*?\)",gse.gsms[gsm].columns.loc['ID_REF','description'])
    #values.append(r1[0][1:-1])
    #values.append(gse.gsms[gsm].columns.loc['ID_REF','description'])

gsm_platform = dict(zip(keys, values))

## Setup paths
currentdir = os.getcwd()
dirlist = currentdir.split('/')
projectdir = '/'.join(dirlist[0:-1])
datadir = os.path.join(projectdir,'data')
outputdir = os.path.join(projectdir,'output')
gene = 'GSE2770'
genedir = os.path.join(datadir,gene + '_RAW')

# create a folder for each platform
for key in celformat.keys():
    platformdir = os.path.join(genedir,key)
    if not os.path.exists(platformdir):
        os.makedirs(platformdir)
        print('Path created: {}'.format(platformdir))
    else:
        print('Path already exist: {}'.format(platformdir))

# Move Corresponding Cel files to Folders
#onlyfiles = [f for f in os.listdir(genedir) if os.path.isfile(os.path.join(genedir, f))]
onlyfiles = [f for f in os.listdir(genedir) if f.endswith('.gz')]

for file in onlyfiles:
    filelist = file.split('.')
    prefix = filelist[0]
    if prefix in gsm_platform:
        platform = gsm_platform[prefix]
        platformdir = os.path.join(genedir,platform)
        src_path = os.path.join(genedir, file)
        dst_path = os.path.join(platformdir, file)
        os.rename(src_path,dst_path)
        print('Move {} to {}'.format(src_path,dst_path))

platforms = ['GPL96','GPL97','GPL8300']

maps_list = []
gene_maps = pd.DataFrame([],columns=['GPL96','GPL97','GPL8300','ENTREZ_GENE_ID'])
gene_maps.set_index('ENTREZ_GENE_ID',inplace=True)
for platform in platforms:
    temp =gse.gpls[platform].table[['ID','ENTREZ_GENE_ID']]
    # Save to file
    filefullpath = os.path.join(datadir,'{}entrez.csv'.format(platform))
    print(filefullpath)
    temp.to_csv(filefullpath, index=False)
    # Single Table
    temp.dropna(axis=0,inplace=True)
    temp.set_index('ENTREZ_GENE_ID',inplace=True)
    maps_list.append(temp)


# create 3 tables by platform

# initialize with platform
gsm_tables = {}
for key,val in celformat.items():
    #df = pd.DataFrame([],columns=['ID','ENTREZ_GENE_ID'])
    temp = gse.gpls[key].table
    df = temp[['ID','ENTREZ_GENE_ID']]
    df.set_index('ID',inplace=True)
    df.drop_duplicates(keep='last',inplace=True)
    gsm_tables[key] = df

# fill each table
for key,val in gsm_platform.items():
    temp = gse.gsms[key].table.copy()
    temp.rename(columns={"ID_REF": "ID"},inplace=True)
    temp.dropna(subset=['ID'],how='any',inplace=True)
    temp.set_index('ID',inplace=True)
    #col1,col2,col3 = '{}.VALUE'.format(key), '{}.ABS_CALL'.format(key), '{}.DETECTION P-VALUE'.format(key)
    col1,col2,col3 = '{}.CEL.gz'.format(key), '{}.CEL.gz.1'.format(key), '{}.CEL.gz.2'.format(key)
    #gsm_tables[val].loc[:,[col1,col2,col3]] = temp[['VALUE','ABS_CALL','DETECTION P-VALUE']]
    # Sort by VALUE and drop duplicated ID
    temp['ENTREZ_GENE_ID']=gsm_tables[val]['ENTREZ_GENE_ID']
    temp.sort_values(by=['ENTREZ_GENE_ID', 'VALUE'],inplace=True)
    temp.drop_duplicates(subset=['ENTREZ_GENE_ID'],keep='last',inplace=True)
    gsm_tables[val][col1] = temp['VALUE']
    gsm_tables[val][col2] = temp['ABS_CALL']
    gsm_tables[val][col3] = temp['DETECTION P-VALUE']

# dropna
for key,val in celformat.items():
    #df = pd.DataFrame([],columns=['ID','ENTREZ_GENE_ID'])
    gsm_tables[key].dropna(subset=['ENTREZ_GENE_ID'],inplace=True)
    gsm_tables[key].set_index('ENTREZ_GENE_ID',inplace=True)
    # gsm_tables[key].dropna(how='all',inplace=True)
    # gsm_tables[key].drop_duplicates(keep='last',inplace=True)

# merge files
df_outer = None
for key,val in celformat.items():
    #df = pd.DataFrame([],columns=['ID','ENTREZ_GENE_ID'])
    print('{}: {}'.format(key,gsm_tables[key].shape))
    if df_outer is None:
        df_outer = gsm_tables[key]
    else:
        df_outer = pd.merge(df_outer,gsm_tables[key],on='ENTREZ_GENE_ID', how='outer')

df_outer.dropna(how='all',inplace=True)
print('full : {}'.format(df_outer.shape))
df_outer.sort_index(inplace=True)
df_outer.to_csv(os.path.join(genedir,'{}_full_table.csv'.format(gsename)))


affy = importr("affy")
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
# initialize with platform
gsm_maps = {}
for key,val in celformat.items():
    #df = pd.DataFrame([],columns=['ID','ENTREZ_GENE_ID'])
    temp = gse.gpls[key].table.copy()
    df = temp[['ID','ENTREZ_GENE_ID']]
    df.set_index('ID',inplace=True)
    #df.drop_duplicates(keep='last',inplace=True)
    gsm_maps[key] = df

gsm_tables_sc500 = {}
for key in celformat.keys():
    platformdir = os.path.join(genedir,key)
    print('Affy Read Path: {}'.format(platformdir))
    if os.path.exists(platformdir):
        outputdf = affyio.readaffydir(platformdir)
    else:
        print('Path not exist: {}'.format(platformdir))
    outputdf['ENTREZ_GENE_ID'] = gsm_maps[key]['ENTREZ_GENE_ID']
    gsm_tables_sc500[key] = outputdf

# Read files
# gpl8300rawdir = '/Users/zhzhao/Dropbox/Helikar/pipelines/data/GSE2770_RAW/GPL8300/'
# outputdf = affyio.readaffydir(gpl8300rawdir)
for key,val in celformat.items():
    #df = pd.DataFrame([],columns=['ID','ENTREZ_GENE_ID'])
    try:
        gsm_tables_sc500[key].dropna(subset=['ENTREZ_GENE_ID'],inplace=True)
        gsm_tables_sc500[key].set_index('ENTREZ_GENE_ID',inplace=True)
    except:
        pass

# Merge table to one
df_outer_sc500 = None
for key,val in celformat.items():
    #df = pd.DataFrame([],columns=['ID','ENTREZ_GENE_ID'])
    print('{}: {}'.format(key,gsm_tables[key].shape))
    if df_outer_sc500 is None:
        df_outer_sc500 = gsm_tables_sc500[key]
    else:
        df_outer_sc500 = pd.merge(df_outer_sc500,gsm_tables_sc500[key],on='ENTREZ_GENE_ID', how='outer')

df_outer_sc500.dropna(how='all',inplace=True)
print('Full: {}'.format(df_outer.shape))
df_outer_sc500.rename(str.lower, axis='columns',inplace=True)


df_clean_sc500 = pd.DataFrame([],index=df_outer_sc500.index)
df_clean_sc500 = df_clean_sc500[~df_clean_sc500.index.duplicated(keep='first')]
# df_outer_sc500.rename(str.lower, axis='columns',inplace=True)
for key,val in gsm_platform.items():
    key_low = key.lower()
    col1,col2,col3 = '{}.cel.gz'.format(key_low), '{}.cel.gz.1'.format(key_low), '{}.cel.gz.2'.format(key_low)
    try:
        temp = df_outer_sc500.loc[:,[col1,col2,col3]]
    except:
        print('{} not in df_outer_sc500'.format(key))
        continue
    temp.sort_values(by=['ENTREZ_GENE_ID',col1],inplace=True)
    temp = temp[~temp.index.duplicated(keep='last')]
    df_clean_sc500[col1] = temp[col1]
    df_clean_sc500[col2] = temp[col2]
    df_clean_sc500[col3] = temp[col3]

df_clean_sc500.to_csv(os.path.join(genedir,'{}_sc500_full_table.csv'.format(gsename)))
