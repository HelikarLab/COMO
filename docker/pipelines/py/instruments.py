#!/usr/bin/python3
import os, time
import pandas as pd
from rpy2.robjects.packages import importr
from rpy2.robjects import r, pandas2ri
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage
from bioservices import BioDBNet



# automatically convert ryp2 dataframe to Pandas dataframe
pandas2ri.activate()
# Initialize Rpy2 for Affy package
affy = importr("affy")
limma = importr("limma")
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

fitaffydir <- function(addr, target){
    crd <- getwd()
    setwd(addr)
    library(limma)
    library(affy)
    targets = readTargets(target)
    mydata = ReadAffy()
    setwd(crd)
    eset = rma(mydata)
    f <- factor(targets$Condition, levels = unique(targets$Condition))
    design <- model.matrix(~0 + f)
    colnames(design) <- levels(f)
    fit = lmFit(eset,design)
    contrast.matrix = makeContrasts("patient-control",levels = design)
    fit2 = contrasts.fit(fit,contrast.matrix)
    fit2= eBayes(fit2)
    data = topTable(fit2,number = "inf")
    # write.table(data,"differentialExpression.txt",sep = "\t")
    return(data)
}
"""
affyio = SignatureTranslatedAnonymousPackage(string, "affyio")


agilentstring="""
readaiglent <- function(addr,targets){
    crd <- getwd()
    setwd(addr)
    # targets <- dir(".", "txt.gz")
    x <- read.maimages(targets, "agilent", green.only = TRUE)
    setwd(crd)
    y <- backgroundCorrect(x,method="normexp") 
    y <- normalizeBetweenArrays(y,method="quantile") 
    yo <- c(y$genes,as.data.frame(y$E))
    ydf <- data.frame(yo)
    return(ydf)
}
fitagilent <- function(addr,target){
    crd <- getwd()
    setwd(addr)    
    targets <- readTargets(target)
    x <- read.maimages(targets, path="somedirectory", source="agilent",green.only=TRUE)
    y <- backgroundCorrect(x, method="normexp", offset=16)
    y <- normalizeBetweenArrays(y, method="quantile")
    y.ave <- avereps(y, ID=y$genes$ProbeName)
    f <- factor(targets$Condition, levels = unique(targets$Condition))
    design <- model.matrix(~0 + f)
    colnames(design) <- levels(f)
    fit <- lmFit(y.ave, design)
    contrast.matrix <- makeContrasts("Treatment1-Treatment2", "Treatment1-Treatment3", "Treatment2-Treatment1", levels=design)
    fit2 = contrasts.fit(fit,contrast.matrix)
    fit2 = eBayes(fit2)
    data = topTable(fit2,number = "inf")
    return(data)
    }
"""

agilentio = SignatureTranslatedAnonymousPackage(agilentstring, "agilentio")

def agilent_raw(datadir, gsms):
    files = os.listdir(datadir)
    txts = []
    gzs = []
    keys = []
    if gsms:
        for gsm in gsms:
            for file in files:
                if gsm in file:
                    gzs.append(file)
                    txts.append(file[0:-3])
                    keys.append(gsm)
    else:
        for file in files:
            gzs.append(file)
            txts.append(file[0:-3])
            keys.append(file.split('_')[0])
    cols = dict(zip(txts,keys))

    targets = pd.DataFrame(gzs,columns=['FileName'],index=txts)
    df_agilent = agilentio.readaiglent(datadir, targets)
    return df_agilent.rename(columns=cols)

def readagilent(datadir, gsms, scalefactor=1.1, quantile=0.95):
    df_raw = agilent_raw(datadir,gsms)
    df = df_raw.drop(columns=['Row', 'Col'])
    df_negctl = df[df['ControlType'] == -1]
    df_cutoff = scalefactor * df_negctl.drop(columns=['ControlType']).quantile(quantile, axis=0)
    df_bool = df.loc[df['ControlType'] == 0, df_cutoff.index.tolist()].gt(df_cutoff, axis=1)
    idx_ones = df_bool[df_bool.all(axis=1)].index
    idx_zeros = df_bool[~df_bool.all(axis=1)].index
    df.loc[idx_ones, 'Express'] = 1
    df.loc[idx_zeros, 'Express'] = 0
    df_results = df.loc[df['ControlType'] == 0, :].copy()
    for gsm in gsms:
        col = '{}.cel.gz.1'.format(gsm.lower())
        df_results.loc[:, col] = 'A'
        df_results.loc[df_bool.loc[:, gsm],col] = 'P'
        col = '{}.cel.gz.2'.format(gsm.lower())
        df_results.loc[:, col] = 1.0-quantile
        col = '{}.cel.gz'.format(gsm.lower())
        df_results.rename(columns={gsm: col}, inplace=True)
    df_results.rename(columns={'ProbeName': 'ID'}, inplace=True)
    df_results.set_index('ID', inplace=True)
    return df_results.drop(['ControlType','SystematicName'],axis=1)


def fetch_entrez_gene_id(input_values, input_db='Agilent ID', output_db = ['Gene ID','Ensembl Gene ID'], delay=30):
    s = BioDBNet()
    # input_db = 'Agilent ID'
    # input_values = df_results.ProbeName.tolist()

    df_maps = pd.DataFrame([],columns=output_db)
    df_maps.index.name=input_db
    i = 0
    # for i in range(0,len(input_values),500):
    while i < len(input_values):
        print('retrieve {}:{}'.format(i,min(i+500,len(input_values))))
        df_test = s.db2db(input_db, output_db, input_values[i:min(i+500,len(input_values))], 9606)
        if isinstance(df_test, pd.DataFrame):
            df_maps = pd.concat([df_maps, df_test], sort=False)
        elif df_test == '414':
            print("bioDBnet busy, try again in {} seconds".format(delay))
            time.sleep(delay)
            continue
        i += 500
    return df_maps