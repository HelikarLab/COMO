#!/usr/bin/python3
import os
import pandas as pd
from rpy2.robjects.packages import importr
from rpy2.robjects import r, pandas2ri
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage



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
"""

agilentio = SignatureTranslatedAnonymousPackage(agilentstring, "agilentio")

def agilent_raw(datadir, gsms):
    files = os.listdir(datadir)
    txts = []
    gzs = []
    if gsms:
        for gsm in gsms:
            for file in files:
                if gsm in file:
                    gzs.append(file)
                    txts.append(file[0:-3])
    else:
        for file in files:
            gzs.append(file)
            txts.append(file[0:-3])

    targets = pd.DataFrame(gzs,columns=['FileName'],index=txts)
    df_agilent = agilentio.readaiglent(datadir, targets)
    return df_agilent

def readagilent(datadir,gsms,scalefactor=1.1):
    df_raw = agilent_raw(datadir,gsms)
    df = df_raw.drop(columns=['Row', 'Col'])
    df_negctl = df[df['ControlType'] == -1]
    df_cutoff = scalefactor * df_negctl.drop(columns=['ControlType']).quantile(0.95, axis=0)
    df_bool = df.loc[df['ControlType'] == 0, df_cutoff.index.tolist()].gt(df_cutoff, axis=1)
    idx_ones = df_bool[df_bool.all(axis=1)].index
    idx_zeros = df_bool[~df_bool.all(axis=1)].index
    df.loc[idx_ones, 'Express'] = 1
    df.loc[idx_zeros, 'Express'] = 0
    df_results = df.loc[df['ControlType'] == 0, ['ProbeName','SystematicName','Express']]
    return df_results
