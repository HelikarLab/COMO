#!/usr/bin/python3
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
readaiglent <- function(addr){
    crd <- getwd()
    setwd(addr)
    x <- read.maimages(dir(".", "txt.gz"), "agilent", green.only = TRUE)
    setwd(crd)
    y <- backgroundCorrect(x,method="normexp") 
    y <- normalizeBetweenArrays(y,method="quantile") 
    return y
}
"""

agilentio = SignatureTranslatedAnonymousPackage(agilentstring, "agilentio")
