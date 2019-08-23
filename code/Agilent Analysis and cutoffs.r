library(limma)
setwd('/Users/zhzhao/Dropbox/Helikar/pipelines/data/GSE43005_RAW')
x <- read.maimages(dir(".", "txt.gz"), "agilent", green.only = TRUE)
y <- backgroundCorrect(x,method="normexp") 
y <- normalizeBetweenArrays(y,method="quantile") 

##Now filter out control probes and low expressed probes. To get an idea of how bright expression probes should be, we compute the 95% percentile of the negative control probes on each array. We keep probes that are at least 10% brighter than the negative controls on at least four arrays (because there are four replicates) #####
neg95 <- apply(y$E[y$genes$ControlType==-1,],2,function(x) quantile(x,p=0.95)) 
cutoff <- matrix(1.1*neg95,nrow(y),ncol(y),byrow=TRUE) 
isexpr <- rowSums(y$E > cutoff) >= 4    # (this number should be the replicates)
table(isexpr) 
y0 <- y[y$genes$ControlType==0 & isexpr,] 

y1 <- y[y$genes$ControlType==0 & !isexpr,] 
