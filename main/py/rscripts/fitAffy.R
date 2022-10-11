# Check if rlogs directory exists, From: https://stackoverflow.com/a/46008094
library("stringr")
username <- Sys.info()["user"]
work_dir <- str_interp("/home/${username}/work")

if (!dir.exists(str_interp("${work_dir}/py/rlogs"))) {
    dir.create(str_interp("${work_dir}/py/rlogs"))
}

zz <- file(file.path("/home", username, "main", "py", "rlogs", "fitAffy.Rout"), open="wt")
sink(zz, type="message")

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
username <- Sys.info()["user"]
work_dir <- str_interp("/home/${username}/work")
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
