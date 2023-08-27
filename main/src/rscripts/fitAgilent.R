suppressPackageStartupMessages(library("stringr"))

# Check if rlogs directory exists, From: https://stackoverflow.com/a/46008094
# Then prevent messy messages from repeatedly writing to juypter
work_dir <- getwd()
r_log_directory <- str_interp("${work_dir}/logs")
if (!dir.exists(r_log_directory)) { dir.create(r_log_directory) }
zz <- file(file.path(r_log_directory, "fitAligent.Rout"), open="wt")
sink(zz, type="message")

readagilent <- function(addr,targets){
    crd <- getwd()
    setwd(addr)
    # targets <- dir(".", "txt.gz")
    x <- read.maimages(targets, "agilent", green.only = TRUE)
    setwd(crd)
    y <- backgroundCorrect(x,method="normexp") 
    y <- normalizeBetweenArrays(y,method="quantile") 
    yo <- c(y$genes,as.data.frame(y$E))
    ydf <- data.frame(yo)
    write.csv(ydf, paste(c(addr, "/ydf_temp.csv"),collapse=""))
    
    return(ydf)
}
fitagilent <- function(addr,target){
    username <- Sys.info()["user"]
    work_dir <- str_interp("/home/${username}/work")
    crd <- getwd()
    setwd(addr)
    targets <- readTargets(target)
    x <- read.maimages(dir(), path = addr, source="agilent",green.only=TRUE)
    y <- backgroundCorrect(x, method="normexp", offset=16)
    y <- normalizeBetweenArrays(y, method="quantile")
    y.ave <- avereps(y, ID=y$genes$ProbeName)
    f <- factor(targets$Condition, levels = unique(targets$Condition))
    design <- model.matrix(~0 + f)
    colnames(design) <- levels(f)
    fit <- lmFit(y.ave, design)
    contrast.matrix <- makeContrasts("control-patient", levels=design)
    fit2 = contrasts.fit(fit,contrast.matrix)
    fit2 = eBayes(fit2)
    data = topTable(fit2,number = "inf")
    
    return(data)
    }
