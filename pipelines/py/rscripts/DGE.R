zz <- file("/home/jupyteruser/work/py/rlogs/DGE.Rout", open="wt")
sink(zz, type="message")

library(DESeq2)
library(edgeR)
readCountMatrix <- function(cmat_file, config_file) {
  conf <- read.xlsx(config_file, sheet=1, colNames=TRUE)
  cmat_whole <- read.csv(cmat_file, header=TRUE)
  cmat_whole[,-1] <- lapply(cmat_whole[,-1], as.numeric)
  #print(cmat_whole)
  cmat_whole <- cmat_whole[rowSums(cmat_whole[,-1])>0,]
  genes <- cmat_whole$ENTREZ_GENE_ID
  samps <- as.character(conf$Sample)
  exps <- as.character(conf$Experiment)
  if ( length(genes) == 0 | is.null(genes) ) {
    print("disease count matrix must have a column headed 'ENTREZ_GENE_ID'")
  }
  SampMetrics <- list()
  for ( i in 1:length(samps) ) {
    entry <- samps[i]
    if ( entry %in% colnames(cmat_whole) ) {
      counts <- cmat_whole[entry]
      group <- exps[i]
      SampMetrics[[group]][[entry]][["Counts"]] <- counts
      SampMetrics[[group]][[entry]][["Entrez"]] <- genes
    } else if ( paste(c("X", entry),collapse="") %in% colnames(cmat_whole) ) {
      entry <- paste(c("X", entry),collapse="")
      counts <- cmat_whole[entry]
      group <- exps[i]
      SampMetrics[[group]][[entry]][["Counts"]] <- counts
      SampMetrics[[group]][[entry]][["Entrez"]] <- genes
    } else {
      print(paste(c(entry, " not found in disease count matrix")))
    }
  }
  return(SampMetrics)
}

dgeAnalysis <- function(SampMetrics, test_name, tissue_name, disease_name) {
  gene_list <- SampMetrics[[1]][[1]][["Entrez"]]
  
  df <- data.frame(Entrez=gene_list)
  #colnames(df)[1] <- "Entrez"
  group_list <- c(rep("control", length(SampMetrics[["control"]])), rep('patient', length(SampMetrics[["patient"]])))
  for ( j in 1:length(SampMetrics[["control"]]) ) {
    df <- cbind(df, SampMetrics[["control"]][[j]][["Counts"]])
  }
  for ( j in 1:length(SampMetrics[["patient"]]) ) {
    df <- cbind(df, SampMetrics[["patient"]][[j]][["Counts"]])
  }
  df[is.na(df)] <- 0
  df["Entrez"] <- NULL
  df <- data.frame(sapply(df, as.numeric))
  dgList <- DGEList(counts=df, genes=gene_list, group=group_list)
  dgList$samples$group <- relevel(dgList$samples$group, ref="control")
  dgList <- calcNormFactors(dgList, method="TMM")
  
  tmm <- cpm(dgList)
  write.csv(tmm, paste(c("/home/jupyteruser/work/data/results/", tissue_name, "/", disease_name, "/TMM_Matrix.csv"),collapse=""))
  
  # MDS Plot
  plotname <- paste(c("/home/jupyteruser/work/data/results/", tissue_name, "/", disease_name, "/", "MDS_plot.jpg"), collapse="")
  title <- paste(c("DGEList Multi-Dimensional Scaling"),collapse="")
  jpeg(plotname)
  lab <- colnames(df)
  plotMDS(dgList, labels=lab, main=title)
  dev.off()
  
  # create design matrix
  designMat <- model.matrix(~0+group, data=dgList$samples)
  colnames(designMat) <- levels(dgList$samples$group)
  
  # BCV plot
  plotname <- paste(c("/home/jupyteruser/work/data/results/", tissue_name, "/", disease_name, "/", "BCV_plot.jpg"), collapse="")
  title <- paste(c("DGEList Biological Coefficient of Variation"),collapse="")
  dgList <- estimateGLMCommonDisp(dgList, design=designMat)
  dgList <- estimateGLMTrendedDisp(dgList, design=designMat)
  dgList <- estimateGLMTagwiseDisp(dgList, design=designMat)
  jpeg(plotname)
  plotBCV(dgList, main=title)
  dev.off()
  
  # GLM approach
  fit <- glmQLFit(dgList, designMat)
  gp <- as.character(unique(dgList$samples$group))
  glen <- length(gp)
  for ( i in 2:glen ) {
    test_cell <- gp[i]
    con <- rep(0, glen)
    con[1] <- -1
    con[i] <- 1
    qlf <- glmQLFTest(fit, contrast=con)
    edgeR_result <- topTags(qlf, n=65000)
    deGenes <- decideTestsDGE(qlf, adjust.method="BH", p.value=0.05, lfc=1)
    deGenes <- rownames(qlf)[as.logical(deGenes)]
    expTab <- edgeR_result$table
    
    # save results
    #expTab['abs_logFC'] <- abs(expTab$logFC)
    names(expTab)[names(expTab) == "PValue"] <- "P.Value"
    names(expTab)[names(expTab) == "genes"] <- "Gene ID"
    #expTab <- expTab[,deGenes==TRUE]
    write.csv(expTab,paste(c("/home/jupyteruser/work/data/results/", tissue_name, "/", disease_name, "/", "DiffExp.csv"),collapse=""))
    
    # smear plot
    plotname <- paste(c("/home/jupyteruser/work/data/results/", tissue_name, "/", disease_name, "/", "smear_plot.jpg"),collapse="")
    title <- paste(c("DGEList Smear Plot ", test_cell),collapse="")
    jpeg(plotname)
    plotSmear(qlf, de.tags=deGenes, main=title)
    abline(h=c(-1, 1), col=2)
    dev.off()
  }
  return(expTab)
}

DGE_main <- function(cmat_file, config_file, tissue_name, disease_name) {
  print("Reading Counts Matrix")
  test_name <- cmat_file
  test_name <-unlist(strsplit(test_name, "_RawCounts"))[1]
  test_list <-unlist(strsplit(test_name, "/"))
  test_name <- test_list[length(test_list)]
  SampMetrics <- readCountMatrix(cmat_file, config_file)
  entrez_all <- SampMetrics[[1]][[1]][["Entrez"]]
  print("Performing DGE")
  data_table <- dgeAnalysis(SampMetrics, test_name, tissue_name, disease_name)
  return(data_table)
}