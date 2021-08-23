zz <- file("/home/jupyteruser/work/py/rlogs/bulk.Rout", open="wt")
sink(zz, type="message")

library(tidyverse)
library(limma)
library(edgeR)
library(genefilter)
library(biomaRt)
library(sjmisc)


readCountMatrix <- function(cmat_file, config_file, info_file, model_name) {
  gene_info <- read.csv(info_file)
  gene_info$size <- (gene_info$end_position-gene_info$start_position)
  conf <- read.xlsx(config_file, model_name, colNames=TRUE)
  cmat_whole <- read.csv(cmat_file, header=TRUE)
  cmat_whole <- cmat_whole[(gene_info$entrezgene_id!="-"),]
  gene_info <- gene_info[(gene_info$entrezgene_id!="-"),]
  genes <- gene_info$entrezgene_id
  # remove version numbers from ensembl id
  for ( j in 1:length(genes) ) {
    r <- genes[j]
    if ( grepl("\\\\.", r) ) {
      gen <- unlist(str_split(r, "\\\\."))[1]
      genes[j] <- gen
    }
  }
  samp_mat <- genes
  SampMetrics <- list()
  insert_sizes <- c()
  samp_names <- c()
  groups <- unique(conf$Group)
  # initialize groups
  for (group in groups) {
    SampMetrics[[group]][["CountMatrix"]] <- genes
    SampMetrics[[group]][["InsertSizes"]] <- c()
    SampMetrics[[group]][["SampNames"]] <- c()
  }
  # add to group count matrices and insert lists
  for ( i in 1:length(conf$SampleName) ) {
    entry <- conf$SampleName[i]
    #group <- unlist(str_split(entry, "_", n=2))[2]
    group <- conf$Group[i]
    insert_size <- conf$InsertSize[i]     
    if ( entry %in% colnames(cmat_whole) ) {
        samp_mat <- SampMetrics[[group]][["CountMatrix"]]
        insert_sizes <- SampMetrics[[group]][["InsertSizes"]]
        samp_names <- SampMetrics[[group]][["SampNames"]]
        samp_mat <- cbind(samp_mat, cmat_whole[,entry])
        insert_sizes <- c(insert_sizes, conf$InsertSize[i])
        samp_names <- c(samp_names, entry)
        samp_mat <- SampMetrics[[group]][["CountMatrix"]] <- samp_mat
        SampMetrics[[group]][["InsertSizes"]] <- insert_sizes
        SampMetrics[[group]][["SampNames"]] <-  samp_names
    } else {
       print(paste(c(entry., " not found in count matrix."),collapse=""))
    }
  }  
  for (group in groups) {
    samp_mat <- SampMetrics[[group]][["CountMatrix"]] 
    #print(samp_mat)
    samp_mat <- samp_mat[,-1]
    samp_mat <- sapply(data.frame(samp_mat), as.numeric)
    samp_mat <- as.matrix(samp_mat)
    colnames(samp_mat) <- SampMetrics[[group]][["SampNames"]]
    SampMetrics[[group]][["CountMatrix"]] <- samp_mat
    SampMetrics[[group]][["NumSamples"]] <- ncol(samp_mat)
    SampMetrics[[group]][["InsertSizes"]] <- insert_sizes
    SampMetrics[[group]][["Entrez"]] <- as.character(gene_info$entrezgene_id)
    SampMetrics[[group]][["GeneSizes"]] <- gene_info$size
  }
  return(SampMetrics)
}


calculateTPM <- function(SampMetrics, cell_types) {
  for ( i in 1:length(SampMetrics) ) {
    count_matrix <- SampMetrics[[i]][["CountMatrix"]]
    gene_size <- SampMetrics[[i]][["GeneSizes"]]
    tpm_matrix <- do.call(cbind, lapply(1:ncol(count_matrix), function(j) {
      rate = log(count_matrix[,j]) - log(gene_size[j])
      denom = log(sum(exp(rate)))
      exp(rate - denom + log(1e6))
    }))
    colnames(tpm_matrix) <- colnames(count_matrix)
    SampMetrics[[i]][["TPM_Matrix"]] <- tpm_matrix
  }
  return(SampMetrics)
}


calculateZscore <- function(SampMetrics, cell_types, norm_tech) {
  for ( i in 1:length(SampMetrics)) {
    if ( norm_tech=="CPM" ) {
      tmat <- SampMetrics[[i]][["CPM_Matrix"]]
    } else if ( norm_tech=="TPM" ) {
      tmat <- SampMetrics[[i]][["TPM_Matrix"]]
    }
    zmat <- matrix(nrow=nrow(tmat), ncol=ncol(tmat))
    #rownames(zmat) <- rownames(tmat)
    for ( j in 1:ncol(tmat) ) {
      tvec <- tmat[,j]
      logvec <- log2(tvec)
      logvec[is.infinite(logvec)] <- NA
      zvec <- scale(logvec, center=TRUE, scale=TRUE)
      zmat[,j] <- zvec
    }
    zmat <- data.frame(zmat)
    colnames(zmat) <- colnames(tmat)
    SampMetrics[[i]][["Zscore"]] <- zmat
  }
  return(SampMetrics)
}


CPM_filter <- function(SampMetrics, filt_options, cell_types) {
  N_exp <- filt_options$pos_rep
  N_top <- filt_option$top_rep
  min.count <- filt_options$min_count
  for ( i in 1:length(SampMetrics) ) {
    counts <- SampMetrics[[i]][["CountMatrix"]]
    ent <- SampMetrics[[i]][["Entrez"]]
    size <- SampMetrics[[i]][["GeneSizes"]]
    lib.size <- colSums(counts)
    MedianLibSize <- median(lib.size)
    if ( min.count=="default" ) {
      CPM.Cutoff <- 10000000/(median(colSums(counts)))
    } else {
      CPM.Cutoff <- min.count/MedianLibSize*1e6
    }
    CPM <- cpm(counts,lib.size=lib.size)
    min.samples <- round(N_exp * ncol(counts))
    top.samples <- round(N_top * ncol(counts))
    test_bools_top <- data.frame(genes=ent)
    for ( j in 1:ncol(CPM) ) {
      cpm_q <- CPM[,j]
      cpm_q <- cpm_q[cpm_q>0]
      q_cutoff_top <- quantile(cpm_q, prob=1-perc_top/100)
      test_bools_top <- cbind(test_bools_top, as.integer(CPM[,j]>q_cutoff_top))
    }
    
    f1 <- genefilter::kOverA(min.samples, CPM.Cutoff)
    flist <- genefilter::filterfun(f1)
    keep <- genefilter::genefilter(CPM, flist)
    SampMetrics[[i]][["Entrez"]] <- ent[keep]
    SampMetrics[[i]][["GeneSizes"]] <- size[keep]
    SampMetrics[[i]][["CountMatrix"]] <- counts[keep,]
    SampMetrics[[i]][["CPM_Matrix"]] <- CPM[keep,]
    
    # top percentile genes
    test_bools_top["genes"] <- NULL
    f1_top <- genefilter::kOverA(top.samples, CPM.Cutoff)
    flist_top <- genefilter::filterfun(f1_top)
    keep_top <- genefilter::genefilter(test_bools_top, flist_top)
    SampMetrics[[i]][["Entrez_top"]] <- ent[keep_top]
  }
  SampMetrics <- calculateZscore(SampMetrics, cell_types, "CPM")
  return(SampMetrics)
}


TPM_quant_filter <- function(SampMetrics, filt_options, cell_types) {
  N_exp <- filt_options$pos_rep
  N_top <- filt_options$top_rep
  perc <- filt_options$percentile
  perc_top <- filt_options$top_percentile
  SampMetrics <- calculateTPM(SampMetrics, cell_types)
  for ( i in 1:length(SampMetrics) ) {
    
    counts <- SampMetrics[[i]][["CountMatrix"]]
    ent <- SampMetrics[[i]][["Entrez"]]
    size <- SampMetrics[[i]][["GeneSizes"]]
    tpm <- SampMetrics[[i]][["TPM_Matrix"]]
    min.samples <- round(N_exp * ncol(tpm))
    top.samples <- round(N_top * ncol(tpm))
    test_bools <- data.frame(gene=ent)
    test_bools_top <- test_bools
    for ( j in 1:ncol(tpm) ) {
      tpm_q <- tpm[,j]
      tpm_q <- tpm_q[tpm_q>0]
      q_cutoff <- quantile(tpm_q, prob=1-perc/100)
      #q_cutoff_top <- quantile(tpm_q, prob=1-perc_top/100)
      #bools <- data.frame(as.integer(tpm[,j]>q_cutoff))
      #bools_top <- data.frame(as.integer(tpm[,j]>q_cutoff_top))
      test_bools <- cbind(test_bools, as.integer(tpm[,j]>q_cutoff))
    }
    test_bools["gene"] <- NULL
    #test_bools_top["gene"] <- NULL
    f1 <- genefilter::kOverA(min.samples, 0.9)
    flist <- genefilter::filterfun(f1)
    keep <- genefilter::genefilter(test_bools, flist)
    SampMetrics[[i]][["Entrez"]] <- ent[keep]
    SampMetrics[[i]][["GeneSizes"]] <- size[keep]
    SampMetrics[[i]][["CountMatrix"]] <- counts[keep,]
    SampMetrics[[i]][["TPM_Matrix"]] <- tpm[keep,]
    f1_top <- genefilter::kOverA(top.samples, 0.9)
    flist_top <- genefilter::filterfun(f1_top)
    keep_top <- genefilter::genefilter(test_bools, flist_top)
    SampMetrics[[i]][["Entrez_top"]] <- ent[keep_top]
  }
  SampMetrics <- calculateZscore(SampMetrics, cell_types, "TPM")
  return(SampMetrics)
}


zFPKM_filter <- function(SampMetrics, filt_options, cell_types) {
  N_exp <- filt_options$pos_rep
  N_top <- filt_options$top_rep
  perc_top <- filt_options$top_percentile
  for ( i in 1:length(SampMetrics) ) {
    cmat <- SampMetrics[[i]][["CountMatrix"]]
    #rownames(cmat) <- make.names(gnames, unique = TRUE)
    fmat <- fpkm(cmat, SampMetrics[[i]][["GeneSizes"]],
                 SampMetrics[[i]][["InsertSizes"]])
    
    ent <- SampMetrics[[i]][["Entrez"]]
    size <- SampMetrics[[i]][["GenesSizes"]]
    
    fdf <- data.frame(fmat)
    zmat <- zFPKM(fdf, assayName="FPKM")
    #inames <- rownames(zmat)
    
    min.samples <- round(N_exp * ncol(zmat))
    top.samples <- round(N_top * ncol(zmat))
    test_bools_top <- data.frame(genes=ent)
    for ( j in 1:ncol(zmat) ) {
      z_q <- zmat[,j]
      z_q <- z_q[z_q>-3]
      q_cutoff_top <- quantile(z_q, prob=1-perc_top/100)
      test_bools_top <- cbind(test_bools_top, as.integer(zmat[,j]>q_cutoff_top))
    }
    
    cutoff <- -3
    f1 <- genefilter::kOverA(min.samples, cutoff)
    flist <- genefilter::filterfun(f1)
    keep <- genefilter::genefilter(zmat, flist)
    SampMetrics[[i]][["zFPKM_Matrix"]] <- zmat[keep,]
    SampMetrics[[i]][["FPKM_Matrix"]] <- fmat
    
    # top percentile genes
    test_bools_top["genes"] <- NULL
    f1_top <- genefilter::kOverA(top.samples, cutoff)
    flist_top <- genefilter::filterfun(f1_top)
    keep_top <- genefilter::genefilter(test_bools_top, flist_top)
    SampMetrics[[i]][["Entrez_top"]] <- ent[keep_top]
  }
  return(SampMetrics)
}


filterCounts <- function(SampMetrics, technique, filt_options, mart, cell_types) {
  switch(technique,
         cpm = CPM_filter(SampMetrics, filt_options, cell_types),
         zFPKM = zFPKM_filter(SampMetrics, filt_options, mart, cell_types),
         quantile = TPM_quant_filter(SampMetrics, filt_options, cell_types))
}



save_bulk_tests <- function(cmat_file, config_file, out_file, info_file, model_name,
                            pos_rep=0.5, pos_samp=0.5, top_rep=0.9, top_samp=0.9, 
                            technique="quantile", quantile=0.9, min_count=10) {
  
  # condense filter options
  filt_options <- list()
  if ( exists("pos_rep") ) {
    filt_options$pos_rep <- pos_rep
  } else {
    filt_options$pos_rep <- 0.2
  }
  if ( exists("pos_samp") ) {
    filt_options$pos_samp <- pos_samp
  } else {
    filt_options$pos_samp <- 0.5
  }
  if ( exists("quantile") ) {
    filt_options$percentile <- quantile
  } else {
    filt_options$percentile <- 50
  }
  if ( exists("min_count") ) {
    filt_options$min_count <- min_count
  } else {
    filt_options$min_count <- 1
  }
  if ( exists("top_rep") ) {
    filt_options$top_rep<- top_rep
  } else {
    filt_options$top_rep <- 0.9
  }
  if ( exists("top_samp") ) {
    filt_options$top_samp <- top_samp
  } else {
    filt_options$top_samp <- 0.9
  }
  print("Reading Counts Matrix")
  SampMetrics <- readCountMatrix(cmat_file, config_file, info_file, model_name)
  entrez_all <- SampMetrics[[1]][["Entrez"]]
  print("Filtering Counts")
  SampMetrics <- filterCounts(SampMetrics, technique, filt_options, mart, cell_types)
  expressedGenes <- c()
  topGenes <- c()
  for ( i in 1:length(SampMetrics) ) {
    expressedGenes <- c(expressedGenes, SampMetrics[[i]][["Entrez"]])
    topGenes <- c(topGenes, SampMetrics[[i]][["Entrez_top"]])
  }
  expMat <- as.data.frame(table(expressedGenes))
  topMat <- as.data.frame(table(topGenes))
  nc <- length(SampMetrics)
  expMat <- cbind(expMat, "Prop"=expMat$Freq/nc)
  topMat <- cbind(topMat, "Prop"=topMat$Freq/nc)
  SampMetrics[["ExpressionMatrix"]] <- expMat
  SampMetrics[["TopMatrix"]] <- topMat
  SampMetrics[["ExpressedGenes"]] <- as.character(expMat$expressedGenes[expMat$Prop>=pos_samp])
  SampMetrics[["TopGenes"]] <- as.character(topMat$topGenes[topMat$Prop>=top_samp])
  #write_table <- data.frame(entrez_all[!is.na(entrez_all)])
  write_table <- data.frame(entrez_all)
  write_table <- cbind(write_table, rep(0, nrow(write_table)))
  write_table <- cbind(write_table, rep(0, nrow(write_table)))
  for ( i in 1:nrow(write_table) ) {
    if (as.character(write_table[i,1]) %in% as.character(SampMetrics$ExpressedGenes)) {
      write_table[i,2] <- 1
    }
    if (as.character(write_table$entrez_all[i]) %in% as.character(SampMetrics$TopGenes)) {
      write_table[i,3] <- 1
    }
  }
  header <- c("ENTREZ_GENE_ID", "expressed", "top")
  #write_table <- rbind(header, write_table)
  colnames(write_table) <- header
  write.csv(write_table, out_file, row.names=FALSE, col.names=FALSE)
}