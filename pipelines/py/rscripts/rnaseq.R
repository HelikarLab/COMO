# prevent messy messages from repeatedly writing to juypter
zz <- file("/home/jupyteruser/work/py/rlogs/rnaseq.Rout", open="wt")
sink(zz, type="message")

library(tidyverse)
library(limma)
library(edgeR)
library(genefilter)
library(biomaRt)
library(sjmisc)
library(zFPKM)
library(countToFPKM)


read_counts_matrix <- function(counts_matrix_file, config_file, info_file, model_name) {

    print(counts_matrix_file) # print counts matrix file name
    print(model_name) # print tissue name
    gene_info <- read.csv(info_file) # read gene info file made in rnaseq_preprocess.py
    gene_info$size <- (gene_info$end_position-gene_info$start_position) # calculate gene size
    conf <- read.xlsx(config_file, model_name, colNames=TRUE) # read configuration sheet
    cmat_whole <- read.csv(counts_matrix_file, header=TRUE) # read counts matrix
    cmat_whole <- cmat_whole[(gene_info$entrezgene_id!="-"),] # remove unnamed genes
    gene_info <- gene_info[(gene_info$entrezgene_id!="-"),]  # remove unnamed genes
    genes <- gene_info$entrezgene_id  # get gene names

    # remove version numbers from ensembl id
    for ( j in 1:length(genes) ) {
        r <- genes[j]
        if ( grepl("\\\\.", r) ) {
            gen <- unlist(str_split(r, "\\\\."))[1]
            genes[j] <- gen
        }
    }

    # delete?
    samp_mat <- genes
    SampMetrics <- list()
    fragment_lengths <- c()
    samp_names <- c()
    groups <- unique(conf$Group)
    layouts <- c()

    # initialize groups
    for (group in groups) {
        SampMetrics[[group]][["CountMatrix"]] <- genes
        SampMetrics[[group]][["FragmentLengths"]] <- c()
        SampMetrics[[group]][["SampNames"]] <- c()
        SampMetrics[[group]][["Layout"]] <- c()
    }
    # add to group count matrices and insert lists
    for ( i in 1:length(conf$SampleName) ) {
        entry <- conf$SampleName[i]
        #group <- unlist(str_split(entry, "_", n=2))[2]
        group <- conf$Group[i]
        insert_size <- conf$FragmentLength[i]
        if ( entry %in% colnames(cmat_whole) ) {
            # get init values
            samp_mat <- SampMetrics[[group]][["CountMatrix"]]
            fragment_lengths <- SampMetrics[[group]][["FragmentLengths"]]
            samp_names <- SampMetrics[[group]][["SampNames"]]
            layouts <- SampMetrics[[group]][["Layout"]]

            # add replicate to values
            samp_mat <- cbind(samp_mat, cmat_whole[,entry])
            fragment_lengths <- c(fragment_lengths, conf$FragmentLength[i])
            layouts <- c(layouts, conf$Layout[i])
            samp_names <- c(samp_names, entry)

            # update values
            samp_mat <- SampMetrics[[group]][["CountMatrix"]] <- samp_mat
            SampMetrics[[group]][["FragmentLengths"]] <- fragment_lengths
            SampMetrics[[group]][["SampNames"]] <-  samp_names
            SampMetrics[[group]][["Layout"]] <- layouts

        } else {
            print(paste(c(entry, " not found in count matrix."),collapse="")) # inform that a sample is missing
        }
    }
    for ( group in groups ) { # for each study/batch group
        samp_mat <- SampMetrics[[group]][["CountMatrix"]]
        samp_mat <- samp_mat[,-1]
        samp_mat <- sapply(data.frame(samp_mat), as.numeric)
        samp_mat <- as.matrix(samp_mat) # convert df to matrix
        colnames(samp_mat) <- SampMetrics[[group]][["SampNames"]] # set column names to sample names
        SampMetrics[[group]][["CountMatrix"]] <- samp_mat # update counts matrix
        SampMetrics[[group]][["NumSamples"]] <- ncol(samp_mat) # set number of samples
        #SampMetrics[[group]][["FragmentLengths"]] <- fragment_lengths
        SampMetrics[[group]][["Entrez"]] <- as.character(gene_info$entrezgene_id) # store entrez ids
        SampMetrics[[group]][["GeneSizes"]] <- gene_info$size # store gene size
    }

    return(SampMetrics)
}


calculate_tpm <- function(SampMetrics) {

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


calculate_fpkm <- function(SampMetrics) {

    for ( i in 1:length(SampMetrics) ) {

        layout = SampMetrics[[i]][["Layout"]] # get layout

        # normalize wiht fpkm if paired, rpkm if single end, kill if something else is specified
        print(layout)
        stopifnot( layout == "paired-end" || layout == "single-end" )
        if ( layout == "paired-end" ) { # fpkm
            count_matrix <- SampMetrics[[i]][["CountMatrix"]]
            gene_size <- SampMetrics[[i]][["GeneSizes"]]
            mean_fragment_lengths <- SampMetrics[[i]][["FragmentLengths"]]
            fpkm_matrix <- do.call(cbind, lapply(1:ncol(count_matrix), function(j) {
                fpkm(count_matrix[,j], gene_size[j], mean_fragment_lengths[j])
            }))
            colnames(fpkm_matrix) <- colnames(count_matrix)
            SampMetrics[[i]][["FPKM_Matrix"]] <- fpkm_matrix

        } else { # rpkm
            count_matrix <- SampMetrics[[i]][["CountMatrix"]]
            gene_size <- SampMetrics[[i]][["GeneSizes"]]
            rpkm_matrix <- do.call(cbind, lapply(1:ncol(count_matrix), function(j) {
                rate = log(count_matrix[,j]) - log(gene_size[j])
                exp(rate - log(sum(count_matrix[,j])) + log(10e9))
            }))
            colnames(rpkm_matrix) <- colnames(count_matrix)
            SampMetrics[[i]][["FPKM_Matrix"]] <- rpkm_matrix

        }
    }

    return(SampMetrics)
}


calculate_z_score <- function(SampMetrics, norm_tech) {

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


cpm_filter <- function(SampMetrics, filt_options) {

    N_exp <- filt_options$replicate_ratio
    N_top <- filt_options$replicate_ratio_high
    min.count <- filt_options$min_count
    for ( i in 1:length(SampMetrics) ) {
        counts <- SampMetrics[[i]][["CountMatrix"]]
        ent <- SampMetrics[[i]][["Entrez"]]
        size <- SampMetrics[[i]][["GeneSizes"]]
        lib.size <- colSums(counts)
        MedianLibSize <- median(lib.size)
        if ( min.count=="default" ) {
            CPM.Cutoff <- 10000000/(median(colSums(counts))) # cite paper!!
        } else {
            CPM.Cutoff <- min.count/MedianLibSize*1e6
        }
        CPM <- cpm(counts,lib.size=lib.size)
        min.samples <- round(N_exp * ncol(counts))
        top.samples <- round(N_top * ncol(counts))
    
        f1 <- genefilter::kOverA(min.samples, CPM.Cutoff)
        flist <- genefilter::filterfun(f1)
        keep <- genefilter::genefilter(CPM, flist)
        SampMetrics[[i]][["Entrez"]] <- ent[keep]
        SampMetrics[[i]][["GeneSizes"]] <- size[keep]
        SampMetrics[[i]][["CountMatrix"]] <- counts[keep,]
        SampMetrics[[i]][["CPM_Matrix"]] <- CPM[keep,]
    
        # top percentile genes
        f1_top <- genefilter::kOverA(top.samples, CPM.Cutoff)
        flist_top <- genefilter::filterfun(f1_top)
        keep_top <- genefilter::genefilter(CPM, flist_top)
        SampMetrics[[i]][["Entrez_top"]] <- ent[keep_top]
    }

    SampMetrics <- calculate_z_score(SampMetrics, "CPM")

    return(SampMetrics)
}


TPM_quant_filter <- function(SampMetrics, filt_options, model_name) {

    N_exp <- filt_options$replicate_ratio
    N_top <- filt_options$replicate_ratio_high
    quant <- filt_options$quantile
    SampMetrics <- calculate_tpm(SampMetrics)

    for ( i in 1:length(SampMetrics) ) {
        counts <- SampMetrics[[i]][["CountMatrix"]]
        ent <- SampMetrics[[i]][["Entrez"]]
        size <- SampMetrics[[i]][["GeneSizes"]]
        tpm <- SampMetrics[[i]][["TPM_Matrix"]]
        tpm_fname <- paste(c("/home/jupyteruser/work/data/results/", model_name, "/TPM_Matrix_S", i, ".csv"),collapse="")
        write_tpm <- cbind(ent, tpm)
        write.csv(write_tpm, tpm_fname)

        min.samples <- round(N_exp * ncol(tpm))
        top.samples <- round(N_top * ncol(tpm))
        test_bools <- data.frame(gene=ent)
        test_bools_top <- test_bools

        for ( j in 1:ncol(tpm) ) {
            tpm_q <- tpm[,j]
            tpm_q <- tpm_q[tpm_q>0]
            q_cutoff <- quantile(tpm_q, prob=1-quant/100)
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

    SampMetrics <- calculate_z_score(SampMetrics, "TPM")

    return(SampMetrics)
}


zfpkm_filter <- function(SampMetrics, filt_options, model_name) {

    N_exp <- filt_options$replicate_ratio # ratio replicates for active
    N_top <- filt_options$replicate_ratio_high # ratio of replicates for high-confidence
    cutoff <- -3

    calculate_fpkm(SampMetrics)

    for ( i in 1:length(SampMetrics) ) {
        fmat <- SampMetrics[[i]][["FPKM_Matrix"]] # get fpkm matrix
        fdf <- data.frame(fmat) # convert to df
        ent <- SampMetrics[[i]][["Entrez"]] # get entrez ids
        size <- SampMetrics[[i]][["GenesSizes"]] # get gene sizes
        zmat <- zFPKM(fdf, assayName="FPKM") # calculate zFPKM
        zfpkm_fname <- paste(c("/home/jupyteruser/work/data/results/", model_name, "/zFPKM_Matrix_S", i, ".csv"),collapse="")
        write_zfpkm <- cbind(ent, zfpkm)
        write.csv(write_zfpkm, zfpkm_fname)

        #inames <- rownames(zmat)
    
        min.samples <- round(N_exp * ncol(zmat)) # min number of samples for active
        top.samples <- round(N_top * ncol(zmat)) # top number of samples for high-confidence

        f1 <- genefilter::kOverA(min.samples, cutoff)
        flist <- genefilter::filterfun(f1)
        keep <- genefilter::genefilter(zmat, flist)
        SampMetrics[[i]][["zFPKM_Matrix"]] <- zmat[keep,]
    
        # top percentile genes
        f1_top <- genefilter::kOverA(top.samples, cutoff)
        flist_top <- genefilter::filterfun(f1_top)
        keep_top <- genefilter::genefilter(zmat, flist_top)

        SampMetrics[[i]][["Entrez_top"]] <- ent[keep_top]
    }

    return(SampMetrics)
}


filter_counts <- function(SampMetrics, technique, filt_options, mart, model_name) {

    switch(technique,
         cpm = cpm_filter(SampMetrics, filt_options),
         zfpkm = zfpkm_filter(SampMetrics, filt_options, model_name),
         quantile = TPM_quant_filter(SampMetrics, filt_options, model_name))
}


save_rnaseq_tests <- function(counts_matrix_file, config_file, out_file, info_file, model_name,
                            replicate_ratio=0.5, sample_ratio=0.5, replicate_ratio_high=0.9, sample_ratio_high=0.9, 
                            technique="quantile", quantile=0.9, min_count=10) {
  
      # condense filter options
      filt_options <- list()
      if ( exists("replicate_ratio") ) {
          filt_options$replicate_ratio <- replicate_ratio
      } else {
          filt_options$replicate_ratio <- 0.2
      }
      if ( exists("sample_ratio") ) {
          filt_options$sample_ratio <- sample_ratio
      } else {
          filt_options$sample_ratio <- 0.5
      }
      if ( exists("quantile") ) {
          filt_options$quantile <- quantile
      } else {
          filt_options$quantile <- 50
      }
      if ( exists("min_count") ) {
          filt_options$min_count <- min_count
      } else {
          filt_options$min_count <- 1
      }
      if ( exists("replicate_ratio_high") ) {
          filt_options$replicate_ratio_high<- replicate_ratio_high
      } else {
          filt_options$replicate_ratio_high <- 0.9
      }
      if ( exists("sample_ratio_high") ) {
          filt_options$sample_ratio_high <- sample_ratio_high
      } else {
          filt_options$sample_ratio_high <- 0.9
      }

      print("Reading Counts Matrix")
      SampMetrics <- read_counts_matrix(counts_matrix_file, config_file, info_file, model_name) # read count matrix
      entrez_all <- SampMetrics[[1]][["Entrez"]] #get entrez ids
      print("Filtering Counts")
      SampMetrics <- filter_counts(SampMetrics, technique, filt_options, mart, model_name) # normalize and filter count
      expressedGenes <- c()
      topGenes <- c()

      for ( i in 1:length(SampMetrics) ) { # get high confidence and expressed genes for each study/batch number
          expressedGenes <- c(expressedGenes, SampMetrics[[i]][["Entrez"]])
          topGenes <- c(topGenes, SampMetrics[[i]][["Entrez_top"]])
      }

      expMat <- as.data.frame(table(expressedGenes)) # convert expression to df
      topMat <- as.data.frame(table(topGenes)) # convert high confidence to df
      nc <- length(SampMetrics) # number of columns
      expMat <- cbind(expMat, "Prop"=expMat$Freq/nc) # calculate proportion of studies/batches expressed
      topMat <- cbind(topMat, "Prop"=topMat$Freq/nc) # calculate proportion of studies/batch high-confidence
      SampMetrics[["ExpressionMatrix"]] <- expMat # store expression matrix for saving
      SampMetrics[["TopMatrix"]] <- topMat # store high confidence matrix for saving

      # get genes which are expressed and high-confidence according to use defined ratio
      SampMetrics[["ExpressedGenes"]] <- as.character(expMat$expressedGenes[expMat$Prop>=sample_ratio])
      SampMetrics[["TopGenes"]] <- as.character(topMat$topGenes[topMat$Prop>=sample_ratio_high])

      # create a table to write gene expression and high confidence to
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