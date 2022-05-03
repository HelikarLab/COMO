# Check if rlogs directory exists, From: https://stackoverflow.com/a/46008094
library("stringr")
username <- Sys.info()["user"]
work_dir <- str_interp("/home/${username}/work")

if (!dir.exists(str_interp("${work_dir}/py/rlogs"))) {
    dir.create(str_interp("${work_dir}/py/rlogs"))
}

# prevent messy messages from repeatedly writing to juypter
zz <- file(file.path("/home", username, "work", "py", "rlogs", "rnaseq.Rout"), open="wt")
sink(zz, type="message")

library(tidyverse)
library(limma)
library(edgeR)
library(genefilter)
library(biomaRt)
library(sjmisc)
library(zFPKM)
library(zFPKM)
library(stringr)
library(readxl)


read_counts_matrix <- function(counts_matrix_file, config_file, info_file, model_name) {

    #print(counts_matrix_file) # print counts matrix file name
    print(model_name) # print tissue name

    conf <- read_excel(config_file, sheet=model_name)# read configuration sheet
    cmat_whole <- read.csv(counts_matrix_file, header=TRUE) %>% arrange(., genes) # read counts matrix

    #print(head(cmat_whole))

    gene_info <- read.csv(info_file) %>%
        mutate(size=end_position-start_position) %>%
        arrange(., ensembl_gene_id) %>%
        filter(.$ensembl_gene_id %in% cmat_whole$genes)

    #print(head(gene_info))
    #print(tail(gene_info))
    #print(head(cmat_whole))
    #print(tail(cmat_whole))

    cmat_whole <- cmat_whole[(gene_info$entrezgene_id!="-"),] # remove unnamed genes
    gene_info <- gene_info[(gene_info$entrezgene_id!="-"),]  # remove unnamed genes
    genes <- gene_info$entrezgene_id  # get gene names
    #print(genes)

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
        SampMetrics[[group]][["SampleNames"]] <- c()
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
            samp_names <- SampMetrics[[group]][["SampleNames"]]
            layouts <- SampMetrics[[group]][["Layout"]]

            # add replicate to values
            samp_mat <- cbind(samp_mat, cmat_whole[,entry])
            fragment_lengths <- c(fragment_lengths, conf$FragmentLength[i])
            layouts <- c(layouts, conf$Layout[i])
            samp_names <- c(samp_names, entry)

            # update values
            samp_mat <- SampMetrics[[group]][["CountMatrix"]] <- samp_mat
            SampMetrics[[group]][["FragmentLengths"]] <- fragment_lengths
            SampMetrics[[group]][["SampleNames"]] <-  samp_names
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
        colnames(samp_mat) <- SampMetrics[[group]][["SampleNames"]] # set column names to sample names
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
        stopifnot( layout[1] == "paired-end" || layout == "single-end" )
        if ( layout == "paired-end" ) { # fpkm
            count_matrix <- SampMetrics[[i]][["CountMatrix"]]
            gene_size <- SampMetrics[[i]][["GeneSizes"]]
            mean_fragment_lengths <- SampMetrics[[i]][["FragmentLengths"]]
            fpkm_matrix <- do.call(cbind, lapply(1:ncol(count_matrix), function(j) {
                eff_len = gene_size - mean_fragment_lengths[j] + 1 # plus one to prevent div by 0
                N = sum(count_matrix[,j])
                exp( log(count_matrix[,j]) + log(1e9) - log(eff_len) - log(N) )
            })) + 1e-9
            colnames(fpkm_matrix) <- colnames(count_matrix)
            SampMetrics[[i]][["FPKM_Matrix"]] <- fpkm_matrix

        } else { # rpkm
            count_matrix <- SampMetrics[[i]][["CountMatrix"]]
            gene_size <- SampMetrics[[i]][["GeneSizes"]]
            rpkm_matrix <- do.call(cbind, lapply(1:ncol(count_matrix), function(j) {
                rate = log(count_matrix[,j]) - log(gene_size[j])
                exp(rate - log(sum(count_matrix[,j])) + log(1e9))
            })) + 1e-9
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


cpm_filter <- function(SampMetrics, filt_options, model_name) {

    N_exp <- filt_options$replicate_ratio
    N_top <- filt_options$replicate_ratio_high
    min.count <- filt_options$min_count
    for ( i in 1:length(SampMetrics) ) {
        study_number <- str_extract_all(SampMetrics[[i]][["SampleNames"]][1][1], "S\\d+")
        counts <- SampMetrics[[i]][["CountMatrix"]]
        ent <- SampMetrics[[i]][["Entrez"]]
        size <- SampMetrics[[i]][["GeneSizes"]]
        lib.size <- colSums(counts)
        CPM <- cpm(counts,lib.size=lib.size)
        cpm_fname <- file.path("home", "jupyteruser", "work", "data", "results",
                                model_name, paste0("CPM_Matrix_", study_number, ".csv"))
        write_cpm <- cbind(ent, CPM)
        write.csv(write_cpm, cpm_fname, row.names=FALSE)

        min.samples <- round(N_exp * ncol(counts))
        top.samples <- round(N_top * ncol(counts))
        test_bools <- data.frame(gene=ent)

        for ( j in 1:ncol(CPM) ) {
            cutoff <- ifelse( min.count=="default",
                     10e6/(median(sum(counts[,j]))),
                     1e6*min.count/(median(sum(counts[,j]))) )
            test_bools <- cbind(test_bools, as.integer(CPM[,j]>cutoff))
        }

        test_bools["gene"] <- NULL
        f1 <- genefilter::kOverA(min.samples, 0.9)
        flist <- genefilter::filterfun(f1)
        keep <- genefilter::genefilter(test_bools, flist)
        SampMetrics[[i]][["Entrez"]] <- ent[keep]
        SampMetrics[[i]][["GeneSizes"]] <- size[keep]
        SampMetrics[[i]][["CountMatrix"]] <- counts[keep,]
        SampMetrics[[i]][["CPM_Matrix"]] <- CPM[keep,]

        f1_top <- genefilter::kOverA(top.samples, 0.9)
        flist_top <- genefilter::filterfun(f1_top)
        keep_top <- genefilter::genefilter(test_bools, flist_top)

        SampMetrics[[i]][["Entrez_hc"]] <- ent[keep_top]
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
        study_number <- str_extract_all(SampMetrics[[i]][["SampleNames"]][1][1], "S\\d+")
        counts <- SampMetrics[[i]][["CountMatrix"]]
        ent <- SampMetrics[[i]][["Entrez"]]
        size <- SampMetrics[[i]][["GeneSizes"]]
        tpm <- SampMetrics[[i]][["TPM_Matrix"]]
        tpm_fname <- file.path("home", "jupyteruser", "work", "data", "results",
                        model_name, paste0("TPM_Matrix_", study_number, ".csv"))
        write_tpm <- cbind(ent, tpm)
        write.csv(write_tpm, tpm_fname, row.names=FALSE)

        min.samples <- round(N_exp * ncol(tpm))
        top.samples <- round(N_top * ncol(tpm))
        test_bools <- data.frame(gene=ent)

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

        SampMetrics[[i]][["Entrez_hc"]] <- ent[keep_top]
    }

    SampMetrics <- calculate_z_score(SampMetrics, "TPM")

    return(SampMetrics)
}



zfpkm_filter <- function(SampMetrics, filt_options, model_name) {

    N_exp <- filt_options$replicate_ratio # ratio replicates for active
    N_top <- filt_options$replicate_ratio_high # ratio of replicates for high-confidence
    cutoff <- -3

    SampMetrics <- calculate_fpkm(SampMetrics)

    for ( i in 1:length(SampMetrics) ) {
        study_number <- str_extract_all(SampMetrics[[i]][["SampleNames"]][1][1], "S\\d+")
        ent <- SampMetrics[[i]][["Entrez"]] # get entrez ids
        fmat <- SampMetrics[[i]][["FPKM_Matrix"]] # get fpkm matrix
        fdf <- data.frame(fmat) # convert to df
        fpkm_fname <- file.path("/home", username, "work", "data", "results",
                                 model_name, paste0("FPKM_Matrix_", study_number, ".csv"))
        write_fpkm <- cbind(ent, fdf)
        write.csv(write_fpkm, fpkm_fname, row.names=FALSE)
        missing_vals <- is.na(fdf) # get NA values from fdf
        fdf[missing_vals] <- 0 # set NA values to zero to prevent error in zfpkm calculation
        zmat <- zFPKM(fdf, assayName="FPKM") # calculate zFPKM
        #SampMetrics[[i]][["zFPKM_Matrix"]] <- zmat
        zmat[missing_vals] <- NA # set NA values back to NA
        zfpkm_fname <- file.path("/home", username, "work", "data", "results",
                                 model_name, paste0("zFPKM_Matrix_", study_number, ".csv"))
        write_zfpkm <- cbind(ent, zmat)
        write.csv(write_zfpkm, zfpkm_fname, row.names=FALSE)

        zfpkm_plot_dir <- file.path("/home", username, "work", "data", "results",
                                    model_name, "figures")

        if ( !file.exists(zfpkm_plot_dir) ) {
            dir.create(zfpkm_plot_dir)
        }
        study_number <- str_extract_all(SampMetrics[[i]][["SampleNames"]][1][1], "S\\d+")
        zfpkm_plotname <- file.path(zfpkm_plot_dir, paste0("zFPKM_plot_", study_number, ".pdf"))
        pdf(zfpkm_plotname)
        zFPKMPlot(fdf, assayName="FPKM")
        dev.off()


        #inames <- rownames(zmat)

        min.samples <- round(N_exp * ncol(zmat)) # min number of samples for active
        top.samples <- round(N_top * ncol(zmat)) # top number of samples for high-confidence

        f1 <- genefilter::kOverA(min.samples, cutoff)
        flist <- genefilter::filterfun(f1)
        keep <- genefilter::genefilter(zmat, flist)
        SampMetrics[[i]][["Entrez"]] <- ent[keep]

        # top percentile genes
        f1_top <- genefilter::kOverA(top.samples, cutoff)
        flist_top <- genefilter::filterfun(f1_top)
        keep_top <- genefilter::genefilter(zmat, flist_top)

        SampMetrics[[i]][["Entrez_hc"]] <- ent[keep_top]

    }

    return(SampMetrics)
}


filter_counts <- function(SampMetrics, technique, filt_options, model_name) {

    switch(technique,
         cpm = cpm_filter(SampMetrics, filt_options, model_name),
         zfpkm = zfpkm_filter(SampMetrics, filt_options, model_name),
         quantile = TPM_quant_filter(SampMetrics, filt_options, model_name))
}


save_rnaseq_tests <- function(counts_matrix_file, config_file, out_file, info_file, model_name,
                            replicate_ratio=0.5, batch_ratio=0.5, replicate_ratio_high=0.9, batch_ratio_high=0.9,
                            technique="quantile", quantile=0.9, min_count=10) {

      # condense filter options
      filt_options <- list()
      if ( exists("replicate_ratio") ) {
          filt_options$replicate_ratio <- replicate_ratio
      } else {
          filt_options$replicate_ratio <- 0.2
      }
      if ( exists("batch_ratio") ) {
          filt_options$batch_ratio <- batch_ratio
      } else {
          filt_options$batch_ratio <- 0.5
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
      if ( exists("batch_ratio_high") ) {
          filt_options$batch_ratio_high <- batch_ratio_high
      } else {
          filt_options$batch_ratio_high <- 0.9
      }

      print("Reading Counts Matrix")
      SampMetrics <- read_counts_matrix(counts_matrix_file, config_file, info_file, model_name) # read count matrix
      entrez_all <- SampMetrics[[1]][["Entrez"]] #get entrez ids
      print("Filtering Counts")
      SampMetrics <- filter_counts(SampMetrics, technique, filt_options, model_name) # normalize and filter count
      expressedGenes <- c()
      topGenes <- c()

      for ( i in 1:length(SampMetrics) ) { # get high confidence and expressed genes for each study/batch number
          expressedGenes <- c(expressedGenes, SampMetrics[[i]][["Entrez"]])
          topGenes <- c(topGenes, SampMetrics[[i]][["Entrez_hc"]])
      }

      expMat <- as.data.frame(table(expressedGenes)) # convert expression to df
      topMat <- as.data.frame(table(topGenes)) # convert high confidence to df
      nc <- length(SampMetrics) # number of columns
      expMat <- cbind(expMat, "Prop"=expMat$Freq/nc) # calculate proportion of studies/batches expressed
      topMat <- cbind(topMat, "Prop"=topMat$Freq/nc) # calculate proportion of studies/batch high-confidence
      SampMetrics[["ExpressionMatrix"]] <- expMat # store expression matrix for saving
      SampMetrics[["TopMatrix"]] <- topMat # store high confidence matrix for saving

      # get genes which are expressed and high-confidence according to use defined ratio
      SampMetrics[["ExpressedGenes"]] <- as.character(expMat$expressedGenes[expMat$Prop>=batch_ratio])
      SampMetrics[["TopGenes"]] <- as.character(topMat$topGenes[topMat$Prop>=batch_ratio_high])

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
