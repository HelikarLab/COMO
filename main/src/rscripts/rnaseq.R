# Check if rlogs directory exists, From: https://stackoverflow.com/a/46008094
library("stringr")
username <- Sys.info()["user"]
work_dir <- str_interp("/home/${username}/main")

if (!dir.exists(str_interp("${work_dir}/py/rlogs"))) {
    dir.create(str_interp("${work_dir}/py/rlogs"))
}

# prevent messy messages from repeatedly writing to juypter
zz <- file(file.path("/home", username, "main", "py", "rlogs", "rnaseq.Rout"), open="wt")
sink(zz, type="message")

library(tidyverse)
library(limma)
library(edgeR)
library(genefilter)
library(biomaRt)
library(sjmisc)
library(zFPKM)
library(stringr)
library(readxl)
library(dplyr)


read_counts_matrix <- function(counts_matrix_filepath, config_filepath, info_filepath, context_name) {
    config_object <- read_excel(config_filepath, sheet=context_name)# read configuration sheet
    counts_matrix <- read.csv(counts_matrix_filepath, header=TRUE) %>% arrange(., genes) # read counts matrix

    gene_info <- read.csv(info_filepath) %>%
        mutate(size=end_position-start_position) %>%
        arrange(., ensembl_gene_id) %>%
        filter(.$ensembl_gene_id %in% counts_matrix$genes)

    counts_matrix <- counts_matrix[(gene_info$entrezgene_id!="-"),] # remove unnamed genes
    gene_info <- gene_info[(gene_info$entrezgene_id!="-"),]  # remove unnamed genes
    genes <- gene_info$entrezgene_id  # get gene names

    # remove version numbers from ensembl id
    for ( i in 1:length(genes) ) {
        r <- genes[i]
        if ( grepl("\\\\.", r) ) {
            gen <- unlist(str_split(r, "\\\\."))[1]
            genes[i] <- gen
        }
    }
    
    SampMetrics <- list()
    groups <- unique(config_object$Group)
    
    
    testing_value <- "naivecd8_S11R1"

    # initialize groups
    for (group in groups) {
        SampMetrics[[group]][["CountMatrix"]] <- genes
        SampMetrics[[group]][["FragmentLengths"]] <- c()
        SampMetrics[[group]][["SampleNames"]] <- c()
        SampMetrics[[group]][["Layout"]] <- c()
    }

    # add to group count matrices and insert lists
    for ( i in 1:length(config_object$SampleName) ) {
        entry <- config_object$SampleName[i]  # CELL-TYPE_S##R##
        group <- config_object$Group[i]
        
        # These values will be added to the SampMetrix object
        matrix_value <- counts_matrix[, entry]
        fragment_length_value <- config_object$FragmentLength[i]
        layout_value <- config_object$Layout[i]
        sample_name_value <- entry
        
        
        # In case insert_size is NULL, set it to 0
        if (is.na(fragment_length_value)) {
            fragment_length_value <- 0
        }
        
        if ( entry %in% colnames(counts_matrix) ) {
            
            # get init values
            sample_matrix <- SampMetrics[[group]][["CountMatrix"]]
            fragment_lengths <- SampMetrics[[group]][["FragmentLengths"]]
            sample_names <- SampMetrics[[group]][["SampleNames"]]
            layouts <- SampMetrics[[group]][["Layout"]]
            
            # add replicate to values
            sample_matrix <- cbind(sample_matrix, matrix_value)
            fragment_lengths <- c(fragment_lengths, fragment_length_value)
            layouts <- c(layouts, layout_value)
            sample_names <- c(sample_names, sample_name_value)
            
            # update values
            SampMetrics[[group]][["CountMatrix"]] <- sample_matrix
            SampMetrics[[group]][["FragmentLengths"]] <- fragment_lengths
            SampMetrics[[group]][["SampleNames"]] <-  sample_names
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
        SampMetrics[[group]][["Entrez"]] <- as.character(gene_info$entrezgene_id) # store entrez ids
        SampMetrics[[group]][["GeneSizes"]] <- gene_info$size # store gene size
        SampMetrics[[group]][["StudyNumber"]] <- group
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

        layout <- SampMetrics[[i]][["Layout"]] # get layout

        # normalize with fpkm if paired, rpkm if single end, kill if something else is specified
        stopifnot( layout[1] == "paired-end" || layout == "single-end" )
        if ( layout == "paired-end" ) { # fpkm
            count_matrix <- SampMetrics[[i]][["CountMatrix"]]
            gene_size <- SampMetrics[[i]][["GeneSizes"]]
            mean_fragment_lengths <- SampMetrics[[i]][["FragmentLengths"]]
            fpkm_matrix <- do.call(cbind, lapply(1:ncol(count_matrix), function(j) {
                eff_len = gene_size - mean_fragment_lengths[j] + 1 # plus one to prevent div by 0
                N = sum(count_matrix[,j])
                exp( log(count_matrix[,j]) + log(1e9) - log(eff_len) - log(N) )
            }))
            fpkm_matrix[is.nan(fpkm_matrix)] <- 0
            colnames(fpkm_matrix) <- colnames(count_matrix)
            SampMetrics[[i]][["FPKM_Matrix"]] <- fpkm_matrix

        } else { # rpkm
            count_matrix <- SampMetrics[[i]][["CountMatrix"]]
            gene_size <- SampMetrics[[i]][["GeneSizes"]]
            rpkm_matrix <- do.call(cbind, lapply(1:ncol(count_matrix), function(j) {
                rate = log(count_matrix[,j]) - log(gene_size[j])
                exp(rate - log(sum(count_matrix[,j])) + log(1e9))
            }))
            rpkm_matrix[is.nan(rpkm_matrix)] <- 0
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


cpm_filter <- function(SampMetrics, filt_options, context_name, prep) {

    N_exp <- filt_options$replicate_ratio
    N_top <- filt_options$replicate_ratio_high
    min.count <- filt_options$min_count
    for ( i in 1:length(SampMetrics) ) {

        study_number <- SampMetrics[[i]][["StudyNumber"]]
        counts <- SampMetrics[[i]][["CountMatrix"]]
        ent <- SampMetrics[[i]][["Entrez"]]
        size <- SampMetrics[[i]][["GeneSizes"]]
        lib.size <- colSums(counts)
        CPM <- cpm(counts,lib.size=lib.size)
        cpm_fname <- file.path("home", "jupyteruser", "main", "data", "results",
                                context_name, prep, paste0("CPM_Matrix_", prep, "_", study_number, ".csv"))
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


TPM_quant_filter <- function(SampMetrics, filt_options, context_name, prep) {

    N_exp <- filt_options$replicate_ratio
    N_top <- filt_options$replicate_ratio_high
    quant <- filt_options$quantile
    SampMetrics <- calculate_tpm(SampMetrics)

    for ( i in 1:length(SampMetrics) ) {

        study_number <- SampMetrics[[i]][["StudyNumber"]]
        counts <- SampMetrics[[i]][["CountMatrix"]]
        ent <- SampMetrics[[i]][["Entrez"]]
        size <- SampMetrics[[i]][["GeneSizes"]]
        tpm <- SampMetrics[[i]][["TPM_Matrix"]]
        tpm_fname <- file.path("home", "jupyteruser", "main", "data", "results",
                        context_name, prep, paste0("TPM_Matrix_", prep, "_", study_number, ".csv"))
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



zfpkm_filter <- function(SampMetrics, filt_options, context_name, prep) {
    N_exp <- filt_options$replicate_ratio # ratio replicates for active
    N_top <- filt_options$replicate_ratio_high # ratio of replicates for high-confidence
    cutoff <- filt_options$min_zfpkm
    
    SampMetrics <- calculate_fpkm(SampMetrics)
    for ( i in seq_along(SampMetrics)) {
        study_number <- SampMetrics[[i]][["StudyNumber"]]
        
        entrez_ids <- SampMetrics[[i]][["Entrez"]] # get entrez ids
        fpkm_matrix <- SampMetrics[[i]][["FPKM_Matrix"]] # get fpkm matrix
        fpkm_df <- data.frame(fpkm_matrix) # convert to df
        fpkm_df[rowSums(fpkm_df[])>0,]
        fpkm_filename <- file.path("/home", username, "main", "data", "results", context_name, prep, paste0("FPKM_Matrix_", prep, "_", study_number, ".csv"))
        write_fpkm <- cbind(entrez_ids, fpkm_df)
        colnames(write_fpkm)[1] <- "ENTREZ_GENE_ID"
        write.csv(write_fpkm, fpkm_filename, row.names=FALSE)
        
        minimums <- fpkm_df == 0
        nas <- is.na(fpkm_df) == 1
        
        # calculate zFPKM
        zmat <- zFPKM(fpkm_df, min_thresh=0, assayName="FPKM")
        zmat[minimums] <- -4 # instead of -inf set to lower limit
        
        zfpkm_fname <- file.path("/home", username, "main", "data", "results", context_name, prep, paste0("zFPKM_Matrix_", prep, "_", study_number, ".csv"))
        write_zfpkm <- cbind(entrez_ids, zmat)
        colnames(write_zfpkm)[1] <- "ENTREZ_GENE_ID"
        write.csv(write_zfpkm, zfpkm_fname, row.names=FALSE)
        
        zfpkm_plot_dir <- file.path("/home", username, "main", "data", "results", context_name, prep, "figures")
        if ( !file.exists(zfpkm_plot_dir) ) {
            dir.create(zfpkm_plot_dir)
        }
        
        zfpkm_plotname <- file.path(zfpkm_plot_dir, paste0("zFPKM_plot_", study_number, ".pdf"))
        pdf(zfpkm_plotname)
        zFPKMPlot(fpkm_df, min_thresh=min(fpkm_df), assayName="FPKM")
        dev.off()
        
        min.samples <- round(N_exp * ncol(zmat)) # min number of samples for active
        top.samples <- round(N_top * ncol(zmat)) # top number of samples for high-confidence
        
		# active genes
        f1 <- genefilter::kOverA(min.samples, cutoff)
        flist <- genefilter::filterfun(f1)
        keep <- genefilter::genefilter(zmat, flist)
        SampMetrics[[i]][["Entrez"]] <- entrez_ids[keep]
        
        # top percentile genes
        f1_top <- genefilter::kOverA(top.samples, cutoff)
        flist_top <- genefilter::filterfun(f1_top)
        keep_top <- genefilter::genefilter(zmat, flist_top)
        SampMetrics[[i]][["Entrez_hc"]] <- entrez_ids[keep_top]
    }

    
    return(SampMetrics)
}


umi_filter <- function(SampMetrics, filt_options, context_name) {
    prep <- "scrna"
    N_exp <- filt_options$replicate_ratio # ratio replicates for active
    N_top <- filt_options$replicate_ratio_high # ratio of replicates for high-confidence
    cutoff <- filt_options$min_zfpkm

    #SampMetrics <- calculate_fpkm(SampMetrics)
    for ( i in 1:length(SampMetrics) ) {
        study_number <- SampMetrics[[i]][["StudyNumber"]]
        ent <- SampMetrics[[i]][["Entrez"]] # get entrez ids
        umat <- SampMetrics[[i]][["CountMatrix"]] # get fpkm matrix
        udf <- data.frame(umat) # convert to df
        udf[rowSums(udf[])>0,]
        minimums <- udf == 0
        nas <- is.na(udf) == 1
        zmat <- zFPKM(udf, min_thresh=0, assayName="UMI") # calculate zFPKM
        zmat[minimums] <- -4 # instead of -inf set to lower limit
        zumi_fname <- file.path("/home", username, "main", "data", "results", context_name, prep, paste0("zUMI_Matrix_", prep, "_", study_number, ".csv"))

        write_zumi <- cbind(ent, zmat)
        colnames(write_zumi)[1] <- "ENTREZ_GENE_ID"
        write.csv(write_zumi, zumi_fname, row.names=FALSE)

        zumi_plot_dir <- file.path("/home", username, "main", "data", "results",
                                    context_name, prep, "figures")

        if ( !file.exists(zumi_plot_dir) ) {
            dir.create(zumi_plot_dir)
        }

        batch_size <- 12
        plot_batches <- ceiling(ncol(udf)/batch_size)

        if ( plot_batches < 2 ) {
            zumi_plotname <- file.path(zumi_plot_dir, paste0("zumi_plot_", study_number, ".pdf"))
            pdf(zumi_plotname)
            zFPKMPlot(udf, min_thresh=min(udf), assayName="UMI")
            dev.off()

        } else {

            for ( j in 1:(plot_batches-1) ) {
                zumi_plotname <- file.path(zumi_plot_dir, paste0("zumi_plot_", study_number, "_", j, ".pdf"))
                pdf(zumi_plotname)
                samps <- c()
                jmin <- (batch_size*(j-1))+1
                jmax <- batch_size*j
                samps <- jmin:jmax

                while ( samps[length(samps)] > ncol(udf) ) {
                  samps <- samps[1:length(samps)-1]
                }

                zFPKMPlot(
                  udf[,samps],
                  min_thresh=0,
                  assayName="UMI"
                )
                dev.off()
            }
        }

        min.samples <- round(N_exp * ncol(zmat)) # min number of samples for active
        top.samples <- round(N_top * ncol(zmat)) # top number of samples for high-confidence

		# active genes
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


filter_counts <- function(SampMetrics, technique, filt_options, context_name, prep) {
    switch(
      technique,
      cpm = cpm_filter(SampMetrics, filt_options, context_name, prep),
      zfpkm = zfpkm_filter(SampMetrics, filt_options, context_name, prep),
      quantile = TPM_quant_filter(SampMetrics, filt_options, context_name, prep),
      umi = umi_filter(SampMetrics, filt_options, context_name)
    )
}

save_rnaseq_tests <- function(
  counts_matrix_file,
  config_file,
  out_file,
  info_file,
  context_name,
  prep="total",
  replicate_ratio=0.5,
  batch_ratio=0.5,
  replicate_ratio_high=0.9,
  batch_ratio_high=0.9,
  technique="quantile",
  quantile=0.9,
  min_count=10,
  min_zfpkm=-3
) {
    
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
    if ( exists("min_zfpkm") ) {
    filt_options$min_zfpkm <- min_zfpkm
    } else {
    filt_options$min_zfpkm <- -3
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

    if ( prep == "scrna" ) {
      technique <- "umi"
      print("Note: Single cell filtration does not normalize and assumes counts are counted with UMI")
    }
    
    SampMetrics <- read_counts_matrix(counts_matrix_file, config_file, info_file, context_name) # read count matrix
    entrez_all <- SampMetrics[[1]][["Entrez"]] #get entrez ids
    
    SampMetrics <- filter_counts(SampMetrics, technique, filt_options, context_name, prep) # normalize and filter count
    expressedGenes <- c()
    topGenes <- c()
    for ( i in 1:length(SampMetrics) ) { # get high confidence and expressed genes for each study/batch number
        expressedGenes <- c(expressedGenes, SampMetrics[[i]][["Entrez"]])
        topGenes <- c(topGenes, SampMetrics[[i]][["Entrez_hc"]])
    }
    SampMetrics <- read_counts_matrix(counts_matrix_file, config_file, info_file, context_name) # read count matrix

    entrez_all <- SampMetrics[[1]][["Entrez"]] #get entrez ids
    
    SampMetrics <- filter_counts(SampMetrics, technique, filt_options, context_name, prep) # normalize and filter count
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
    for ( i in 1:nrow(write_table) ) {
        if (as.character(write_table[i,1]) %in% as.character(SampMetrics$ExpressedGenes)) {
            write_table[i,2] <- 1
        }
        if (as.character(write_table$entrez_all[i]) %in% as.character(SampMetrics$TopGenes)) {
            write_table[i,3] <- 1
        }
    }
    header <- c("ENTREZ_GENE_ID", "expressed", "high")
    #write_table <- rbind(header, write_table)
    colnames(write_table) <- header
    write.csv(write_table, out_file, row.names=FALSE, col.names=FALSE)
}
