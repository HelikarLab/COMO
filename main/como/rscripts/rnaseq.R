suppressPackageStartupMessages(library("biomaRt"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("edgeR"))
suppressPackageStartupMessages(library("genefilter"))
suppressPackageStartupMessages(library("limma"))
suppressPackageStartupMessages(library("readr"))
suppressPackageStartupMessages(library("readxl"))
suppressPackageStartupMessages(library("sjmisc"))
suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("zFPKM"))

# Check if rlogs directory exists, From: https://stackoverflow.com/a/46008094
# Then prevent messy messages from repeatedly writing to juypter
work_dir <- getwd()
r_log_directory <- stringr::str_interp("${work_dir}/logs")
if (!dir.exists(r_log_directory)) { dir.create(r_log_directory) }
zz <- file(file.path(r_log_directory, "rnaseq.Rout"), open = "wt")
sink(zz, type = "message")

# Create a variable x set equal to 1


read_counts_matrix <- function(counts_matrix_filepath, config_filepath, info_filepath, context_name) {
    config_object <- readxl::read_excel(config_filepath, sheet = context_name)
    counts_matrix <- data.frame(readr::read_csv(file = counts_matrix_filepath)) %>% dplyr::arrange(., ensembl_gene_id)
    gene_info <- data.frame(readr::read_csv(info_filepath)) %>%
      dplyr::mutate(size = end_position - start_position) %>%  # Calculate gene size
      dplyr::arrange(.$ensembl_gene_id) %>%
      dplyr::filter(.$entrez_gene_id != "-") %>%  # Remove un-named genes
      dplyr::filter(.$ensembl_gene_id %in% counts_matrix$ensembl_gene_id)  # Only keep genes that are in counts_matrix

    # Select rows where "genes" in counts_matrix are also present in gene_info
    counts_matrix <- counts_matrix[(counts_matrix$ensembl_gene_id %in% gene_info$ensembl_gene_id),]
    entrez_genes <- gene_info$entrez_gene_id

    # remove version numbers from ensembl id
    for (i in seq_along(entrez_genes)) {
        row <- entrez_genes[i]
        if (grepl("\\\\.", row)) {
            gen <- unlist(stringr::str_split(row, "\\\\."))[1]
            entrez_genes[i] <- gen
        }
    }


    # initialize groups and pre-allocate space
    sample_metrics <- list()
    num_samples <- length(config_object$sample_name)
    groups <- unique(config_object$study)
    for (group in groups) {
        sample_metrics[[group]] <- list(
          CountMatrix = matrix(data = NA, nrow = nrow(counts_matrix), ncol = num_samples),
          FragmentLengths = numeric(length = num_samples),
          SampleNames = character(length = num_samples),
          Layout = character(length = num_samples)
        )
    }

    # add to group count matrices and insert lists
    for (i in seq_along(config_object$sample_name)) {
        entry <- config_object$sample_name[i]  # CELL-TYPE_S##R##
        group <- config_object$study[i]

        # Test if `entry` is not in counts_matrix columns
        if (!(entry %in% colnames(counts_matrix))) {
            print(paste(entry, "not found in count matrix."))
            next
        }

        # In case insert_size is NULL, set it to 0
        fragment_length_value <- ifelse(is.na(config_object$fragment_length[i]), 0, config_object$fragment_length[i])

        # update sample_metrics with current replicate values
        col_index <- which(colnames(counts_matrix) == entry)  # Use column index for faster access
        sample_metrics[[group]][["CountMatrix"]][, i] <- as.numeric(counts_matrix[, col_index])
        sample_metrics[[group]][["FragmentLengths"]][i] <- fragment_length_value
        sample_metrics[[group]][["SampleNames"]][i] <- entry
        sample_metrics[[group]][["Layout"]][i] <- config_object$layout[i]
    }

    for (group in groups) { # for each study/batch group
        samp_mat <- sample_metrics[[group]][["CountMatrix"]]
        samp_mat <- apply(samp_mat, FUN = as.numeric, MARGIN = 2, simplify = FALSE)  # '2' means apply over columns
        write.csv(samp_mat, file = "/home/joshl/projects/COMO/main/samp_mat.csv")
        samp_mat <- t(samp_mat)

        colnames(samp_mat) <- sample_metrics[[group]][["SampleNames"]] # set column names to sample names
        sample_metrics[[group]][["CountMatrix"]] <- samp_mat # update counts matrix
        sample_metrics[[group]][["NumSamples"]] <- ncol(samp_mat) # set number of samples
        sample_metrics[[group]][["Entrez"]] <- as.character(gene_info$entrez_gene_id) # store entrez ids
        sample_metrics[[group]][["GeneSizes"]] <- gene_info$size # store gene size
        sample_metrics[[group]][["StudyNumber"]] <- group
    }
    return(sample_metrics)
}


calculate_tpm <- function(sample_metrics) {

    for (i in seq_along(sample_metrics)) {
        count_matrix <- sample_metrics[[i]][["CountMatrix"]]
        gene_size <- sample_metrics[[i]][["GeneSizes"]]
        tpm_matrix <- do.call(cbind, lapply(seq_len(ncol(count_matrix)), function(j) {
            rate <- log(count_matrix[, j]) - log(gene_size[j])
            denom <- log(sum(exp(rate)))
            exp(rate - denom + log(1e6))
        }))
        colnames(tpm_matrix) <- colnames(count_matrix)
        sample_metrics[[i]][["TPM_Matrix"]] <- tpm_matrix
    }
    
    return(sample_metrics)
}


calculate_fpkm <- function(sample_metrics) {
    
    for (i in seq_along(sample_metrics)) {
        
        layout <- sample_metrics[[i]][["Layout"]] # get layout

        if (layout[1] == "paired-end") { # fpkm
            # Print components of sample_metrics[[i]], which is of type list
            count_matrix <- sample_metrics[[i]][["CountMatrix"]]
            gene_size <- sample_metrics[[i]][["GeneSizes"]]
            mean_fragment_lengths <- sample_metrics[[i]][["FragmentLengths"]]

            fpkm_matrix <- matrix(nrow = nrow(count_matrix), ncol = ncol(count_matrix))
            print(class(count_matrix))
            print(head(colSums(count_matrix)))
            for (j in seq_len(ncol(count_matrix))) {
                eff_len <- gene_size - mean_fragment_lengths[j] + 1 # plus one to prevent div by 0

                print(paste("j:", j))
                print(paste("ncol:", ncol(count_matrix)))
                N <- sum(count_matrix[, j])

                fpkm_matrix[, j] <- exp(log(count_matrix[, j]) + log(1e9) - log(eff_len) - log(N))
            }
            fpkm_matrix[is.nan(fpkm_matrix)] <- 0
            colnames(fpkm_matrix) <- colnames(count_matrix)
            sample_metrics[[i]][["FPKM_Matrix"]] <- fpkm_matrix

        } else if (layout[1] == "single-end" ) { # rpkm
            count_matrix <- sample_metrics[[i]][["CountMatrix"]]
            gene_size <- sample_metrics[[i]][["GeneSizes"]]
            rpkm_matrix <- do.call(cbind, lapply(seq_len(ncol(count_matrix)), function(j) {
                rate <- log(count_matrix[, j]) - log(gene_size[j])
                exp(rate - log(sum(count_matrix[, j])) + log(1e9))
            }))
            rpkm_matrix[is.nan(rpkm_matrix)] <- 0
            colnames(rpkm_matrix) <- colnames(count_matrix)
            sample_metrics[[i]][["FPKM_Matrix"]] <- rpkm_matrix
        } else {
            stop(paste0("Invalid layout for sample '", sample_metrics[[i]][["StudyNumber"]], "'. Must be paired-end or single-end"))
        }
    }
    return(sample_metrics)
}

calculate_z_score <- function(sample_metrics, norm_tech) {
    
    for (i in seq_along(sample_metrics)) {
        if (norm_tech == "CPM") {
            tmat <- sample_metrics[[i]][["CPM_Matrix"]]
        } else if (norm_tech == "TPM") {
            tmat <- sample_metrics[[i]][["TPM_Matrix"]]
        }
        zmat <- matrix(nrow = nrow(tmat), ncol = ncol(tmat))
        
        for (j in seq_len(ncol(tmat))) {
            tvec <- tmat[, j]
            logvec <- log2(tvec)
            logvec[is.infinite(logvec)] <- NA
            zvec <- scale(logvec, center = TRUE, scale = TRUE)
            zmat[, j] <- zvec
        }
        zmat <- data.frame(zmat)
        colnames(zmat) <- colnames(tmat)
        sample_metrics[[i]][["Zscore"]] <- zmat
    }
    
    return(sample_metrics)
}


cpm_filter <- function(sample_metrics, filt_options, context_name, prep) {
    
    N_exp <- filt_options$replicate_ratio
    N_top <- filt_options$replicate_ratio_high
    min.count <- filt_options$min_count
    for (i in seq_along(sample_metrics)) {
        
        study_number <- sample_metrics[[i]][["StudyNumber"]]
        counts <- sample_metrics[[i]][["CountMatrix"]]
        ent <- sample_metrics[[i]][["Entrez"]]
        size <- sample_metrics[[i]][["GeneSizes"]]
        lib.size <- colSums(counts)
        CPM <- edgeR::cpm(counts, lib.size = lib.size)
        cpm_fname <- file.path(work_dir, "data", "results", context_name, prep, paste0("CPM_Matrix_", prep, "_", study_number, ".csv"))
        write_cpm <- cbind(ent, CPM)
        write.csv(write_cpm, cpm_fname, row.names = FALSE)
        
        min.samples <- round(N_exp * ncol(counts))
        top.samples <- round(N_top * ncol(counts))
        test_bools <- data.frame(gene = ent)
        
        for (j in seq_len(ncol(CPM))) {
            cutoff <- ifelse(min.count == "default",
                             10e6 / (median(sum(counts[, j]))),
                             1e6 * min.count / (median(sum(counts[, j]))))
            test_bools <- cbind(test_bools, as.integer(CPM[, j] > cutoff))
        }
        
        test_bools["gene"] <- NULL
        f1 <- genefilter::kOverA(min.samples, 0.9)
        flist <- genefilter::filterfun(f1)
        keep <- genefilter::genefilter(test_bools, flist)
        sample_metrics[[i]][["Entrez"]] <- ent[keep]
        sample_metrics[[i]][["GeneSizes"]] <- size[keep]
        sample_metrics[[i]][["CountMatrix"]] <- counts[keep,]
        sample_metrics[[i]][["CPM_Matrix"]] <- CPM[keep,]
        
        f1_top <- genefilter::kOverA(top.samples, 0.9)
        flist_top <- genefilter::filterfun(f1_top)
        keep_top <- genefilter::genefilter(test_bools, flist_top)
        
        sample_metrics[[i]][["Entrez_hc"]] <- ent[keep_top]
    }
    
    sample_metrics <- calculate_z_score(sample_metrics, "CPM")
    
    return(sample_metrics)
}


TPM_quant_filter <- function(sample_metrics, filt_options, context_name, prep) {
    
    N_exp <- filt_options$replicate_ratio
    N_top <- filt_options$replicate_ratio_high
    quant <- filt_options$quantile
    sample_metrics <- calculate_tpm(sample_metrics)
    
    for (i in seq_along(sample_metrics)) {
        
        study_number <- sample_metrics[[i]][["StudyNumber"]]
        counts <- sample_metrics[[i]][["CountMatrix"]]
        ent <- sample_metrics[[i]][["Entrez"]]
        size <- sample_metrics[[i]][["GeneSizes"]]
        tpm <- sample_metrics[[i]][["TPM_Matrix"]]
        tpm_fname <- file.path(work_dir, "data", "results", context_name, prep, paste0("TPM_Matrix_", prep, "_", study_number, ".csv"))
        write_tpm <- cbind(ent, tpm)
        write.csv(write_tpm, tpm_fname, row.names = FALSE)
        
        min.samples <- round(N_exp * ncol(tpm))
        top.samples <- round(N_top * ncol(tpm))
        test_bools <- data.frame(gene = ent)
        
        for (j in seq_len(ncol(tpm))) {
            tpm_q <- tpm[, j]
            tpm_q <- tpm_q[tpm_q > 0]
            q_cutoff <- quantile(tpm_q, prob = 1 - quant / 100)
            #q_cutoff_top <- quantile(tpm_q, prob=1-perc_top/100)
            #bools <- data.frame(as.integer(tpm[,j]>q_cutoff))
            #bools_top <- data.frame(as.integer(tpm[,j]>q_cutoff_top))
            test_bools <- cbind(test_bools, as.integer(tpm[, j] > q_cutoff))
        }
        
        test_bools["gene"] <- NULL
        f1 <- genefilter::kOverA(min.samples, 0.9)
        flist <- genefilter::filterfun(f1)
        keep <- genefilter::genefilter(test_bools, flist)
        sample_metrics[[i]][["Entrez"]] <- ent[keep]
        sample_metrics[[i]][["GeneSizes"]] <- size[keep]
        sample_metrics[[i]][["CountMatrix"]] <- counts[keep,]
        sample_metrics[[i]][["TPM_Matrix"]] <- tpm[keep,]
        
        f1_top <- genefilter::kOverA(top.samples, 0.9)
        flist_top <- genefilter::filterfun(f1_top)
        keep_top <- genefilter::genefilter(test_bools, flist_top)
        
        sample_metrics[[i]][["Entrez_hc"]] <- ent[keep_top]
    }
    
    sample_metrics <- calculate_z_score(sample_metrics, "TPM")
    
    return(sample_metrics)
}


zfpkm_filter <- function(sample_metrics, filt_options, context_name, prep) {
    N_exp <- filt_options$replicate_ratio # ratio replicates for active
    N_top <- filt_options$replicate_ratio_high # ratio of replicates for high-confidence
    cutoff <- filt_options$min_zfpkm

    sample_metrics <- calculate_fpkm(sample_metrics)
    for (i in seq_along(sample_metrics)) {
        study_number <- sample_metrics[[i]][["StudyNumber"]]

        entrez_ids <- sample_metrics[[i]][["Entrez"]] # get entrez ids
        fpkm_matrix <- sample_metrics[[i]][["FPKM_Matrix"]] # get fpkm matrix
        fpkm_df <- data.frame(fpkm_matrix) # convert to df
        fpkm_df[rowSums(fpkm_df[]) > 0,]
        fpkm_filename <- file.path(work_dir, "data", "results", context_name, prep, paste0("FPKM_Matrix_", prep, "_", study_number, ".csv"))
        write_fpkm <- cbind(entrez_ids, fpkm_df)
        colnames(write_fpkm)[1] <- "ENTREZ_GENE_ID"
        write.csv(write_fpkm, fpkm_filename, row.names = FALSE)

        minimums <- fpkm_df == 0
        na_values <- is.na(fpkm_df) == 1

        # calculate zFPKM
        zmat <- zFPKM::zFPKM(fpkm_df, min_thresh = 0, assayName = "FPKM")
        zmat[minimums] <- -4 # instead of -inf set to lower limit
        
        zfpkm_fname <- file.path(work_dir, "data", "results", context_name, prep, paste0("zFPKM_Matrix_", prep, "_", study_number, ".csv"))
        write_zfpkm <- dplyr::bind_cols(entrez_ids, zmat)
        colnames(write_zfpkm)[1] <- "ENTREZ_GENE_ID"
        write.csv(write_zfpkm, zfpkm_fname, row.names = FALSE)
        
        zfpkm_plot_dir <- file.path(work_dir, "data", "results", context_name, prep, "figures")
        # zfpkm_plot_dir <- file.path("/home", username, "main", "data", "results", context_name, prep, "figures")
        if (!file.exists(zfpkm_plot_dir)) {
            dir.create(zfpkm_plot_dir)
        }
        
        zfpkm_plotname <- file.path(zfpkm_plot_dir, paste0("zFPKM_plot_", study_number, ".pdf"))
        pdf(zfpkm_plotname)
        zFPKM::zFPKMPlot(fpkm_df, min_thresh = min(fpkm_df), assayName = "FPKM")
        dev.off()
        
        min.samples <- round(N_exp * ncol(zmat)) # min number of samples for active
        top.samples <- round(N_top * ncol(zmat)) # top number of samples for high-confidence
        
        # active genes
        f1 <- genefilter::kOverA(min.samples, cutoff)
        flist <- genefilter::filterfun(f1)
        keep <- genefilter::genefilter(zmat, flist)
        sample_metrics[[i]][["Entrez"]] <- entrez_ids[keep]
        
        # top percentile genes
        f1_top <- genefilter::kOverA(top.samples, cutoff)
        flist_top <- genefilter::filterfun(f1_top)
        keep_top <- genefilter::genefilter(zmat, flist_top)
        sample_metrics[[i]][["Entrez_hc"]] <- entrez_ids[keep_top]
    }

    return(sample_metrics)
}


umi_filter <- function(sample_metrics, filt_options, context_name) {
    prep <- "scrna"
    N_exp <- filt_options$replicate_ratio # ratio replicates for active
    N_top <- filt_options$replicate_ratio_high # ratio of replicates for high-confidence
    cutoff <- filt_options$min_zfpkm


    #sample_metrics <- calculate_fpkm(sample_metrics)
    for (i in seq_along(sample_metrics)) {
        study_number <- sample_metrics[[i]][["StudyNumber"]]

        entrez_ids <- sample_metrics[[i]][["Entrez"]] # get entrez ids
        count_matrix <- data.frame(sample_metrics[[i]][["CountMatrix"]])
        # Convert count matrix to numeric
        count_matrix <- sapply(count_matrix, as.numeric)
        count_matrix[rowSums(count_matrix[]) > 0,]
        minimums <- count_matrix == 0

        zmat <- zFPKM::zFPKM(count_matrix, min_thresh = 0, assayName = "UMI")
        zmat[minimums] <- -4 # instead of -inf set to lower limit

        write_zumi <- dplyr::bind_cols(entrez_ids, zmat)
        # write_zumi <- cbind(ent, zmat)
        colnames(write_zumi)[1] <- "ENTREZ_GENE_ID"
        zumi_fname <- file.path(work_dir, "data", "results", context_name, prep, paste0("zUMI_Matrix_", prep, "_", study_number, ".csv"))
        write.csv(write_zumi, zumi_fname, row.names = FALSE)

        zumi_plot_dir <- file.path(work_dir, "data", "results", context_name, prep, "figures")
        # zumi_plot_dir <- file.path("/home", username, "main", "data", "results", context_name, prep, "figures")

        if (!file.exists(zumi_plot_dir)) {
            dir.create(zumi_plot_dir)
        }

        batch_size <- 12
        plot_batches <- ceiling(ncol(count_matrix) / batch_size)

        if (plot_batches < 2) {
            zumi_plotname <- file.path(zumi_plot_dir, paste0("zumi_plot_", study_number, ".pdf"))
            pdf(zumi_plotname)
            zFPKM::zFPKMPlot(count_matrix, min_thresh = min(count_matrix), assayName = "UMI")
            dev.off()

        } else {

            for (j in 1:(plot_batches - 1)) {
                zumi_plotname <- file.path(zumi_plot_dir, paste0("zumi_plot_", study_number, "_", j, ".pdf"))
                pdf(zumi_plotname)
                samps <- c()
                jmin <- (batch_size * (j - 1)) + 1
                jmax <- batch_size * j
                samps <- jmin:jmax

                while (samps[length(samps)] > ncol(count_matrix)) {
                    samps <- samps[seq_along(samps) - 1]
                }

                zFPKM::zFPKMPlot(
                  count_matrix[, samps],
                  min_thresh = 0,
                  assayName = "UMI"
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
        sample_metrics[[i]][["Entrez"]] <- entrez_ids[keep]

        # top percentile genes
        f1_top <- genefilter::kOverA(top.samples, cutoff)
        flist_top <- genefilter::filterfun(f1_top)
        keep_top <- genefilter::genefilter(zmat, flist_top)

        sample_metrics[[i]][["Entrez_hc"]] <- entrez_ids[keep_top]

    }

    return(sample_metrics)
}


filter_counts <- function(sample_metrics, technique, filt_options, context_name, prep) {
    switch(
      technique,
      cpm = cpm_filter(sample_metrics, filt_options, context_name, prep),
      zfpkm = zfpkm_filter(sample_metrics, filt_options, context_name, prep),
      quantile = TPM_quant_filter(sample_metrics, filt_options, context_name, prep),
      umi = umi_filter(sample_metrics, filt_options, context_name)
    )
}

save_rnaseq_tests <- function(
    counts_matrix_file,
    config_filepath,
    out_file,
    info_file,
    context_name,
    prep = "total",
    replicate_ratio = 0.5,
    batch_ratio = 0.5,
    replicate_ratio_high = 0.9,
    batch_ratio_high = 0.9,
    technique = "quantile",
    quantile = 0.9,
    min_count = 10,
    min_zfpkm = -3
) {
    
    # condense filter options
    filt_options <- list()
    filt_options$replicate_ratio <- replicate_ratio
    filt_options$batch_ratio <- batch_ratio
    filt_options$quantile <- quantile
    filt_options$min_count <- min_count
    filt_options$min_zfpkm <- min_zfpkm
    filt_options$replicate_ratio_high <- replicate_ratio_high
    filt_options$batch_ratio_high <- batch_ratio_high
    
    if (prep == "scrna") {
        technique <- "umi"
        print("Note: Single cell filtration does not normalize and assumes counts are counted with UMI")
    }

    print("1")
    sample_metrics <- read_counts_matrix(counts_matrix_file, config_filepath, info_file, context_name) # read count matrix
    print("2")
    entrez_all <- sample_metrics[[1]][["Entrez"]] #get entrez ids
    
    sample_metrics <- filter_counts(sample_metrics, technique, filt_options, context_name, prep) # normalize and filter count
    print("3")
    expressedGenes <- c()
    topGenes <- c()
    for (i in seq_along(sample_metrics)) { # get high confidence and expressed genes for each study/batch number
        expressedGenes <- c(expressedGenes, sample_metrics[[i]][["Entrez"]])
        topGenes <- c(topGenes, sample_metrics[[i]][["Entrez_hc"]])
    }
    print("4")

    expMat <- as.data.frame(table(expressedGenes)) # convert expression to df
    topMat <- as.data.frame(table(topGenes)) # convert high confidence to df
    nc <- length(sample_metrics) # number of columns
    expMat <- cbind(expMat, "Prop" = expMat$Freq / nc) # calculate proportion of studies/batches expressed
    topMat <- cbind(topMat, "Prop" = topMat$Freq / nc) # calculate proportion of studies/batch high-confidence
    sample_metrics[["ExpressionMatrix"]] <- expMat # store expression matrix for saving
    sample_metrics[["TopMatrix"]] <- topMat # store high confidence matrix for saving
    print("8")
    # get genes which are expressed and high-confidence according to use defined ratio
    sample_metrics[["ExpressedGenes"]] <- as.character(expMat$expressedGenes[expMat$Prop >= batch_ratio])
    sample_metrics[["TopGenes"]] <- as.character(topMat$topGenes[topMat$Prop >= batch_ratio_high])
    print("9")

    # create a table to write gene expression and high confidence to
    write_table <- data.frame(entrez_all)
    write_table <- cbind(write_table, rep(0, nrow(write_table)))
    for (i in seq_len(nrow(write_table))) {
        if (as.character(write_table[i, 1]) %in% as.character(sample_metrics$ExpressedGenes)) {
            write_table[i, 2] <- 1
        }
        if (as.character(write_table$entrez_all[i]) %in% as.character(sample_metrics$TopGenes)) {
            write_table[i, 3] <- 1
        }
    }
    print("10")
    header <- c("ENTREZ_GENE_ID", "expressed", "high")
    #write_table <- rbind(header, write_table)
    colnames(write_table) <- header
    write.csv(write_table, out_file, row.names = FALSE, col.names = FALSE)
    print("11")
    print(head(write_table))
}

save_rnaseq_tests(
  counts_matrix_file = "/home/joshl/projects/COMO/main/data/data_matrices/naiveB/gene_counts_matrix_total_naiveB.csv",
  config_filepath = "/home/joshl/projects/COMO/main/data/config_sheets/trnaseq_data_inputs_auto.xlsx",
  out_file = "/home/joshl/projects/COMO/rnaseq_r_out.csv",
  info_file = "/home/joshl/projects/COMO/main/data/gene_info.csv",
  context_name = "naiveB",
  technique = "zfpkm"
)