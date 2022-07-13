# Check if rlogs directory exists, From: https://stackoverflow.com/a/46008094
library("stringr")
username <- Sys.info()["user"]
work_dir <- str_interp("/home/${username}/work")

if (!dir.exists(str_interp("${work_dir}/py/rlogs"))) {
    dir.create(str_interp("${work_dir}/py/rlogs"))
}

# prevent messy messages from repeatedly writing to juypter

zz <- file(file.path("/home", username, "work", "py", "rlogs", "generate_counts_matrix.Rout"), open="wt")
sink(zz, type="message")

library(tidyverse)

# fetch and organize MADRID_input
organize_gene_counts_files <- function(data_dir) {
    # fetch and organize MADRID_input
    # accepts path to data directory, normalization technique

    study_metrics <- list() # list storing all the samples
    tissue_name <- basename(data_dir) # get tissue name from directory name

    # collect paths to  files
    count_dir_i <- file.path(data_dir, "geneCounts") # MADRID_inputs/tissueName/geneCounts
    count_dir <- list.dirs(path = count_dir_i, full.names = TRUE, recursive = FALSE) %>%
      Filter(function(x) !any(grepl(".ipynb_checkpoints", x)), .) # MADRID_inputs/tissueName/geneCounts/SX

    strand_dir_i <- file.path(data_dir, "strandedness") # MADRID_inputs/tissueName/strandedness
    strand_dir <- list.dirs(path = strand_dir_i, full.names = TRUE, recursive = FALSE) # tissueName/strandedness/SX

    # for each study, collect gene count files, frag files, insert size files, layouts, and strand info
    for ( j in 1:length(count_dir) ) {
        cnt_d <- count_dir[j] # study's gene count folder
        sname <- basename(cnt_d) # get study name from directory name
        entry <- list() # list specific to a study, top level of SampleMetrics
        str_d <- file.path(data_dir, "strandedness", sname) # strandedness file directory

        #counts_glob <- paste(c(cnt_d, "/*.tab"), collapse="")
        counts_glob <- file.path(cnt_d, "*.tab")
        counts_files <- Sys.glob(counts_glob) # outputs of STAR

        n_replicates <- length(counts_files) # number of replicates for this study
        replicate_names <- rep(0,n_replicates) # initialize empty array to store rep names
        strandedness_files <- rep(0, n_replicates) # initialize empty array to store

        for ( i in 1:n_replicates ) { # for each replicate, get name of replicate ex. S1R1
            replicate_file <- str_match(counts_files[i], "geneCounts/\\s*(.*?)\\s*.tab")[,2] # file path
            replicate_names[i] <- basename(replicate_file) # file name only
            rname <- sub(pattern = "(.*)\\..*$", replacement = "\\1", basename(replicate_file)) # rem file extension
            strandedness_file <- paste0(rname, "_strandedness.txt") # use rep name to get matching strand
            strandedness_files[i] <- file.path(str_d, strandedness_file) # path to strand files

        }

        entry[["CountFiles"]] <- counts_files # assign gene count file paths to study entry
        entry[["NumSamples"]] <- n_replicates # assign number of replicates to study entry
        entry[["SampleNames"]] <- replicate_names # assign sample names to study entry
        entry[["StrandFiles"]] <- strandedness_files # assign strandedness files
        study_metrics[[sname]] <- entry # assign study entry to the list of studies

      }

      return(study_metrics)
}


prepare_sample_counts <- function(counts_file, counts_files, strand_file) {
    
    # prepare a gene count file to be added to the count matrix
    if ( grepl("R\\d+r1", counts_file, ignore.case=FALSE) ) { # if first run in a multirun sample
        rm_run_number <- unlist(strsplit(counts_file, "r1"))[1] # remove run number from sample id
        #split_directory <- unlist(strsplit(rm_run_number, "/")) # split from path
        #basename <- split_directory[length(split_directory)] # get sample id
        bname <- basename(rm_run_number)
        search_str <- paste0(bname, "r\\d+")
        run_files <- counts_files[grepl(search_str, counts_files)] # search for other files that are part of same run
        samp_count <- NULL

        
        for ( f in run_files) { # for each run associated with this replicate,
            run_count <- read.delim(f) # read gene counts
            run_len <- length(run_count[,1]) # number of genes
            run_count <- run_count[4:run_len,] # only count mapped genes
            run_count <- na.omit(run_count) # omit na if applicable
            genes <- run_count[,1] # list of genes

            strand_read <- gsub("[\r\n]", "", readChar(strand_file, file.info(strand_file)$size))

            stopifnot( (strand_read=="NONE") ||
                       (strand_read=="SECOND_READ_TRANSCRIPTION_STRAND") ||
                       (strand_read=="FIRST_READ_TRANSCRIPTION_STRAND")
                     )
            
            if ( strand_read=="FIRST_READ_TRANSCRIPTION_STRAND" ) { # forward
                cl <- 3
            } else if ( strand_read=="SECOND_READ_TRANSCRIPTION_STRAND" ) { # reverse
                cl <- 4
            } else { # unstranded
                cl <- 2
            }

            #cl <- which.max(colSums(run_count[,2:4])) + 1 # max counts is most likely the correct strandedness

            run_count <- data.frame(cbind(run_count[,1], run_count[,cl])) # get counts from correct strandedness
            run_count[,1] <- genes # force rewrite genes bc rpy2 does not handle previous step properly
            colnames(run_count) <- c("genes", "counts") # rename columns to descriptive headers

            if ( is.null(samp_count) ) { # if haven't declared sample counts yet
                samp_count <- run_count # assign samp_counts to this run count if first one
            }

            else { # otherwise merge
                samp_count <- merge(samp_count, run_count, by="genes", all=TRUE)
                samp_count[is.na(samp_count)] <- 0
            }
        }

        genes <- samp_count["genes"]
        samp_count["genes"] <- NULL
        samp_count <- sapply(samp_count, as.numeric)
        samp_sum <- rowSums(samp_count) # sum the counts from each run
        samp_count <- data.frame(cbind(genes, samp_sum))
        samp_count[,1] <- genes # rewrite bc rpy2 is weird
        colnames(samp_count)[2] <- "counts"

        return(samp_count)
    }

    else if ( grepl("R\\d+r", counts_file)==TRUE ) { # multirun files handled when the first run is iterated on, skip
        return("skip")
    }

    else { # if not multirun, handle normally
        samp_count <- read.delim(counts_file) # read file
        samp_len <- length(samp_count[,1]) # number of genes
        samp_count <- samp_count[4:samp_len,] # only count mapped genes
        samp_count <- na.omit(samp_count) # remove NAs
        genes <- samp_count[,1] # get gene names

        strand_read <- gsub("[\r\n]", "", readChar(strand_file, file.info(strand_file)$size))
        
        stopifnot( (strand_read=="NONE") ||
                   (strand_read=="SECOND_READ_TRANSCRIPTION_STRAND") ||
                   (strand_read=="FIRST_READ_TRANSCRIPTION_STRAND")
                 )

        if ( strand_read=="FIRST_READ_TRANSCRIPTION_STRAND" ) { # forward
            cl <- 3
        } else if ( strand_read=="SECOND_READ_TRANSCRIPTION_STRAND" ) { # reverse
            cl <- 4
        } else { # unstranded
            cl <- 2
        }
        #cl <- which.max(colSums(samp_count[,2:4])) + 1 # max counts is most likely the correct strandedness

        samp_count <- data.frame(cbind(samp_count[,1], samp_count[,cl])) # df from gene names and correct counts
        samp_count[,1] <- genes # force rewrite genes bc rpy2 does not handle previous step properly
        colnames(samp_count) <- c("genes", "counts") # rename columns for merging

        return(samp_count)
    }
}


create_counts_matrix <- function(counts_files, replicate_names, n_replicates, strand_files) {
    i_adjust <- 0  # adjusted index, subtracts number of multiruns processed
    counts <- prepare_sample_counts(counts_files[1], counts_files, strand_files[1]) # get first column of counts to add
    
    if ( grepl("R\\d+r1", replicate_names[1], ignore.case=FALSE) ) { # if first of a set of multiruns
        colnames(counts)[2] <- unlist(strsplit(replicate_names[1], "r\\d+"))[1] # remove run number tag

    } else {
        colnames(counts)[2] <- replicate_names[1] # if not multirun, header is okay
    }

    
    for ( i in 2:n_replicates) { # iterate through rest of replicates after handling the first and initializing matrix
        new_cnts <- prepare_sample_counts(counts_files[i], counts_files, strand_files[i]) # get next column of counts to add
        options(warn=-1) # turn off warnings
        
        if ( new_cnts=="skip" ) { # handle skipped multirun calls
            i_adjust <- i_adjust+1
            next
        }

        counts <- merge(counts, new_cnts, by="genes", all=TRUE) # merge new count list with matrix
        counts[is.na(counts)] <- 0 # if gene lists differ (could be the case if aligment platforms differ) make na's 0

        if ( grepl("R\\d+r1", replicate_names[i], ignore.case=FALSE) ) { # remove run number from multirun name
            samp_name <- unlist(strsplit(replicate_names[i], "r\\d+"))[1]
        }
        else {
            samp_name <- replicate_names[i] # if not multirun, sample name is okay
        }

        colnames(counts)[i+1-i_adjust] <- samp_name # set column name of added counts

    }

    return(counts)
}


generate_counts_matrix_main <- function(data_dir, out_dir) {

    print("Organizing Files")
    study_metrics <- organize_gene_counts_files(data_dir)
    print("Creating counts matrix")

    for ( i in 1:length(study_metrics) ) { # for each study
        res <- create_counts_matrix(study_metrics[[i]][["CountFiles"]],
                             study_metrics[[i]][["SampleNames"]],
                             study_metrics[[i]][["NumSamples"]],
                             study_metrics[[i]][["StrandFiles"]]
                             ) # create count matrix for study

        study_metrics[[i]][["CountMatrix"]] <- res # store count matrix

        study_metrics[[i]][["NumSamples"]] <- ncol(study_metrics[[i]][["CountMatrix"]]) # store number of samples

        if ( i == 1 ) { # if first study
            full_count_matrix <- study_metrics[[i]][["CountMatrix"]] # initialize supermatrix
        } else { # for successive studies
            add_mat <- study_metrics[[i]][["CountMatrix"]]
            full_count_matrix <- merge(full_count_matrix, add_mat, by="genes", all=TRUE) # merge study matrix to super
        }

    }

    full_count_matrix["genes"] <- sub("\\.\\d+", "", full_count_matrix$genes)
    file_split <- unlist(strsplit(data_dir, "/"))
    file_name <- file.path(out_dir, paste0("gene_counts_matrix_full_", file_split[length(file_split)], ".csv"))
    write.csv(full_count_matrix, file_name, row.names=FALSE)
    cat("Count Matrix written at ", file_name, "\n")
}
