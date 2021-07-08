#!/usr/bin/python3

import os, time, sys
import pandas as pd
import getopt
from rpy2.robjects.packages import importr, SignatureTranslatedAnonymousPackage
#from rpy2.robjects import r, pandas2ri
#import rpy2.robjects as ro
#from rpy2.robjects.conversion import localconverter
#from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage

#pandas2ri.activate()

#limma = importr("limma")
tidyverse = importr("tidyverse")
#edgeR = importr("edgeR")
#genefilter = importr("genefilter")
#biomaRt = importr("biomaRt")
#sjmisc = importr("sjmisc")

# automatically convert ryp2 dataframe to Pandas dataframe
string="""
library(tidyverse)
organizeFiles <- function(data_dir, technique) {
  SampMetrics <- list()
  countFiles <- list()
  fragmentFiles <- list()
  N_samples <- list()
  Names <- list()
  count_dir_i <- paste(c(data_dir, "/geneCounts"),collapse="")
  tab_dir <- list.dirs(path = count_dir_i, full.names = TRUE, recursive = FALSE)
  frag_dir_i <- paste(c(data_dir, "/fragLengths"),collapse="")
  frag_dir <- list.dirs(path = frag_dir_i, full.names = TRUE, recursive = FALSE)
  for (j in 1:length(tab_dir) ) {
    s <- tab_dir[j]
    f <- frag_dir[j]
    sname <- unlist(strsplit(s,"/"))
    sname <- sname[length(sname)]
    entry <- list()
    if ( technique=="zFPKM" ) {
      fragGlob <- paste(c(f, "/*.frag"), collapse="")
      fragFiles <- Sys.glob(fragGlob)
      entry[["FragmentFiles"]] <- fragFiles
    }
    else {
      entry[["FragmentFiles"]] <- NA
    }
    fileGlob <- paste(c(s, "/*.tab"), collapse="")
    cntFiles <- Sys.glob(fileGlob)
    n_samples <- length(cntFiles)
    sample_names <- rep(0,n_samples)
    for (i in 1:n_samples) {
      samp_file <- str_match(cntFiles[i], "bulk_geneCounts/\\\\s*(.*?)\\\\s*.tab")[,2]
      sample_names[i] <- unlist(strsplit(samp_file,"/"))[2]
    }
    entry[["CountFiles"]] <- cntFiles
    entry[["NumSamples"]] <- n_samples
    entry[["SampleNames"]] <- sample_names
    SampMetrics[[sname]] <- entry
  }
  return(SampMetrics)
}

prepSampCnts <- function(cntFile,files) {
  if ( grepl("R\\\\d+r1", cntFile, ignore.case=FALSE) ) {
    basename <- unlist(strsplit(cntFile, "r1"))[1]
    basename <- unlist(strsplit(basename, "/"))
    basename <- basename[length(basename)]
    search_str <- paste(c(basename, "r\\\\d+"),collapse="")
    run_files <- files[grepl(search_str, files)]
    samp_count <- NULL
    for ( f in run_files) {
      run_count <- read.delim(f)
      run_len <- length(run_count[,1])
      run_count <- run_count[4:run_len,]
      run_count <- na.omit(run_count)
      cl <- which.max(colSums(run_count[,2:4])) + 1
      run_count <- data.frame(cbind(run_count[,1], run_count[,cl]))
      colnames(run_count) <- c("genes", "counts")
      if ( is.null(samp_count) ) {
        samp_count <- run_count
      }
      else {
        samp_count <- merge(samp_count, run_count, by="genes", all=TRUE)
        samp_count[is.na(samp_count)] <- 0
      }
    }
    rnames <- samp_count["genes"]
    samp_count["genes"] <- NULL
    samp_count <- sapply(samp_count, as.numeric)
    samp_sum <- rowSums(samp_count)
    samp_count <- data.frame(cbind(rnames, samp_sum))
    colnames(samp_count)[2] <- "counts"
    return(samp_count)
  }
  else if ( grepl("R\\\\d+r", cntFile)==TRUE ) {
    return("skip")
  }
  else {
    samp_count <- read.delim(cntFile)
    samp_len <- length(samp_count[,1])
    samp_count <- samp_count[4:samp_len,]
    samp_count <- na.omit(samp_count)
    cl <- which.max(colSums(samp_count[,2:4])) + 1
    samp_count <- data.frame(cbind(samp_count[,1], samp_count[,cl]))
    colnames(samp_count) <- c("genes", "counts")
    return(samp_count)
  }
}

createCountMatrix <- function(files, sample_names, n_samples,
                              fragFiles, insertSizes, technique) {
  if ( technique=="zFPKM" ) {
    if ( grepl("r1",fragFiles[1],ignore.case=FALSE ) ) {
      groupSize <- c(fragFiles[1])
    }
    else {
      groupSize <- c()
      re_insertSizes <- c(insertSizes[1])
    }
  }
  i_adjust <- 0
  counts <- prepSampCnts(files[1])
  colnames(counts)[2] <- sample_names[1]
  for ( i in 2:n_samples) {
    new_cnts <- prepSampCnts(files[i], files)
    options(warn=-1)
    if ( new_cnts=="skip" ) {
      i_adjust <- i_adjust+1
      if ( technique=="zFPKM" ) {
        groupSize <- c(groupSize, insertSizes[i])
      }
      next
    }
    counts <- merge(counts, new_cnts, by="genes", all=TRUE)
    counts[is.na(counts)] <- 0
    if ( grepl("R\\\\d+r1", sample_names[i], ignore.case=FALSE) ) {
      samp_name <- unlist(strsplit(sample_names[i], "r"))[1]
      if ( technique=="zFPKM" ) {
        if ( length(groupSize) > 1 ) {
          re_insertSizes[i-i_adjust] <- mean(groupSize)
          groupSize <- c(insertSizes[i])
        }
        else {
          re_insertSizes[i-i_adjust] <- insertSizes[i]
        }
      }
    }
    else {
      samp_name <- sample_names[i]
      if ( technique=="zFPKM" ) {
        if ( length(groupSize) > 1 ) {
          re_insertSizes[i_adjust] <- mean(groupSize)
          groupSize <- c()
        }
        else {
          re_insertSizes[i-i_adjust] <- insertSizes[i]
        }
      }
    }
    colnames(counts)[i+1-i_adjust] <- samp_name
  }
  if ( technique=="zFPKM" ) {
    if ( length(groupSize) > 1 ) {
      re_insertSizes[i-i_adjust] <- mean(groupSize)
      groupSize <- c()
    }
  }
  rnames <- counts[,"genes"]
  counts <- counts[,-1]
  counts <- sapply(counts, as.numeric)
  counts <- as.matrix(counts)
  row.names(counts) <- rnames
  if ( technique=="zFPKM" ) {
    res <-list("counts"=counts, "insertSizes"=re_insertSizes)
  }
  else {
    res <- counts
  }
  return(res)
}

genCountMatrix_main <- function(data_dir, out_dir, technique="quantile") {
  filt_options <- list()
  if ( exists("N_rep") ) {
    filt_options$N_rep <- N_rep
  } else {
    filt_options$N_rep <- 0.2
  }
  if ( exists("N_samp") ) {
    filt_options$N_samp <- N_samp
  } else {
    filt_options$N_samp <- 0.5
  }
  if ( exists("percentile") ) {
    filt_options$percentile <- percentile
  } else {
    filt_options$percentile <- 0.5
  }
  if ( exists("min_count") ) {
    filt_options$min_count <- min_count
  } else {
    filt_options$min_count <- 1
  }

  SampMetrics <- organizeFiles(data_dir, technique)


  for ( i in 1:length(SampMetrics) ) {
    if ( technique=="zFPKM" ) {
      j<-0
      insertSizes<-rep(0, SampMetrics[[i]][["NumSamples"]])
      for ( file in SampMetrics[[i]][["FragmentFiles"]] ) {
        j<-j+1
        lines <- readLines(file, n= 50)
        read_idx <- grep("METRICS CLASS", lines)
        size <- as.numeric(unlist(strsplit(lines[read_idx+2],"\t"))[6])
        insertSizes[j] <- as.numeric(unlist(strsplit(lines[read_idx+2],"\t"))[6])
      }
      SampMetrics[[i]][["InsertSizes"]] <- insertSizes
    }
    else {
      SampMetrics[[i]][["InsertSizes"]] <- NA
    }
  }
  for ( i in 1:length(SampMetrics) ) {
    res <- createCountMatrix(SampMetrics[[i]][["CountFiles"]],
                             SampMetrics[[i]][["SampleNames"]],
                             SampMetrics[[i]][["NumSamples"]],
                             SampMetrics[[i]][["FragmentFiles"]],
                             SampMetrics[[i]][["InsertSizes"]],
                             technique)
    if ( technique=="zFPKM" ) {
      SampMetrics[[i]][["CountMatrix"]] <- res$counts
      SampMetrics[[i]][["InsertSizes"]] <- res$insertSizes
    }
    else {
      SampMetrics[[i]][["CountMatrix"]] <- res
    }
    SampMetrics[[i]][["NumSamples"]] <- ncol(SampMetrics[[i]][["CountMatrix"]])
    if ( i == 1 ) {
      full_count_matrix <- data.frame(SampMetrics[[i]][["CountMatrix"]])
    } else {
      add_mat <-data.frame(SampMetrics[[i]][["CountMatrix"]])
      full_count_matrix <- merge(full_count_matrix, add_mat,
                               by="row.names", all=TRUE)
      row.names(full_count_matrix) <- full_count_matrix$Row.names
      full_count_matrix["Row.names"] <- NULL
    }
  }
  write.csv(full_count_matrix, paste(c(out_dir, "BulkRNAseqDataMatrix.csv"), collapse=""))
}
"""
genCountMatrixio = SignatureTranslatedAnonymousPackage(string, "genCountMatrixio")

def main(argv):
    try:
        opts, args = getopt.getopt(argv, "hi:o:t:", ["input_dir=", "output_dir=", "technique="])
    except getopt.GetoptError:
        print('python3 proteomics_gen.py -i <input directory> -o <output_directory> -t <technique>')
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            print('python3 proteomics_gen.py -i <input_directory> -o <output_directory> -t <technique>')
            sys.exit()
        elif opt in ("-i", "--input_dir"):
            input_dir = arg
        elif opt in ("-o", "--output_dir"):
            output_dir = arg
        elif opt in ("-t", "--technique"):
            technique = arg
    print('Input directory is "{}"'.format(input_dir))
    print('Output file is "{}"'.format(output_dir))
    print('Active gene determination technique is "{}"'.format(technique))

    genCountMatrixio.genCountMatrix_main(input_dir, output_dir, technique)


if __name__ == "__main__":
    print(sys.argv)
    main(sys.argv[1:])
