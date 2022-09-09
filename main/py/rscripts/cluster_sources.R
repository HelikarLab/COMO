# This file will cluster transcriptomic and proteomic sources
# It will perform the following functions
#   1. reads the zFPKMs (or TPMs or CPMs) for each replicate
#               *binarizes them if given that arg (see RNAseq.R using the kOverA function)
#   2. clusters all replicates for all batches/studies across all contexts/cell types within one data source at a time
#   3. clusters all replicates for all batches/studies across all data sources within one cell type at a time
#   4. clusters all batches/studies across all contexts/cell types within one data source at a time
#   5. clusters all batches/studies across all data sources within one cell type at a time (edited)

library(stringr)
library(tools)

get_study_value <- function(file_path) {
    # This function will get the S## value from the file path
    
    # Find "S##.csv", but use a lookahead to avoid "taking" the ".csv" portion
    match <- stringr::str_extract(string = toupper(file_path), pattern = "S\\d{1,2}(?=\\.CSV)")

    if (!is.na(match))
        match <- toupper(match)

    return(match)
}

get_replicate_files <- function(results_directory, context_names, source_type, use_trna, use_mrna) {

    # This function will get the file paths of zFPKMs for each replicate
    all_context_files <- list()
     for (context_name in context_names) {
        lower_context_name <- tolower(context_name)
        source_type <- tolower(source_type)
        
        current_context_files <- list.files(file.path(results_directory, lower_context_name), full.names = TRUE, recursive = TRUE)
        context_files <- c()
        
        for (file in current_context_files) {
            file_name <- tolower(file)
            
            # Check the current file meets our criteria
            if (
              tools::file_ext(file_name) == "csv" &&                      # Ensure the current file is a CSV file
              grepl(lower_context_name, file_name) &&  # Test if the current file is part of the current context (i.e., naiveB, immNK)
              grepl(source_type, file_name) &&   # Test if the current file has source type (zFPKM, TPM, CPM)
              # Create a "group" that works if either trna or mrna is TRUE
              (
                (use_trna && grepl("total", file_name)) ||         # Test if the current file is a total-rna file
                (use_mrna && grepl("mrna", file_name))             # Test if the current file is an mRNA (polyA) file
              )
            ) {
                context_files <- append(context_files, file)
            }
        }
        
        # Only append new list if context_files has at least one item
        if (length(context_files) > 0)
            all_context_files[[context_name]] <- context_files
    }
    
    # Return list if it has at least one item, otherwise return "NA"
    if (length(all_context_files) > 0) {
        return(all_context_files)
    } else {
        return(NA)
    }
}

read_matrix_values <- function(study_files) {
    # This function is responsible for reading in the matrix values found within the replicate files
    # It takes the list of replicate files and returns a list of lists of matrix values
    replicate_dataframes <- list()
    context_names <- names(study_files)

    for (context in context_names) {
        index <- 1
        context_files <- study_files[[context]]
        context_dataframe <- c()
        for (file in context_files) {
            dataframe <- read.csv(file, header = TRUE)

            study <- get_study_value(file)
            context_dataframe[[study]] <- dataframe

            index <- index + 1
        }

        if (length(context_dataframe) > 0) {
            replicate_dataframes[[context]] <- context_dataframe
        }
    }
    if (length(replicate_dataframes) > 0) {
        return(replicate_dataframes)
    } else {
        return(NA)
    }
}


cluster_sources_main <- function(
  results_directory,
  context_names,
  source_type,
  use_trna,
  use_mrna,
  binarize_data
) {
    
    study_files <- get_replicate_files(results_directory = results_directory, context_names = context_names,  source_type = source_type, use_trna = use_trna, use_mrna = use_mrna)
    study_dataframes <- read_matrix_values(study_files = study_files)
    print("DONE")
}


results_directory <- "/Users/joshl/docker/madrid/local_files/results"
context_names <- list("immNK", "naiveB")
source_type <- "zFPKM"
use_trna <- TRUE
use_mrna <- TRUE
binarize_data <- FALSE
cluster_sources_main(results_directory = results_directory, context_names = context_names, source_type = source_type, use_trna = use_trna, use_mrna = use_mrna, binarize_data = binarize_data)
