suppressPackageStartupMessages(library("stringr"))
suppressPackageStartupMessages(library("dplyr"))
suppressPackageStartupMessages(library("tidyverse"))
suppressPackageStartupMessages(library("ggplot2"))

# Check if rlogs directory exists, From: https://stackoverflow.com/a/46008094
# Then prevent messy messages from repeatedly writing to juypter
work_dir <- getwd()
r_log_directory <- stringr::str_interp("${work_dir}/logs")
if (!dir.exists(r_log_directory)) { dir.create(r_log_directory) }
zz <- file(file.path(r_log_directory, "combine_distributions.Rout"), open = "wt")
sink(zz, type = "message")


get_batch_name <- function(x) {
  basename(x)
  return(substring(basename(x), 1, nchar(basename(x)) - 4))
}

parse_contexts_zfpkm <- function(wd, contexts, prep) {
  
  batches <- list()
  for (context in contexts) {
    files <- Sys.glob(file.path(wd, context, prep, paste0("zFPKM_Matrix_", prep, "_*.csv")))
    batches[[context]] <- unlist(lapply(files, get_batch_name))
  }
  
  return(batches)
}


parse_contexts_zumi <- function(wd, contexts, prep) {
  
  batches <- list()
  for (context in contexts) {
    files <- Sys.glob(file.path(wd, context, prep, paste0("zUMI_Matrix_", prep, "_*.csv")))
    batches[[context]] <- unlist(lapply(files, get_batch_name))
  }
  
  return(batches)
}


parse_contexts_zscore_prot <- function(wd, contexts) {
  
  batches <- list()
  # files <- Sys.glob(file.path(wd, "*", "proteomics", "protein_zscore_Matrix_*.csv"))
  for (context in contexts) {
    files <- Sys.glob(file.path(wd, context, "proteomics", "protein_zscore_Matrix_*.csv"))
    batches[[context]] <- unlist(lapply(files, get_batch_name))
  }
  
  return(batches)
}


merge_batch <- function(wd, context, batch) {
  files <- Sys.glob(file.path(wd, paste0("*", batch, "*")))
  nrep <- c()
  stopifnot(length(files) > 0)
  
  for (f in files) {
    zmat <- read.table(f, strip.white = T, header = T, sep = ",", row.names = NULL) %>% # read expression matrix
      mutate(across(colnames(.)[-1], as.numeric)) %>% # ensure expression values are numbers
      mutate(across(colnames(.)[1], as.character)) %>% # ensure ENTREZ_GENE_IDs are character
      group_by(ENTREZ_GENE_ID) %>%
      summarise_each(funs(max)) %>% # if multiple of same ENTREZ_GENE_ID, take max value
      na.omit(.) %>%
      as.data.frame(.)
    
    nrep <- c(nrep, ncol(zmat) - 1)
    entrez_gene <- zmat[, "ENTREZ_GENE_ID"]
    rep_names <- colnames(zmat)
    zmat <- do.call(cbind, lapply(2:ncol(zmat), function(j) {
      repz <- zmat[, j]
    })) %>%
      cbind(as.character(entrez_gene), .) %>%
      as.data.frame(.) %>%
      na.omit(.)
    
    colnames(zmat) <- rep_names
    
    stack_df <- do.call(rbind, lapply(2:ncol(zmat), function(j) {
      repz <- as.numeric(as.character(zmat[, j]))
      cbind(ENTREZ_GENE_ID = zmat[, "ENTREZ_GENE_ID"], zscore = repz, source = rep(colnames(zmat)[j], length(repz)))
    })) %>% as.data.frame(.)
    
    stack_df$zscore <- as.numeric(as.character(stack_df$zscore))
    
    plot_name_pdf <- file.path(
      wd,
      "figures",
      paste0(
        "plot_", context, "_",
        substring(
          basename(f),
          1,
          nchar(basename(f)) - 4
        ),
        ".pdf"
      )
    )
    
    plot_name_png <- file.path(
      wd,
      "figures",
      paste0(
        "plot_",
        context,
        "_",
        substring(
          basename(f),
          1,
          nchar(basename(f)) - 4
        ),
        ".png"
      )
    )
    
    
    pdf(plot_name_pdf)
    png(
      plot_name_png,
      res = 1200,
      units = "in",
      width = 3.25,
      height = 3.25,
      type = "cairo"
    )
    # label <- colnames(stack_df)[-1]
    simplified_plot <- ifelse(length(unique(stack_df$source)) > 10, TRUE, FALSE)
    plot <- ggplot(stack_df, aes(zscore, color = source)) +
      geom_density() +
      theme(text = element_text(size = 12, family = "sans"))
    if (simplified_plot) {
      plot <- plot + theme(legend.position = "none")
    }
    max_dens <- 0
    # get y upper limit by finding density of peak at z = 0
    for (source in unique(plot$data$source)) {
      source_z <- plot$data$zscore[plot$data$source == source] %>% .[!is.na(.)]
      source_densx <- density(source_z)$x
      source_densy <- density(source_z)$y
      idx <- min(max(which(source_densx <= 0)), min(which(source_densx == 0)))
      max_d <- source_densy[idx]
      if (max_d > max_dens) {
        max_dens <- max_d
      }
    }
    # plot <- plot + ylim(0, 1.5*max_dens)
    dev.off()
  }
  
  return(list(zmat, nrep))
}


combine_batch_zdistro <- function(wd, context, batch, zmat) {
  plot_name_pdf <- file.path(wd, "figures", paste0("plot_", context, "_", batch, "_combine_distro", ".pdf"))
  plot_name_png <- file.path(wd, "figures", paste0("plot_", context, "_", batch, "_combine_distro", ".png"))
  
  weighted_z <- function(x) {
    floor_score <- -6
    ceil_score <- 6
    x <- as.numeric(x)
    numer <- sum(x)
    denom <- sqrt(length(x))
    result <- numer / denom
    if (result < floor_score) { result <- floor_score }
    if (result > ceil_score) { result <- ceil_score }
    return(result)
  }
  
  if (ncol(zmat) > 2) {
    combine_z <- apply(zmat[, -1], 1, weighted_z)
    merge_df <- cbind(zmat, combined = combine_z)
    combine_z <- cbind(ENTREZ_GENE_ID = as.character(zmat[, "ENTREZ_GENE_ID"]), combine_z)
    
    stack_df <- do.call(rbind, lapply(2:ncol(merge_df), function(j) {
      repz <- as.numeric(as.character(merge_df[, j]))
      cbind(ENTREZ_GENE_ID = merge_df[, "ENTREZ_GENE_ID"], zscore = repz, source = rep(colnames(merge_df)[j], length(repz)))
    })) %>% as.data.frame(.)
    stack_df$zscore <- as.numeric(as.character(stack_df$zscore))
    
    simplified_plot <- ifelse(length(unique(stack_df$source)) > 10, TRUE, FALSE)
    label <- colnames(stack_df)[-1]
    pdf(plot_name_pdf)
    png(
      plot_name_png,
      res = 1200,
      units = "in",
      width = 3.25,
      height = 3.25,
      type = "cairo"
    )
    
    if (simplified_plot) {
      #p <- p + theme(legend.position = "none")
      stack_df <- stack_df[stack_df$source == "combined",]
    }
    
    p <- ggplot(stack_df, aes(zscore, color = source)) +
      geom_density() +
      theme(text = element_text(size = 12, family = "sans"))
    
    max_dens <- 0
    # get y upper limit by finding density of peak at z = 0
    for (source in unique(p$data$source)) {
      source_z <- p$data$zscore[p$data$source == source] %>% .[!is.na(.)] #%>% .[!is.nan(.)]
      source_densx <- density(source_z)$x
      source_densy <- density(source_z)$y
      idx <- min(max(which(source_densx <= 0.5)), min(which(source_densx == 0.5)))
      max_d <- source_densy[idx]
      if (max_d > max_dens) { max_dens <- max_d }
    }
    #p <- p + ylim(0, 1.5*max_dens)
    print(p)
    dev.off()
    
  } else {
    combine_z <- zmat
  }
  
  return(as.data.frame(combine_z))
}


combine_context_zdistro <- function(wd, context, n_reps, zmat) {
  plot_name_pdf <- file.path(wd, "figures", paste0(
    "plot_", context, "_combine_batches_distro", ".pdf"))
  plot_name_png <- file.path(wd, "figures", paste0(
    "plot_", context, "_combine_batches_distro", ".png"))
  
  weighted_z <- function(x, n_reps) {
    floor_score <- -6
    ceil_score <- 6
    x <- as.numeric(x)
    nas <- sort(unique(c(which(is.nan(x)), which(is.na(x)))))
    weights <- c()
    for (i in seq_along(n_reps)) { weights <- c(weights, (n_reps[i]) / sum(n_reps)) }
    if (length(nas) > 0) {
      x <- x[-nas]
      weights <- weights[-nas]
    }
    numer <- sum(weights * x)
    denom <- sqrt(sum(weights^2))
    result <- numer / denom
    if (result < floor_score) { result <- floor_score }
    if (result > ceil_score) { result <- ceil_score }
    return(result)
  }
  
  if (ncol(zmat) > 2) {
    combine_z <- apply(zmat[, -1], 1, weighted_z, n_reps = n_reps)
    merge_df <- cbind(zmat, combined = combine_z)
    combine_z <- cbind(ENTREZ_GENE_ID = as.character(zmat[, "ENTREZ_GENE_ID"]), combine_z)
    
    stack_df <- do.call(rbind, lapply(2:ncol(merge_df), function(j) {
      repz <- as.numeric(as.character(merge_df[, j]))
      cbind(ENTREZ_GENE_ID = merge_df[, 1], zscore = repz, source = rep(colnames(merge_df)[j], length(repz)))
    })) %>% as.data.frame(.)
    stack_df$zscore <- as.numeric(as.character(stack_df$zscore))
    
    label <- colnames(stack_df)[-1]
    pdf(plot_name_pdf)
    png(
      plot_name_png,
      res = 1200,
      units = "in",
      width = 3.25,
      height = 3.25,
      type = "cairo"
    )
    p <- ggplot(stack_df, aes(zscore, color = source)) +
      geom_density() +
      theme(text = element_text(size = 12, family = "sans"))
    
    max_dens <- 0
    # get y upper limit by finding density of peak at z = 0
    for (source in unique(p$data$source)) {
      source_z <- p$data$zscore[p$data$source == source] %>% .[!is.na(.)] #%>% .[!is.nan(.)]
      source_densx <- density(source_z)$x
      source_densy <- density(source_z)$y
      idx <- min(max(which(source_densx <= 0)), min(which(source_densx == 0)))
      max_d <- source_densy[idx]
      if (max_d > max_dens) { max_dens <- max_d }
    }
    #p <- p + ylim(0, 1.5*max_dens)
    print(p)
    dev.off()
    
  } else {
    combine_z <- zmat
    colnames(combine_z) <- c("ENTREZ_GENE_ID", "combine_z")
  }
  
  return(as.data.frame(combine_z))
}


combine_omics_zdistros <- function(
  wd,
  context,
  comb_batches_z_trna,
  comb_batches_z_mrna,
  comb_batches_z_scrna,
  comb_batches_z_prot,
  tweight,
  mweight,
  sweight,
  pweight,
  keep_gene_scores = TRUE) {
  
  
  fig_path <- file.path(wd, context, "figures")
  if (!file.exists(fig_path)) {
    dir.create(fig_path, recursive = TRUE)
  }
  plot_name_pdf <- file.path(fig_path, paste0("plot_", context, "_combine_omics_distro", ".pdf"))
  plot_name_png <- file.path(fig_path, paste0("plot_", context, "_combine_omics_distro", ".png"))
  
  weights <- c()
  names <- c()
  dfs <- list()
  counter <- 0
  if (tweight > 0) {
    counter <- counter + 1
    weights <- c(weights, tweight)
    names <- c(names, "total")
    dfs[[counter]] <- comb_batches_z_trna
  }
  if (mweight > 0) {
    counter <- counter + 1
    weights <- c(weights, mweight)
    names <- c(names, "polyA")
    dfs[[counter]] <- comb_batches_z_mrna
  }
  if (sweight > 0) {
    counter <- counter + 1
    weights <- c(weights, sweight)
    names <- c(names, "singleCell")
    dfs[[counter]] <- comb_batches_z_scrna
  }
  if (pweight > 0) {
    counter <- counter + 1
    weights <- c(weights, pweight)
    names <- c(names, "proteome")
    dfs[[counter]] <- comb_batches_z_prot
  }
  
  weighted_z <- function(x, weights) {
    floor_score <- -6
    ceil_score <- 10
    x <- as.numeric(x)
    
    nas <- which(is.na(x))
    if (length(nas) > 0) {
      x <- x[-nas]
      weights <- weights[-nas]
    }
    weights <- weights / sum(weights)
    numer = sum(weights * x)
    denom = sqrt(sum(weights^2))
    result <- numer / denom
    if (result < floor_score) { result <- floor_score }
    if (result > ceil_score) { result <- ceil_score }
    return(result)
  }
  
  for (i in 1:counter) {
    add_df <- dfs[[i]]
    colnames(add_df)[2] <- names[i]
    if (i == 1) { zmat <- add_df }
    else { zmat <- full_join(zmat, add_df, by = "ENTREZ_GENE_ID", copy = TRUE) }
  }
  
  if (ncol(zmat) > 2) {
    combine_z <- apply(zmat[, -1], 1, weighted_z, weights = weights)
  } else {
    combine_z = zmat[, -1]
  }
  
  merge_df <- cbind(zmat, combined = combine_z)
  combine_z <- cbind(ENTREZ_GENE_ID = as.character(zmat[, "ENTREZ_GENE_ID"]), combine_z)
  
  stack_df <- do.call(rbind, lapply(2:ncol(merge_df), function(j) {
    repz <- as.numeric(as.character(merge_df[, j]))
    cbind(ENTREZ_GENE_ID = merge_df[, 1], zscore = repz, source = rep(colnames(merge_df)[j], length(repz)))
  })) %>% as.data.frame(.)
  stack_df$zscore <- as.numeric(as.character(stack_df$zscore))
  
  label <- colnames(stack_df)[-1]
  pdf(plot_name_pdf)
  png(
    plot_name_png,
    res = 1200,
    units = "in",
    width = 3.25,
    height = 3.25,
    type = "cairo"
  )
  
  p <- ggplot(stack_df, aes(zscore, color = source)) +
    geom_density() +
    theme(text = element_text(size = 12, family = "sans"))
  
  max_dens <- 0
  # get y upper limit by finding density of peak at z = 0
  for (source in unique(p$data$source)) {
    source_z <- p$data$zscore[p$data$source == source] %>% .[!is.na(.)] #%>% .[!is.nan(.)]
    source_densx <- density(source_z)$x
    source_densy <- density(source_z)$y
    idx <- min(max(which(source_densx <= 0)), min(which(source_densx == 0)))
    max_d <- source_densy[idx]
    if (max_d > max_dens) { max_dens <- max_d }
  }
  #p <- p + ylim(0, 1.5*max_dens)
  print(p)
  dev.off()
  
  #if ( keep_gene_scores ) {
  #    combine_z <- rbind(combine_z, comb_batches_z_rna[!(comb_batches_z_rna$ENTREZ_GENE_ID %in% combine_z$ENTREZ_GENE_ID), "ENTREZ_GENE_ID"])
  #}
  
  return(combine_z)
}


combine_zscores_main <- function(
  working_dir,
  context_names,
  global_use_mrna,
  global_use_trna,
  global_use_scrna,
  global_use_proteins,
  keep_gene_scores,
  global_trna_weight,
  global_mrna_weight,
  global_scrna_weight,
  global_protein_weight
) {
  
  figure_output_dir <- file.path(working_dir, "figures")
  if (!file.exists(figure_output_dir)) { dir.create(figure_output_dir) }
  
  global_trna_batches <- parse_contexts_zfpkm(working_dir, context_names, "total")
  global_mrna_batches <- parse_contexts_zfpkm(working_dir, context_names, "mrna")
  global_scrna_batches <- parse_contexts_zumi(working_dir, context_names, "scrna")
  global_protein_batches <- parse_contexts_zscore_prot(working_dir, context_names)
  
  for (context in context_names) {
    context_use_trna <- global_use_trna
    context_use_mrna <- global_use_mrna
    context_use_scrna <- global_use_scrna
    context_use_proteins <- global_use_proteins
    
    context_trna_weight <- global_trna_weight
    context_mrna_weight <- global_mrna_weight
    context_scrna_weight <- global_scrna_weight
    context_protein_weight <- global_protein_weight
    
    context_trna_batch <- global_trna_batches[[context]]
    context_mrna_batch <- global_mrna_batches[[context]]
    context_scrna_batch <- global_scrna_batches[[context]]
    context_protein_batch <- global_protein_batches[[context]]
    
    
    if (length(context_trna_batch) == 0 & global_use_trna) {
      context_use_trna <- FALSE
      print(paste0("No total RNA-seq zFPKM Matrix files found for ", context, ". Will not use for this context."))
    }
    
    if (length(context_mrna_batch) == 0 & global_use_mrna) {
      context_use_mrna <- FALSE
      print(paste0("No polyA RNA-seq zFPKM Matrix files found for ", context, ". Will not use for this context."))
    }
    
    if (length(context_scrna_batch) == 0 & global_use_scrna) {
      context_use_scrna <- FALSE
      print(paste0("No SC RNA-seq zFPKM Matrix files found for ", context, ". Will not use for this context."))
    }
    
    if (length(context_protein_batch) == 0 & global_use_proteins) {
      context_use_proteins <- FALSE
      print(paste0("No proteomics z-score Matrix files found for ", context, ". Will not use for this context."))
    }
    
    if (context_use_trna) {
      print("Will merge total RNA-seq distributions")
      trna_workdir <- file.path(working_dir, context, "total")
      num_reps <- c()
      count <- 0
      for (batch in context_trna_batch) {
        res <- merge_batch(trna_workdir, context, batch)
        zmat <- res[[1]]
        num_reps <- c(num_reps, res[[2]])
        comb_z <- combine_batch_zdistro(trna_workdir, context, batch, zmat)
        colnames(comb_z) <- c("ENTREZ_GENE_ID", batch)
        if (!count) { merge_z <- comb_z }
        else { merge_z <- full_join(merge_z, comb_z, by = "ENTREZ_GENE_ID") }
        count <- count + 1
      }
      
      comb_batches_z_trna <- combine_context_zdistro(trna_workdir, context, num_reps, merge_z)
      filename <- file.path(trna_workdir, paste0("combined_zFPKM_", context, ".csv"))
      write.csv(comb_batches_z_trna, filename, row.names = FALSE)
      
      if (!context_use_proteins & !context_use_mrna & !context_use_scrna) {
        filename <- file.path(working_dir, context, "total", paste0("model_scores_", context, ".csv"))
        write.csv(comb_batches_z_trna, filename, row.names = FALSE)
      }
      
    } else { comb_batches_z_trna <- NA }
    
    
    if (context_use_mrna) {
      print("Will merge polyA enriched RNA-seq distributions")
      mrna_workdir <- file.path(working_dir, context, "mrna")
      num_reps <- c()
      count <- 0
      for (batch in context_mrna_batch) {
        res <- merge_batch(mrna_workdir, context, batch)
        zmat <- res[[1]]
        num_reps <- c(num_reps, res[[2]])
        comb_z <- combine_batch_zdistro(mrna_workdir, context, batch, zmat)
        colnames(comb_z) <- c("ENTREZ_GENE_ID", batch)
        if (!count) { merge_z <- comb_z }
        else { merge_z <- full_join(merge_z, comb_z, by = "ENTREZ_GENE_ID") }
        count <- count + 1
      }
      
      comb_batches_z_mrna <- combine_context_zdistro(mrna_workdir, context, num_reps, merge_z)
      filename <- file.path(mrna_workdir, paste0("combined_zFPKM_", context, ".csv"))
      write.csv(comb_batches_z_mrna, filename, row.names = FALSE)
      
      if (!context_use_proteins & !context_use_trna & !context_use_scrna) {
        filename <- file.path(mrna_workdir, paste0("model_scores_", context, ".csv"))
        write.csv(comb_batches_z_mrna, filename, row.names = FALSE)
      }
      
    } else { comb_batches_z_mrna <- NA }
    
    
    if (context_use_scrna) {
      print("Will merge single-cell RNA-seq distributions")
      scrna_workdir <- file.path(working_dir, context, "scrna")
      num_reps <- c()
      count <- 0
      for (batch in context_scrna_batch) {
        res <- merge_batch(scrna_workdir, context, batch)
        zmat <- res[[1]]
        num_reps <- c(num_reps, res[[2]])
        comb_z <- combine_batch_zdistro(scrna_workdir, context, batch, zmat)
        colnames(comb_z) <- c("ENTREZ_GENE_ID", batch)
        if (!count) { merge_z <- comb_z }
        else { merge_z <- full_join(merge_z, comb_z, by = "ENTREZ_GENE_ID") }
        count <- count + 1
      }
      
      comb_batches_z_scrna <- combine_context_zdistro(scrna_workdir, context, num_reps, merge_z)
      filename <- file.path(scrna_workdir, paste0("combined_zFPKM_", context, ".csv"))
      write.csv(comb_batches_z_scrna, filename, row.names = FALSE)
      
      if (!context_use_proteins & !context_use_trna & !context_use_mrna) {
        filename <- file.path(scrna_workdir, paste0("model_scores_", context, ".csv"))
        write.csv(comb_batches_z_scrna, filename, row.names = FALSE)
      }
      
    } else { comb_batches_z_scrna <- NA }
    
    
    if (context_use_proteins) {
      print("Will merge protein abundance distributions")
      protein_workdir <- file.path(working_dir, context, "proteomics")
      num_reps <- c()
      count <- 0
      for (batch in context_protein_batch) {
        res <- merge_batch(protein_workdir, context, batch)
        zmat <- res[[1]]
        num_reps <- c(num_reps, res[[2]])
        comb_z <- combine_batch_zdistro(protein_workdir, context, batch, zmat)
        colnames(comb_z) <- c("ENTREZ_GENE_ID", batch)
        if (!count) { merge_z <- comb_z }
        else { merge_z <- full_join(merge_z, comb_z, by = "ENTREZ_GENE_ID") }
        count <- count + 1
      }
      
      comb_batches_z_prot <- combine_context_zdistro(protein_workdir, context, num_reps, merge_z)
      filename <- file.path(protein_workdir, paste0("combined_zscore_proteinAbundance_", context, ".csv"))
      write.csv(comb_batches_z_prot, filename, row.names = FALSE)
      
      if (!context_use_mrna & !context_use_trna & !context_use_scrna) {
        filename <- file.path(protein_workdir, paste0("model_scores_", context, ".csv"))
        write.csv(comb_batches_z_prot, filename, row.names = FALSE)
      }
      
    } else { comb_batches_z_prot <- NA }
    
    if (!context_use_trna) { context_trna_weight <- 0 }
    if (!context_use_mrna) { context_mrna_weight <- 0 }
    if (!context_use_scrna) { context_scrna_weight <- 0 }
    if (!context_use_proteins) { context_protein_weight <- 0 }
    
    comb_omics_z <- combine_omics_zdistros(
      working_dir,
      context,
      comb_batches_z_trna,
      comb_batches_z_mrna,
      comb_batches_z_scrna,
      comb_batches_z_prot,
      context_trna_weight,
      context_mrna_weight,
      context_scrna_weight,
      context_protein_weight
    )
    
    filename <- file.path(working_dir, context, paste0("model_scores_", context, ".csv"))
    write.csv(comb_omics_z, filename, row.names = FALSE)
  }
  
}