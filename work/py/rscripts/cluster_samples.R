# prevent messy messages from repeatedly writing to juypter
zz <- file(file.path("/home", "jovyan", "work", "py", "rlogs", "rnaseq.Rout"), open="wt")
sink(zz, type="message")


library(ggplot2)
library(ggrepel)
library(tidyverse)
#library(FactoMineR)
library(uwot)


make_logical_matrix <- function(wd, technique, context_names) {
  ### organize logical matrix
  files <- c()
  for ( context in context_names ) {
      if (technique=="zfpkm") {
          files <- c(files, Sys.glob(file.path(wd, context, "**" ,"zFPKM_Matrix_*.csv")))
      } else if ( technique=="quantile" ) {
          files <- c(files, Sys.glob(file.path(wd, context, "**", "TPM_Matrix_*.csv")))
      } else if ( technique=="cpm" ) {
          files <- c(files, Sys.glob(file.path(wd, context, "**", "CPM_Matrix_*.csv")))
      } else {
          print("Invalid technique. Must be zfpkm, quantile, or cpm")
          stop()
      }
  }
  cnt <- 0
  for ( f in files ) {
    new_matrix <- read.table(f, strip.white=T, header=T, sep=",", row.names=NULL) %>% # read expression matrix
        mutate(across(colnames(.)[-1], as.numeric)) %>% # ensure expression values are numbers
        mutate(across(colnames(.)[1], as.character)) %>% # ensure entrez IDs are character
        group_by(ent) %>%
        summarise_each(funs(max)) # if multiple of same entrez, take max value

    #if ( "X" %in% colnames(new_matrix) ) {
    #  new_matrix <- new_matrix %>% select(-X)
    #}
    if (!cnt) { merge_matrix <- new_matrix }
    else { merge_matrix <- merge_matrix %>% left_join(new_matrix, by="ent") }
    cnt <- cnt+1
  }
  if (technique=="zfpkm") {
    cutoff = -3
    logical_matrix <- do.call(cbind, lapply(2:ncol(merge_matrix), function(j) {
    as.integer(merge_matrix[,j] >  cutoff)
    })) %>% as.data.frame(.) %>% cbind(merge_matrix["ent"], .) %>% na.omit(.)
  } else if ( technique=="quantile" ) {
    logical_matrix <- do.call(cbind, lapply(2:ncol(merge_matrix), function(j) {
    cutoff = quantile(merge_matrix[,j], prob=1-quantile/100)
    merge_matrix[,j] >  cutoff
    })) %>% as.data.frame(.) %>% cbind(merge_matrix["ent"], .) %>% na.omit(.)
  } else if ( technique=="cpm" ) {
    logical_matrix <- do.call(cbind, lapply(2:ncol(merge_matrix), function(j) {
    cutoff = ifelse(min_count=="default",
                    10e6/(median(sum(merge_matrix[,j]))),
                    1e6*min_count/(median(sum(merge_matrix[,j]))) )
    merge_matrix[,j] >  cutoff
    })) %>% as.data.frame(a.) %>% cbind(merge_matrix["ent"], .) %>% na.omit(.)
  }
  colnames(logical_matrix) <- colnames(merge_matrix)
  rsums <- rowSums(logical_matrix[,-1])
  logical_matrix <- logical_matrix %>%
    .[rowSums(.[,-1]) < ncol(.)-1,] %>% .[rowSums(.[,-1]) > 0,]

  logical_matrix <- t(logical_matrix[,-1]) # tranpose

  return(logical_matrix)
}


parse_contexts <- function(logical_matrix) {
  contexts <- c()
  batches <- c()
  contexts <- lapply(row.names(logical_matrix), function(r) {
    c(contexts, unlist(strsplit(r, "_S"))[1])
  })
  batches <- lapply(row.names(logical_matrix), function(r) {
    c(batches, unlist(strsplit(r, "R\\d+"))[1])
  })

  contexts <- unlist(contexts)
  batches <- unique(unlist(batches))

  return(list(contexts, batches))
}


# MCA
plot_MCA_replicates <- function(logical_matrix, contexts, wd, label) {
  mca_results <- MCA(logical_matrix, graph=F)
  d <- as.data.frame(mca_results[["ind"]][["coord"]][,1:2])
  d["contexts"] <- as.data.frame(contexts)
  colnames(d) <- c("x", "y", "contexts")
  fig_path <- file.path(wd, "figures")
  if ( !file.exists(fig_path) ) { dir.create(fig_path) }
  plotname <- file.path(fig_path, "mca_plot_replicates.pdf")
  pdf(plotname)
  p <- ggplot(d, ggplot2::aes(x=x, y=y, label=row.names(d), color=contexts)) +
    geom_point(alpha=0.7) +
    geom_text_repel(max.overlaps = Inf) +
    labs(x="Dim 1", y="Dim 2")
    if ( !label ) { p <- remove_geom(p, "GeomTextRepel") }
  print(p)
  dev.off()
}


plot_UMAP_replicates <- function(logical_matrix, contexts, wd, label, n_neigh, min_dist) {
  n_neigh <- ifelse(n_neigh=="default", as.integer(length(contexts)), n_neigh)
  if ( n_neigh<2) {
    print("Cannot cluster replicates if n nearest neighbors is < 1!")
    stop()
  }

  fac_matrix <- do.call(cbind, lapply(1:ncol(logical_matrix), function(n) {
    as.numeric(logical_matrix[,n])
    }))
  coords <- data.frame(umap(fac_matrix, n_neighbors=n_neigh, metric="euclidean",
                            min_dist=min_dist)) %>%
    cbind(., contexts)
  row.names(coords) <- row.names(logical_matrix)
  colnames(coords) <- c("x", "y", "contexts")
  fig_path <- file.path(wd, "figures")
  if ( !file.exists(fig_path) ) { dir.create(fig_path) }
  plotname <- file.path(fig_path, "umap_plot_replicates.pdf")
  pdf(plotname)

  if ( label ) {
      p <- ggplot(coords, ggplot2::aes(x=x, y=y, label=row.names(coords), color=contexts)) +
        geom_point(alpha=0.7) +
        geom_text_repel(max.overlaps = Inf) +
        labs(x="Dim 1", y="Dim 2")
   } else {
      p <- ggplot(coords, ggplot2::aes(x=x, y=y, color=contexts)) +
        geom_point(alpha=0.7) +
        labs(x="Dim 1", y="Dim 2")}
  print(p)
  dev.off()
}


plot_replicates <- function(logical_matrix, contexts, wd, clust_algo, label, n_neigh="default", min_dist=0.01) {
  switch(tolower(clust_algo),
         mca = plot_MCA_replicates(logical_matrix, contexts, wd, label),
         umap = plot_UMAP_replicates(logical_matrix, contexts, wd, label, n_neigh, min_dist)
  )
}


make_batch_logical_matrix <- function(logical_matrix, batches, ratio) {
  logical_matrix <- t(logical_matrix)
  log_mat_batch <- data.frame(genes=row.names(logical_matrix))
  for ( batch in batches ) {
    batch_log_mat <- logical_matrix[,grep(paste0("^",batch),colnames(logical_matrix))]
    if ( nrow(batch_log_mat) < 2 ) {
        log_mat_batch <- cbind(log_mat_batch, batch_log_mat)
    } else {
        log_mat_batch <- cbind(log_mat_batch, (rowSums(batch_log_mat)/ncol(batch_log_mat))>ratio)
    }
  }
  log_mat_batch["genes"] <- NULL
  colnames(log_mat_batch) <- batches
  log_mat_batch <- t(log_mat_batch)

  return(log_mat_batch)
}


plot_MCA_batches <- function(log_mat_batch, batches, wd, label) {
  contexts <- c()
  contexts <- lapply(batches, function(r) {
    c(contexts, unlist(strsplit(r, "_S"))[1])
  }) %>% unlist(.)

  mca_results <- MCA(log_mat_batch, graph=F)
  d <- as.data.frame(mca_results[["ind"]][["coord"]][,1:2])
  d["contexts"] <- as.data.frame(contexts)
  colnames(d) <- c("x", "y", "contexts")
  fig_path <- file.path(wd, "figures")
  if ( !file.exists(fig_path) ) { dir.create(fig_path) }
  plotname <- file.path(fig_path, "mca_plot_batches.pdf")
  pdf(plotname)
  p <- ggplot(d, ggplot2::aes(x=x, y=y, label=row.names(d), color=contexts)) +
    geom_point(alpha=0.7) +
    geom_text_repel(max.overlaps = Inf) +
    labs(x="Dim 1", y="Dim 2")
  print(p)
  dev.off()

}


plot_UMAP_batches <- function(log_mat_batch, batches, wd, label, n_neigh, min_dist) {

  n_neigh <- ifelse(n_neigh=="default", as.integer(length(batches)), n_neigh)
  if ( n_neigh<2) {
    print("Cannot cluster batches if n nearest neighbors is < 1!")
    stop()
  }

  contexts <- c()
  contexts <- lapply(batches, function(r) {
    c(contexts, unlist(strsplit(r, "_S"))[1])
  }) %>% unlist(.)

  fac_matrix <- do.call(cbind, lapply(1:ncol(log_mat_batch), function(n) {
    as.numeric(log_mat_batch[,n])
    }))
  coords <- data.frame(umap(fac_matrix, n_neighbors=n_neigh, metric="euclidean",
                            min_dist=min_dist))
  row.names(coords) <- row.names(log_mat_batch)
  colnames(coords) <- c("x", "y")
  fig_path <- file.path(wd, "figures")
  if ( !file.exists(fig_path) ) { dir.create(fig_path) }
  plotname <- file.path(fig_path, "umap_plot_batches.pdf")
  pdf(plotname)
  p <- ggplot(coords, ggplot2::aes(x=x, y=y, label=row.names(coords), color=contexts)) +
    geom_point(alpha=0.7) +
    geom_text_repel(max.overlaps = Inf) +
    labs(x="Dim 1", y="Dim 2")
  print(p)
  dev.off()

}


plot_batches <- function(log_mat_batch, batches, wd, clust_algo, label, n_neigh="default", min_dist=0.01) {
  switch(tolower(clust_algo),
         mca = plot_MCA_batches(log_mat_batch, batches, wd, label),
         umap = plot_UMAP_batches(log_mat_batch, batches, wd, label, n_neigh, min_dist)
  )
}


make_context_logical_matrix <- function(log_mat_batch, contexts, ratio) {
  contexts <- unique(contexts)
  log_mat_batch <- t(log_mat_batch)
  log_mat_context <- data.frame(genes=row.names(log_mat_batch))
  for ( context in contexts ) {
    context_log_mat <- as.data.frame(log_mat_batch[,grep(paste0("^",context),colnames(log_mat_batch))])
    if ( nrow(context_log_mat) < 2 ) {
        log_mat_context <- cbind(log_mat_context, context_log_mat)
    } else {
        log_mat_context <- cbind(log_mat_context, (rowSums(context_log_mat)/ncol(context_log_mat))>ratio)
    }
  }
  log_mat_context["genes"] <- NULL
  colnames(log_mat_context) <- contexts
  log_mat_context <- t(log_mat_context)

  return(log_mat_context)
}


plot_MCA_contexts <- function(log_mat_context, contexts, wd) {
  contexts <- unique(contexts)
  mca_results <- MCA(log_mat_context, graph=F)
  d <- as.data.frame(mca_results[["ind"]][["coord"]][,1:2])
  d["contexts"] <- as.data.frame(contexts)
  colnames(d) <- c("x", "y", "contexts")
  fig_path <- file.path(wd, "figures")
  if ( !file.exists(fig_path) ) { dir.create(fig_path) }
  plotname <- file.path(fig_path, "mca_plot_contexts.pdf")
  pdf(plotname)
  p <- ggplot(d, ggplot2::aes(x=x, y=y, label=row.names(d), color=contexts)) +
    geom_point(alpha=0.7) +
    geom_text_repel(max.overlaps = Inf) +
    labs(x="Dim 1", y="Dim 2")
  print(p)
  dev.off()

}


plot_UMAP_contexts <- function(log_mat_context, contexts, wd, label, n_neigh, min_dist) {
  contexts <- unique(contexts)
  n_neigh <- ifelse(n_neigh=="default", as.integer(length(contexts)), n_neigh)
  if ( n_neigh<2) {
    print("Cannot cluster contexts if n nearest neighbors is < 1!")
    stop()
  }
  fac_matrix <- do.call(cbind, lapply(1:ncol(log_mat_context), function(n) {
    as.numeric(log_mat_context[,n])
    }))
  coords <- data.frame(umap(fac_matrix, n_neighbors=n_neigh, metric="euclidean",
                            min_dist=min_dist))
  row.names(coords) <- row.names(log_mat_context)
  colnames(coords) <- c("x", "y")
  fig_path <- file.path(wd, "figures")
  if ( !file.exists(fig_path) ) { dir.create(fig_path) }
  plotname <- file.path(fig_path, "umap_plot_contexts.pdf")
  pdf(plotname)
  p <- ggplot(coords, ggplot2::aes(x=x, y=y, label=row.names(coords), color=contexts)) +
    geom_point(alpha=0.7) +
    geom_text_repel(max.overlaps = Inf) +
    labs(x="Dim 1", y="Dim 2")
  print(p)
  dev.off()


}


plot_contexts <- function(log_mat_context, wd, clust_algo, contexts, label,
                          n_neigh="default", min_dist=0.01) {
  switch(tolower(clust_algo),
         mca = plot_MCA_contexts(log_mat_context, contexts, wd, label),
         umap = plot_UMAP_contexts(log_mat_context, contexts, wd, label, n_neigh, min_dist)
  )
}


cluster_samples_main <- function(wd, context_names, technique, clust_algo, label, min_dist=0.01,
                                n_neigh_rep="default", n_neigh_batch="default", n_neigh_cont="default",
                                rep_ratio=0.5, batch_ratio=0.5, quantile=25, min_count="default", seed=12345) {
    set.seed(seed)

    print("Making logical matrix")
    logical_matrix <- make_logical_matrix(wd, technique, context_names)

    res <- parse_contexts(logical_matrix)
    contexts <- res[[1]]
    batches <- res[[2]]
    print("Clustering replicate level filtered data")
    plot_replicates(logical_matrix, contexts, wd, clust_algo, label, n_neigh=n_neigh_rep, min_dist=min_dist)

    print("Clustering batch level filtered data")
    log_mat_batch <- make_batch_logical_matrix(logical_matrix, batches, rep_ratio)
    rm(logical_matrix)
    plot_batches(log_mat_batch, batches, wd, clust_algo, label, n_neigh=n_neigh_batch, min_dist=min_dist)

    print("Clustering context level filtered data")
    log_mat_context <- make_context_logical_matrix(log_mat_batch, contexts, batch_ratio)
    rm(log_mat_batch)
    plot_contexts(log_mat_context, wd, clust_algo, contexts, label, n_neigh=n_neigh_cont, min_dist=min_dist)
}
