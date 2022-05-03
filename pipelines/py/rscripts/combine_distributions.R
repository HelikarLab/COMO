# prevent messy messages from repeatedly writing to juypter
zz <- file(file.path("/home", "jupyteruser", "work", "py", "rlogs", "merge_distributions.Rout"), open="wt")
sink(zz, type="message")

library("tidyverse")
library(ggplot2)


parse_contexts_zfpkm <- function(wd, contexts) {

    get_batch_name <- function(x) {
        basename(x)
        return(substring(basename(x), 1, nchar(basename(x))-4))
    }

    batches = list()
    for ( context in contexts ) {
        files <- Sys.glob(file.path(wd, context, "zFPKM_Matrix_*.csv"))
        batches[[context]] <- unlist(lapply(files, get_batch_name))
    }

    return(batches)
}


parse_contexts_zscore_prot <- function(wd, contexts) {

    get_batch_name <- function(x) {
        basename(x)
        return(substring(basename(x), 1, nchar(basename(x))-4))
    }

    batches = list()
    for ( context in contexts ) {
        files <- Sys.glob(file.path(wd, context, "zscore_prot_Matrix_*.csv"))
        batches[[context]] <- unlist(lapply(files, get_batch_name))
    }

    return(batches)
}




merge_batch <- function(wd, context, batch, zmat) {
    print(paste0("Merging ", batch))
    files <- Sys.glob(file.path(wd, context, paste0("*", batch, "*")))
    print(files)
    cnt <- 0
    nrep <- c()
    stopifnot(length(files)>0)

    print("file loop")
    for ( f in files ) {
        zmat <- read.table(f, strip.white=T, header=T, sep=",", row.names=NULL) %>% # read expression matrix
          mutate(across(colnames(.)[-1], as.numeric)) %>% # ensure expression values are numbers
          mutate(across(colnames(.)[1], as.character)) %>% # ensure entrez IDs are character
          group_by(ent) %>%
          summarise_each(funs(max)) %>% # if multiple of same entrez, take max value
          na.omit(.) %>%
          as.data.frame(.)

        nrep <- c(nrep, ncol(zmat)-1)
        ent <- zmat[,"ent"]
        rep_names <- colnames(zmat)
        zmat <- do.call(cbind, lapply(2:ncol(zmat), function(j) {
          repz <- zmat[,j]
        })) %>% cbind(as.character(ent), .) %>% as.data.frame(.) %>% na.omit(.)

        colnames(zmat) <- rep_names

        stack_df <- do.call(rbind, lapply(2:ncol(zmat), function(j) {
          repz <- as.numeric(as.character(zmat[,j]))
          cbind(ent=zmat[,"ent"], zscore=repz, sample=rep(colnames(zmat)[j],length(repz)))
        })) %>% as.data.frame(.)
        stack_df$zscore <- as.numeric(as.character(stack_df$zscore))

        plot_name = file.path(wd, context, "figures", paste0(
          "plot_", context, "_", substring(basename(f), 1, nchar(basename(f))-4), ".pdf"))
        label <- colnames(stack_df)[-1]
        pdf(plot_name)
        p <- ggplot(stack_df, aes(zscore, color=sample)) +
            geom_density()
        print(p)
        dev.off()
    }

    return(list(zmat, nrep))
}


combine_batch_zdistro <- function(wd, context, batch, zmat) {
    print(paste0("Combining ", batch))
    plot_name = file.path(wd, context, "figures", paste0(
        "plot_", context, "_", batch, "_combine_distro", ".pdf"))

    weighted_z <- function(x) {
        x <- as.numeric(x)
        weights = rep((length(x))/length(x), length(x))
        #weights = rep(1, length(x))
        numer = sum(weights*x)
        denom = sqrt(length(x))
        numer/denom
    }

    if ( ncol(zmat)>2 ) {
        combine_z <- apply(zmat[,-1], 1, weighted_z)
        merge_df <- cbind(zmat, combined=combine_z)
        combine_z <- cbind(ent=as.character(zmat[,"ent"]), combine_z)

        stack_df <- do.call(rbind, lapply(2:ncol(merge_df), function(j) {
            repz <- as.numeric(as.character(merge_df[,j]))
            cbind(ent=merge_df[,"ent"], zscore=repz, sample=rep(colnames(merge_df)[j],length(repz)))
        })) %>% as.data.frame(.)
        stack_df$zscore <- as.numeric(as.character(stack_df$zscore))

        label <- colnames(stack_df)[-1]
        pdf(plot_name)
        p <- ggplot(stack_df, aes(zscore, color=sample)) +
            geom_density()
        print(p)
        dev.off()

    } else {
        combine_z <- zmat
    }

    return(as.data.frame(combine_z))
}


combine_context_zdistro <- function(wd, context, n_reps, zmat) {
    print(paste0("Combining ", context, "Z-score distributions. "))
    plot_name = file.path(wd, context, "figures", paste0(
      "plot_", context, "_combine_batches_distro", ".pdf"))

    weighted_z <- function(x) {
        x <- as.numeric(x)
        weights <- c()
        for ( i in 1:length(n_reps) ) {
            weights <- c(weights, (n_reps[i])/sum(n_reps))
        }
        numer = sum(weights*x)
        denom = sqrt(sum(weights^2))
        numer/denom
    }

    if ( ncol(zmat)>2 ) {
        combine_z <- apply(zmat[,-1], 1, weighted_z)
        merge_df <- cbind(zmat, combined=combine_z)
        combine_z <- cbind(ENTREZ_GENE_ID=as.character(zmat[,"ent"]), combine_z)

        stack_df <- do.call(rbind, lapply(2:ncol(merge_df), function(j) {
            repz <- as.numeric(as.character(merge_df[,j]))
            cbind(ent=merge_df[,1], zscore=repz, sample=rep(colnames(merge_df)[j],length(repz)))
        })) %>% as.data.frame(.)
        stack_df$zscore <- as.numeric(as.character(stack_df$zscore))

        label <- colnames(stack_df)[-1]
        pdf(plot_name)
        p <- ggplot(stack_df, aes(zscore, color=sample)) +
            geom_density()
        print(p)
        dev.off()

    } else {
        combine_z <- zmat
        colnames(combine_z) <- c("ENTREZ_GENE_ID", "combine_z")
    }

    return(combine_z)
}


combine_zscores_main <- function(wd, contexts, use_rna, use_proteins) {
    # TODO: proteomics combine!
    fig_path <- file.path(wd, "figures")
    if ( !file.exists(fig_path) ) { dir.create(fig_path) }

    if ( use_rna ) {
        batches_rna <- parse_contexts_zfpkm(wd, contexts)
        batches_prot <- parse_contexts_zscore_prot(wd, contexts)
        for ( context in contexts ) {
            cont_batches <- batches_rna[[context]]
            nreps <- c()
            cnt <- 0
            for ( batch in cont_batches ) {
                print(paste0("combining ", batch, " for context ", context, "'s RNA-seq"))
                res <- merge_batch(wd, context, batch)
                zmat <- res[[1]]
                nreps <- c(nreps, res[[2]])
                comb_z <- combine_batch_zdistro(wd, context, batch, zmat)
                colnames(comb_z) <- c("ent", batch)
                if (!cnt) { merge_z <- comb_z }
                else { merge_z <- merge_z %>% inner_join(comb_z, by="ent") }
                cnt <- cnt+1
            }

            comb_batches_z_rna <- combine_context_zdistro(wd, context, nreps, merge_z)
            filename <- file.path(wd, context, paste0("combined_zFPKM_", context, ".csv"))
            write.csv(comb_batches_z, filename, row.names=FALSE)
            if ( !use_proteins ) {
                filename <- file.path(wd, context, paste0("model_scores_", context, ".csv"))
                write.csv(comb_batches_z_rna, filename, row.names=FALSE)
            }
        }
    }

    if ( use_proteins ) {
        for ( context in contexts ) {
            cont_batches <- batches_prot[[context]]
            nreps <- c()
            cnt <- 0
            for ( batch in cont_batches ) {
                print(paste0("combining ", batch, " for context ", context, "'s protein abundance"))
                res <- merge_batch(wd, context, batch)
                zmat <- res[[1]]
                nreps <- c(nreps, res[[2]])
                comb_z <- combine_batch_zdistro(wd, context, batch, zmat)
                colnames(comb_z) <- c("ent", batch)
                if (!cnt) { merge_z <- comb_z }
                else { merge_z <- merge_z %>% inner_join(comb_z, by="ent") }
                cnt <- cnt+1
            }

            comb_batches_z_prot <- combine_context_zdistro(wd, context, nreps, merge_z)
            filename <- file.path(wd, context, paste0("combined_zscore_proteinAbundance_", context, ".csv"))
            write.csv(comb_batches_z_prot, filename, row.names=FALSE)
            if ( !use_proteins ) {
                filename <- file.path(wd, context, paste0("model_scores_", context, ".csv"))
                write.csv(comb_batches_z_prot, filename, row.names=FALSE)
            }
        }
    }

    if ( use_rna & use_proteins ) {
        filename <- file.path(wd, context, paste0("model_scores_", context, ".csv"))
        write.csv(comb_batches_z, filename, row.names=FALSE)

    }
}