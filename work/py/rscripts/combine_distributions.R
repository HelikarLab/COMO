username <- Sys.info()["user"]
work_dir <- str_interp("/home/${username}/work")

if (!dir.exists(str_interp("${work_dir}/py/rlogs"))) {
    dir.create(str_interp("${work_dir}/py/rlogs"))
}

# prevENTREZ_GENE_ID messy messages from repeatedly writing to juypter
zz <- file(file.path("/home", username, "work", "py", "rlogs", "combine_distributions.Rout"), open="wt")
sink(zz, type="message")


library(tidyverse)
library(ggplot2)


parse_contexts_zfpkm <- function(wd, contexts, prep) {

    get_batch_name <- function(x) {
        basename(x)
       return(substring(basename(x), 1, nchar(basename(x))-4))
    }

    batches = list()
    for ( context in contexts ) {
        files <- Sys.glob(file.path(wd, context, prep, paste0("zFPKM_Matrix_", prep, "_*.csv")))
        batches[[context]] <- unlist(lapply(files, get_batch_name))
    }

    return(batches)
}


parse_contexts_zumi <- function(wd, contexts, prep) {

    get_batch_name <- function(x) {
        basename(x)
       return(substring(basename(x), 1, nchar(basename(x))-4))
    }

    batches = list()
    for ( context in contexts ) {
        files <- Sys.glob(file.path(wd, context, prep, paste0("zUMI_Matrix_", prep, "_*.csv")))
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
    files <- Sys.glob(file.path(wd, "*", "proteomics", "protein_zscore_Matrix_*.csv"))
    for ( context in contexts ) {
        files <- Sys.glob(file.path(wd, context, "proteomics", "protein_zscore_Matrix_*.csv"))
        batches[[context]] <- unlist(lapply(files, get_batch_name))
    }

    return(batches)
}


merge_batch <- function(wd, context, batch) {
    print(paste0("Merging ", batch))
    files <- Sys.glob(file.path(wd, paste0("*", batch, "*")))
    cnt <- 0
    nrep <- c()
    stopifnot(length(files)>0)

    for ( f in files ) {
        zmat <- read.table(f, strip.white=T, header=T, sep=",", row.names=NULL) %>% # read expression matrix
            mutate(across(colnames(.)[-1], as.numeric)) %>% # ensure expression values are numbers
            mutate(across(colnames(.)[1], as.character)) %>% # ensure ENTREZ_GENE_IDrez IDs are character
            group_by(ENTREZ_GENE_ID) %>%
            summarise_each(funs(max)) %>% # if multiple of same ENTREZ_GENE_ID, take max value
            na.omit(.) %>%
            as.data.frame(.)

        nrep <- c(nrep, ncol(zmat)-1)
        entrez_gene <- zmat[,"ENTREZ_GENE_ID"]
        rep_names <- colnames(zmat)
        zmat <- do.call(cbind, lapply(2:ncol(zmat), function(j) {
            repz <- zmat[,j]
        })) %>%
            cbind(as.character(entrez_gene), .) %>%
            as.data.frame(.) %>%
            na.omit(.)

        colnames(zmat) <- rep_names

        stack_df <- do.call(rbind, lapply(2:ncol(zmat), function(j) {
            repz <- as.numeric(as.character(zmat[,j]))
            cbind(ENTREZ_GENE_ID=zmat[,"ENTREZ_GENE_ID"], zscore=repz, source=rep(colnames(zmat)[j],length(repz)))
        })) %>% as.data.frame(.)

        stack_df$zscore <- as.numeric(as.character(stack_df$zscore))

        plot_name_pdf = file.path(wd, "figures", paste0(
            "plot_", context, "_", substring(basename(f), 1, nchar(basename(f))-4), ".pdf"))
		
		plot_name_png = file.path(wd, "figures", paste0(
            "plot_", context, "_", substring(basename(f), 1, nchar(basename(f))-4), ".png"))
			
        simplified_plot <- ifelse(length(unique(stack_df$source)) > 10, TRUE, FALSE)
        label <- colnames(stack_df)[-1]
        pdf(plot_name_pdf)
		png(
			plot_name_png,
			res=1200,
			units="in",
			width=3.25,
			height=3.25
		)
        p <- ggplot(stack_df, aes(zscore, color=source)) +
            geom_density() +
			theme(text=element_text(size=12,  family="sans"))
			
        if ( simplified_plot ) {
            p <- p + theme(legend.position = "none")
        }
        max_dens <- 0
        # get y upper limit by finding density of peak at z = 0
        for ( source in unique(p$data$source) ) {
            source_z <- p$data$zscore[p$data$source == source] %>% .[!is.na(.)] #%>% .[!is.nan(.)]
            source_densx <- density(source_z)$x
            source_densy <- density(source_z)$y
            idx <- min(max(which(source_densx <= 0)), min(which(source_densx == 0)))
            max_d <- source_densy[idx]
            if ( max_d > max_dens ) { max_dens <- max_d }
        }
        #p <- p + ylim(0, 1.5*max_dens)
        print(p)
        dev.off()
    }

    return(list(zmat, nrep))
}


combine_batch_zdistro <- function(wd, context, batch, zmat) {
    print(paste0("Combining ", batch))
    plot_name_pdf = file.path(wd, "figures", paste0("plot_", context, "_", batch, "_combine_distro", ".pdf"))
	plot_name_png = file.path(wd, "figures", paste0("plot_", context, "_", batch, "_combine_distro", ".png"))

    weighted_z <- function(x) {
        floor_score <- -6
        ceil_score <- 6
        x <- as.numeric(x)
        numer = sum(x)
        denom = sqrt(length(x))
        result = numer/denom
        if ( result < floor_score ) { result <- floor_score }
        if ( result > ceil_score ) { result <- ceil_score }
        return(result)
    }

    if ( ncol(zmat)>2 ) {
        combine_z <- apply(zmat[,-1], 1, weighted_z)
        merge_df <- cbind(zmat, combined=combine_z)
        combine_z <- cbind(ENTREZ_GENE_ID=as.character(zmat[,"ENTREZ_GENE_ID"]), combine_z)

        stack_df <- do.call(rbind, lapply(2:ncol(merge_df), function(j) {
            repz <- as.numeric(as.character(merge_df[,j]))
            cbind(ENTREZ_GENE_ID=merge_df[,"ENTREZ_GENE_ID"], zscore=repz, source=rep(colnames(merge_df)[j],length(repz)))
        })) %>% as.data.frame(.)
        stack_df$zscore <- as.numeric(as.character(stack_df$zscore))

        simplified_plot <- ifelse(length(unique(stack_df$source)) > 10, TRUE, FALSE)
        label <- colnames(stack_df)[-1]
        pdf(plot_name_pdf)
		png(
			plot_name_png,
			res=1200,
			units="in",
			width=3.25,
			height=3.25
		)
		
        if ( simplified_plot ) {
            #p <- p + theme(legend.position = "none")
            stack_df <- stack_df[stack_df$source=="combined",]
        }

        p <- ggplot(stack_df, aes(zscore, color=source)) +
            geom_density() +
			theme(text=element_text(size=12,  family="sans"))
			
        max_dens <- 0
        # get y upper limit by finding density of peak at z = 0
        for ( source in unique(p$data$source) ) {
            source_z <- p$data$zscore[p$data$source == source] %>% .[!is.na(.)] #%>% .[!is.nan(.)]
            source_densx <- density(source_z)$x
            source_densy <- density(source_z)$y
            idx <- min(max(which(source_densx <= 0.5)), min(which(source_densx == 0.5)))
            max_d <- source_densy[idx]
            if ( max_d > max_dens ) { max_dens <- max_d }
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
    print(paste0("Combining ", context, " Z-score distributions. "))
    plot_name_pdf = file.path(wd, "figures", paste0(
        "plot_", context, "_combine_batches_distro", ".pdf"))
	plot_name_png = file.path(wd, "figures", paste0(
        "plot_", context, "_combine_batches_distro", ".png"))

    weighted_z <- function(x, n_reps) {
        floor_score <- -6
        ceil_score <- 6
        x <- as.numeric(x)
        nas <- sort(   unique(  c( which(is.nan(x)), which(is.na(x)) )  )   )
        weights <- c()
        for ( i in 1:length(n_reps) ) { weights <- c(weights, (n_reps[i])/sum(n_reps)) }
        if ( length(nas) > 0 ) {
            x <- x[-nas]
            weights <- weights[-nas]
        }
        numer = sum(weights*x)
        denom = sqrt(sum(weights^2))
        result <- numer/denom
        if ( result < floor_score ) { result <- floor_score }
        if ( result > ceil_score ) { result <- ceil_score }
        return(result)
    }

    if ( ncol(zmat)>2 ) {
        combine_z <- apply(zmat[,-1], 1, weighted_z, n_reps=n_reps)
        merge_df <- cbind(zmat, combined=combine_z)
        combine_z <- cbind(ENTREZ_GENE_ID=as.character(zmat[,"ENTREZ_GENE_ID"]), combine_z)

        stack_df <- do.call(rbind, lapply(2:ncol(merge_df), function(j) {
            repz <- as.numeric(as.character(merge_df[,j]))
            cbind(ENTREZ_GENE_ID=merge_df[,1], zscore=repz, source=rep(colnames(merge_df)[j],length(repz)))
        })) %>% as.data.frame(.)
        stack_df$zscore <- as.numeric(as.character(stack_df$zscore))

        label <- colnames(stack_df)[-1]
        pdf(plot_name_pdf)
		png(
			plot_name_png,
			res=1200,
			units="in",
			width=3.25,
			height=3.25
		)
        p <- ggplot(stack_df, aes(zscore, color=source)) +
            geom_density() +
			theme(text=element_text(size=12,  family="sans"))

        max_dens <- 0
        # get y upper limit by finding density of peak at z = 0
        for ( source in unique(p$data$source) ) {
            source_z <- p$data$zscore[p$data$source == source] %>% .[!is.na(.)] #%>% .[!is.nan(.)]
            source_densx <- density(source_z)$x
            source_densy <- density(source_z)$y
            idx <- min(max(which(source_densx <= 0)), min(which(source_densx == 0)))
            max_d <- source_densy[idx]
            if ( max_d > max_dens ) { max_dens <- max_d }
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
    keep_gene_scores=TRUE) {

    print(paste0("Combining -omics distributions for ", context))

    fig_path <- file.path(wd, context, "figures")
    if ( !file.exists(fig_path) ) { dir.create(fig_path) }
    plot_name_pdf = file.path(fig_path, paste0("plot_", context, "_combine_omics_distro", ".pdf"))
	plot_name_png = file.path(fig_path, paste0("plot_", context, "_combine_omics_distro", ".png"))
	
    weights <- c()
    names <- c()
    dfs <- list()
    counter <- 0
    if ( tweight > 0 ) {
        counter <- counter+1
        weights <- c(weights, tweight)
        names <- c(names, "total")
        dfs[[counter]] <- comb_batches_z_trna
    }
    if ( mweight > 0 ) {
        counter <- counter+1
        weights <- c(weights, mweight)
        names <- c(names, "polyA")
        dfs[[counter]] <- comb_batches_z_mrna
    }
    if ( sweight > 0 ) {
        counter <- counter+1
        weights <- c(weights, sweight)
        names <- c(names, "singleCell")
        dfs[[counter]] <- comb_batches_z_scrna
    }
    if ( pweight > 0 ) {
        counter <- counter+1
        weights <- c(weights, pweight)
        names <- c(names, "proteome")
        dfs[[counter]] <- comb_batches_z_prot
    }

    weighted_z <- function(x, weights) {
        floor_score <- -6
        ceil_score <- 10
        x <- as.numeric(x)

        nas <- which(is.na(x))
        if ( length(nas) > 0 ) {
            x <- x[-nas]
            weights <- weights[-nas]
        }
        weights <- weights/sum(weights)
        numer = sum(weights*x)
        denom = sqrt(sum(weights^2))
        result <- numer/denom
        if ( result < floor_score ) { result <- floor_score }
        if ( result > ceil_score ) { result <- ceil_score }
        return(result)
    }

    for ( i in 1:counter ) {
        add_df <- dfs[[i]]
        colnames(add_df)[2] <- names[i]
        if ( i == 1 ) { zmat <- add_df }
        else { zmat <- full_join(zmat, add_df, by="ENTREZ_GENE_ID", copy=TRUE) }
    }

    if ( ncol(zmat)>2 ) {
        combine_z <- apply(zmat[,-1], 1, weighted_z, weights=weights)
    } else {
        combine_z = zmat[,-1]
    }
	
    merge_df <- cbind(zmat, combined=combine_z)
    combine_z <- cbind(ENTREZ_GENE_ID=as.character(zmat[,"ENTREZ_GENE_ID"]), combine_z)

    stack_df <- do.call(rbind, lapply(2:ncol(merge_df), function(j) {
        repz <- as.numeric(as.character(merge_df[,j]))
        cbind(ENTREZ_GENE_ID=merge_df[,1], zscore=repz, source=rep(colnames(merge_df)[j],length(repz)))
    })) %>% as.data.frame(.)
    stack_df$zscore <- as.numeric(as.character(stack_df$zscore))

    label <- colnames(stack_df)[-1]
    pdf(plot_name_pdf)
	png(
		plot_name_png,
		res=1200,
		units="in",
		width=3.25,
		height=3.25
	)

    p <- ggplot(stack_df, aes(zscore, color=source)) +
        geom_density() +
		theme(text=element_text(size=12,  family="sans"))
		
    max_dens <- 0
    # get y upper limit by finding density of peak at z = 0
    for ( source in unique(p$data$source) ) {
        source_z <- p$data$zscore[p$data$source == source] %>% .[!is.na(.)] #%>% .[!is.nan(.)]
        source_densx <- density(source_z)$x
        source_densy <- density(source_z)$y
        idx <- min(max(which(source_densx <= 0)), min(which(source_densx == 0)))
        max_d <- source_densy[idx]
        if ( max_d > max_dens ) { max_dens <- max_d }
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
    wd,
    contexts,
    use_mrna_flag,
    use_trna_flag,
    use_scrna_flag,
    use_proteins_flag,
    keep_gene_scores,
    tweight_master=1,
    mweight_master=1,
    sweight_master=1,
    pweight_master=2) {

    fig_path <- file.path(wd, "figures")
    if ( !file.exists(fig_path) ) { dir.create(fig_path) }

    tbatches <- parse_contexts_zfpkm(wd, contexts, "total")
    mbatches <- parse_contexts_zfpkm(wd, contexts, "mrna")
    sbatches <- parse_contexts_zumi(wd, contexts, "scrna")
    pbatches <- parse_contexts_zscore_prot(wd, contexts)

    for ( context in contexts ) {
        use_mrna <- use_mrna_flag
        use_trna <- use_trna_flag
        use_scrna <- use_scrna_flag
        use_proteins <- use_proteins_flag

        tweight <- tweight_master
        mweight <- mweight_master
        sweight <- sweight_master
        pweight <- pweight_master

        cont_tbatches <- tbatches[[context]]
        cont_mbatches <- mbatches[[context]]
        cont_sbatches <- sbatches[[context]]
        cont_pbatches <- pbatches[[context]]

        print(paste0("Context: ", context))

        if ( length(cont_tbatches) == 0 & use_trna_flag ) {
            use_trna <- FALSE
            print(paste0("No total RNA-seq zFPKM Matrix files found for ", context, " will not use for this context."))
        }

        if ( length(cont_mbatches) == 0 & use_mrna_flag ) {
            use_mrna <- FALSE
            print(paste0("No polyA RNA-seq zFPKM Matrix files found for ", context, ". Will not use for this context."))
        }

        if ( length(cont_sbatches) == 0 & use_scrna_flag ) {
            use_scrna <- FALSE
            print(paste0("No SC RNA-seq zFPKM Matrix files found for ", context, ". Will not use for this context."))
        }

        if ( length(cont_pbatches) == 0 & use_proteins_flag ) {
            use_proteins <- FALSE
            print(paste0("No proteomics z-score Matrix files found for ", context, ". Will not use for this context."))
        }

        if ( use_trna ) {
            print("Merging total RNA-seq distributions...")
            twd <- file.path(wd, context, "total")
            nreps <- c()
            cnt <- 0
            for ( batch in cont_tbatches ) {
                print(paste0("Batch: ", batch))
                res <- merge_batch(twd, context, batch)
                zmat <- res[[1]]
                nreps <- c(nreps, res[[2]])
                comb_z <- combine_batch_zdistro(twd, context, batch, zmat)
                colnames(comb_z) <- c("ENTREZ_GENE_ID", batch)
                if (!cnt) { merge_z <- comb_z }
                else { merge_z <- full_join(merge_z, comb_z, by="ENTREZ_GENE_ID") }
                cnt <- cnt+1
            }

            comb_batches_z_trna <- combine_context_zdistro(twd, context, nreps, merge_z)
            filename <- file.path(twd, paste0("combined_zFPKM_", context, ".csv"))
            write.csv(comb_batches_z_trna, filename, row.names=FALSE)

            if ( !use_proteins & !use_mrna & !use_scrna) {
                filename <- file.path(wd, context, "total", paste0("model_scores_", context, ".csv"))
                write.csv(comb_batches_z_trna, filename, row.names=FALSE)
             }

        } else { comb_batches_z_trna <- NA }


        if ( use_mrna ) {
            print("Merging polyA enriched RNA-seq distributions...")
            mwd <- file.path(wd, context, "mrna")
            nreps <- c()
            cnt <- 0
            for ( batch in cont_mbatches ) {
                print(paste0("Batch: ", batch))
                res <- merge_batch(mwd, context, batch)
                zmat <- res[[1]]
                nreps <- c(nreps, res[[2]])
                comb_z <- combine_batch_zdistro(mwd, context, batch, zmat)
                colnames(comb_z) <- c("ENTREZ_GENE_ID", batch)
                if (!cnt) { merge_z <- comb_z }
                else { merge_z <- full_join(merge_z, comb_z, by="ENTREZ_GENE_ID") }
                cnt <- cnt+1
            }

            comb_batches_z_mrna <- combine_context_zdistro(mwd, context, nreps, merge_z)
            filename <- file.path(mwd, paste0("combined_zFPKM_", context, ".csv"))
            write.csv(comb_batches_z_mrna, filename, row.names=FALSE)

            if ( !use_proteins & !use_trna & !use_scrna) {
                filename <- file.path(mwd, paste0("model_scores_", context, ".csv"))
                write.csv(comb_batches_z_mrna, filename, row.names=FALSE)
            }

        } else { comb_batches_z_mrna <- NA }


        if ( use_scrna ) {
            print("Merging single-cell RNA-seq distributions...")
            swd <- file.path(wd, context, "scrna")
            nreps <- c()
            cnt <- 0
            for ( batch in cont_sbatches ) {
                print(paste0("Batch: ", batch))
                res <- merge_batch(swd, context, batch)
                zmat <- res[[1]]
                nreps <- c(nreps, res[[2]])
                comb_z <- combine_batch_zdistro(swd, context, batch, zmat)
                colnames(comb_z) <- c("ENTREZ_GENE_ID", batch)
                if (!cnt) { merge_z <- comb_z }
                else { merge_z <- full_join(merge_z, comb_z, by="ENTREZ_GENE_ID") }
                cnt <- cnt+1
            }

            comb_batches_z_scrna <- combine_context_zdistro(swd, context, nreps, merge_z)
            filename <- file.path(swd, paste0("combined_zFPKM_", context, ".csv"))
            write.csv(comb_batches_z_scrna, filename, row.names=FALSE)

            if ( !use_proteins & !use_trna & !use_mrna) {
                filename <- file.path(swd, paste0("model_scores_", context, ".csv"))
                write.csv(comb_batches_z_scrna, filename, row.names=FALSE)
            }

        } else { comb_batches_z_scrna <- NA }


        if ( use_proteins ) {
            print("Merging protein abundance distributions...")
            pwd = file.path(wd, context, "proteomics")
            nreps <- c()
            cnt <- 0
            for ( batch in cont_pbatches ) {
                print(paste0("Batch: ", batch))
                res <- merge_batch(pwd, context, batch)
                zmat <- res[[1]]
                nreps <- c(nreps, res[[2]])
                comb_z <- combine_batch_zdistro(pwd, context, batch, zmat)
                colnames(comb_z) <- c("ENTREZ_GENE_ID", batch)
                if (!cnt) { merge_z <- comb_z }
                else { merge_z <- full_join(merge_z, comb_z, by="ENTREZ_GENE_ID") }
                cnt <- cnt+1
            }

            comb_batches_z_prot <- combine_context_zdistro(pwd, context, nreps, merge_z)
            filename <- file.path(pwd, paste0("combined_zscore_proteinAbundance_", context, ".csv"))
            write.csv(comb_batches_z_prot, filename, row.names=FALSE)

            if ( !use_mrna & !use_trna & !use_scrna ) {
                filename <- file.path(pwd, paste0("model_scores_", context, ".csv"))
                write.csv(comb_batches_z_prot, filename, row.names=FALSE)
            }

        } else { comb_batches_z_prot <- NA }

        if ( !use_trna ) { tweight <- 0 }
        if ( !use_mrna ) { mweight <- 0 }
        if ( !use_scrna) { sweight <- 0 }
        if ( !use_proteins ) { pweight <- 0 }

        comb_omics_z <- combine_omics_zdistros(
            wd,
            context,
            comb_batches_z_trna,
            comb_batches_z_mrna,
            comb_batches_z_scrna,
            comb_batches_z_prot,
            tweight,
            mweight,
            sweight,
            pweight
        )

        filename <- file.path(wd, context, paste0("model_scores_", context, ".csv"))
        write.csv(comb_omics_z, filename, row.names=FALSE)
    }

}