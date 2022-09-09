username <- Sys.info()["user"]
work_dir <- str_interp("/home/${username}/work")

if (!dir.exists(str_interp("${work_dir}/py/rlogs"))) {
    dir.create(str_interp("${work_dir}/py/rlogs"))
}

# prevent messy messages from repeatedly writing to juypter
zz <- file(file.path("/home", username, "main", "py", "rlogs", "protein_transform.Rout"), open="wt")
sink(zz, type="message")


library(ggplot2)
library(tidyverse)
library(zoo)


z_result <- function(z_vector, density, mu, stdev, max_y) {
    z_res <- list(
        z = z_vector,
        d = density,
        m = mu,
        s = stdev,
        max_y = max_y
    )

    class(z_vector) <- append(class(z_res), "zFPKM")
    return(z_res)
}


z_score_calc <- function(abundance, min_thresh) {
    if (!is.numeric(abundance)) {
        stop("argument 'abundance' must be numeric")
    }
    abundance <- abundance + 1
    log_abundance_filt <- log(abundance[abundance>min_thresh], base=2)
    log_abundance <- log(abundance, base=2)
    d <- density(log_abundance_filt)

    # calculate rolling average
    perc <- as.integer(0.01*length(d[["y"]]) + 1) # 10% roll avg interval
    d[["roll_y"]] <- zoo::rollmean(d[["y"]], perc)

    # from https://stats.stackexchange.com/questions/22974/how-to-find-local-peaks-valleys-in-a-series-of-data
    find_maxima <- function (x, m = 1){
        shape <- diff(sign(diff(as.numeric(x), na.pad = FALSE)))
        pks <- sapply(which(shape < 0), FUN = function(i){
          z <- i - m + 1
          z <- ifelse(z > 0, z, 1)
          w <- i + m + 1
          w <- ifelse(w < length(x), w, length(x))
          if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
        })
        pks <- unlist(pks)
        pks
    }

    local_maxes <- find_maxima(d[["roll_y"]])
    fit_max <- max(local_maxes) + as.integer(perc/2)

    #mu <- d[["x"]][which.max(d[["y"]])]
    #max_y <- max(d[["y"]])

     # Set the maximum point in the density as the mean for the fitted Gaussian
     mu <- d[["x"]][fit_max] # get max with respect to x) local maxima of rolling
     max_y <- d[["y"]][fit_max]
     cnt <- 0

     while  ( (max_y < 0.5*max(d[["y"]])) && (cnt < 5) ) { # while selected local max y is less than 50% of actual maximum
        cnt <- cnt + 1
        perc <- as.integer((0.05-(cnt*0.01))*length(d[["y"]]) + 1) # rm 1 percent from roll avg interval per iteration

        #d[["roll_y"]] <- filter(data.frame(d[["y"]]), f_2perc, sides=2)
        d[["roll_y"]] <- zoo::rollmean(d[["y"]], perc)
        local_maxes[local_maxes < max(local_maxes)]
        fit_max <- max(local_maxes) + as.integer(perc/2)
        # Set the maximum point in the density as the mean for the fitted Gaussian
        #mu <- d[["x"]][which.max(d[["y"]])]
        mu <- d[["x"]][fit_max] # get max with respect to x) local maxima of rolling
        max_y <- d[["y"]][fit_max]
     }

     if ( (max_y < 0.1*max(d[["y"]])) ) {
        mu <- d[["x"]][which.max(d[["y"]])]
        max_y <- max(d[["y"]]) # if doesnt work use regular zFPKM calculation
     }



    # standard deviation from right side
    #U <- mean(log_abundance[log_abundance > mu])
    #stdev <- (U - mu) * sqrt(pi / 2)
    # standard deviation from left side
    U <- mean(log_abundance[log_abundance < mu & log_abundance > 0])
    stdev <- (mu - U) * sqrt(pi / 2)

    # Compute ztransform
    z <- (log_abundance - mu) / stdev

    result <- z_result(z, d, mu, stdev, max_y)

    return(result)
}


plot_gaussian_fit <- function(results, FacetTitles=TRUE, PlotXfloor) {

    df <- data.frame()

    for (name in names(results)) {
        result <- results[[name]]
        d <- result[["d"]]
        mu <- result[["m"]]
        stdev <- result[["s"]]
        max_y <- result[["max_y"]]

        fitted <- dnorm(d[["x"]], mean=mu, sd=stdev)
        max_abundance <- max_y
        max_fit <- max(fitted)

        scale_fit <- fitted * (max_abundance / max_fit)

        new_df <- data.frame(sample_name=name, log_abundance=d[["x"]], abundance_density=d[["y"]],
                             fitted_density_scaled=scale_fit)

        df <- df %>% dplyr::bind_rows(new_df)
    }

    df_stack <- df %>% tidyr::gather(source, density, -c(log_abundance, sample_name))
    labels <- unique(df_stack$sample_name)

    maximum_x = max(df_stack[["log_abundance"]])
    #maximum_y = max(d[["y"]])

    p <- ggplot2::ggplot(df_stack, ggplot2::aes(x=log_abundance, y=density, color=source)) +
        #ggplot2::facet_wrap(~ sample_name) +
        ggplot2::facet_wrap(vars(sample_name)) +
        ggplot2::geom_line(alpha=0.7) +
        ggplot2::theme_bw() +
        ggplot2::labs(x="log2(abundance)", y="[scaled] density")  +
        ggplot2::theme(legend.position="top") +
        ggplot2::xlim(PlotXfloor, maximum_x)

    print(p)
}


z_transform <- function(abundance_df, min_thresh) {

    abundance_df <- rm_infinite(abundance_df)
    z_df <- data.frame(row.names=row.names(abundance_df))
    outputs <- list()
    for (c in colnames(abundance_df)) {
        output <- z_score_calc(abundance_df[, c], min_thresh)
        z_df[, c] <- output[["z"]]
        outputs[[c]] <- output
    }

    return(list(outputs, z_df))
}


z_score_plot <- function(abundance_df, min_thresh, FacetTitles=FALSE, PlotXfloor=-20) {
    plot_gaussian_fit(z_transform(abundance_df, min_thresh)[[1]], FacetTitles, PlotXfloor)
}


rm_infinite <- function(abundance) {
    # Remove FPKM rows containing all NaN values. These are most likely a result
    # of effective lengths = 0 when calculating FPKM.
    abundance <- as.data.frame(abundance)
    return(abundance[which(!apply(abundance, 1, function(r) all(is.nan(r) | is.infinite(r)))), ])
}


protein_transform_main <- function(abundance_matrix, out_dir, group_name) {
    dir.create(file.path(out_dir, "figures"), showWarnings = FALSE)
    prot <- as.data.frame(read.csv(abundance_matrix))
    prot[is.na(prot)] <- 0
    min_thresh <- min(prot>0)
    #row.names(prot) <- prot$ENTREZ_GENE_ID
    #prot["ENTREZ_GENE_ID"] <- NULL
    pdf(file.path(out_dir, "figures", paste0("fit_proteinAbundance_", group_name, ".pdf")))
    z_score_plot(prot[,-1], 0, out_dir)
    dev.off()
    #if ( i ==1 ) {stop()}
    minimums <- prot == 0
    nas <- is.na(prot)== 1
    z_transformed_abundances <- cbind(prot[,1], z_transform(prot[,c(-1)], 0)[[2]])
    z_transformed_abundances[minimums] <- -4
    #z_transformed_abundances[nas] <- -4
    colnames(z_transformed_abundances)[1] <- "ENTREZ_GENE_ID"
    out_file <- file.path(out_dir, paste0("protein_zscore_Matrix_", group_name, ".csv"))
    write.csv(z_transformed_abundances, out_file, row.names=FALSE)
}

