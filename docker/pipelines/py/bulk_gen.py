#!/usr/bin/python3

import os, time, sys
import pandas as pd
import getopt
from rpy2.robjects.packages import importr
from rpy2.robjects import r, pandas2ri
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage
from project import configs
import unidecode

pandas2ri.activate()

limma = importr("limma", lib_loc="/home/jupyteruser/rlibs"
tidyverse = importr("tidyverse", lib_loc="/home/jupyteruser/rlibs")
edgeR = importr("edgeR", lib_loc="/home/jupyteruser/rlibs")
genefilter = importr("genefilter", lib_loc="/home/jupyteruser/rlibs")
biomaRt = importr("biomaRt", lib_loc="/home/jupyteruser/rlibs")
sjmisc = importr("sjmisc", lib_loc="/home/jupyteruser/rlibs")

# automatically convert ryp2 dataframe to Pandas dataframe
string="""
library(tidyverse)
library(limma)
library(edgeR)
library(genefilter)
library(biomaRt)
library(sjmisc)

readCountMatrix <- function(cmat_file, config_file) {
  SampMetrics <- list()
  conf <- read.table(config_file, header=TRUE, sep=",")
  cmat_whole <- read.table(cmat_file, header=TRUE, sep=",")
  genes <- cmat_whole[,1]
  i=0
  insert_sizes <- c()
  samp_names <- c()
  for ( entry in conf$SampleName ) {
    i<-i+1
    if ( i==1 ) {
        #print(paste(c("reading counts for ", entry),collapse=""))
    } else if ( grepl("GROUP", entry) ) {
      if ( i>2 ) {
        rnames <- samp_mat[,"samp_mat"]
        samp_mat <- samp_mat[,-1]
        samp_mat <- sapply(data.frame(samp_mat), as.numeric)
        samp_mat <- as.matrix(samp_mat)
        rownames(samp_mat) <- rnames
        colnames(samp_mat) <- samp_names
        SampMetrics[[group]][["CountMatrix"]] <- samp_mat
        SampMetrics[[group]][["NumSamples"]] <- ncol(samp_mat)-1
        SampMetrics[[group]][["InsertSizes"]] <- insert_sizes
        insert_sizes <- c()
        samp_names <- c()
      }
      group <- unlist(str_split(entry, "_", n=2))[2]
      samp_mat <- genes
    } else if ( entry %in% colnames(cmat_whole) ) {
      samp_mat <- cbind(samp_mat, cmat_whole[,entry])
      insert_sizes <- c(insert_sizes, conf$InsertSize[i])
      samp_names <- c(samp_names, entry)
    } else {
      print(paste(c("config file entry ", entry, " not found in count matrix file."),
                    collapse=""))
    }
  }
  rnames <- samp_mat[,"samp_mat"]
  samp_mat <- samp_mat[,-1]
  samp_mat <- sapply(data.frame(samp_mat), as.numeric)
  samp_mat <- as.matrix(samp_mat)
  rownames(samp_mat) <- rnames
  colnames(samp_mat) <- samp_names

  # remove version numbers from ensembl id
  for ( j in 1:length(rownames(samp_mat)) ) {
    r <- rownames(samp_mat)[j]
    if ( grepl("\\\\.", r) ) {
      gen <- unlist(str_split(r, "\\\\."))[1]
      rownnames(samp_mat)[j] <- gen
    }
  }
  SampMetrics[[group]][["CountMatrix"]] <- samp_mat
  SampMetrics[[group]][["NumSamples"]] <- ncol(samp_mat)-1
  SampMetrics[[group]][["InsertSizes"]] <- insert_sizes

  return(SampMetrics)
}

calculateTPM <- function(SampMetrics, cell_types) {
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

calculateZscore <- function(SampMetrics, cell_types, norm_tech) {
  for ( i in 1:length(SampMetrics)) {
    if ( norm_tech=="CPM" ) {
      tmat <- SampMetrics[[i]][["CPM_Matrix"]]
    } else if ( norm_tech=="TPM" ) {
      tmat <- SampMetrics[[i]][["TPM_Matrix"]]
    }
    zmat <- matrix(nrow=nrow(tmat), ncol=ncol(tmat))
    rownames(zmat) <- rownames(tmat)
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

CPM_filter <- function(SampMetrics, filt_options, cell_types) {
  N <- filt_options$pos_rep
  perc_top <- filt_options$top_percentile
  min.count <- filt_options$min_count
  for ( i in 1:length(SampMetrics) ) {
    counts <- SampMetrics[[i]][["CountMatrix"]]
    ens <- SampMetrics[[i]][["ENSEMBL"]]
    ent <- SampMetrics[[i]][["Entrez"]]
    size <- SampMetrics[[i]][["GeneSizes"]]
    symb <- SampMetrics[[i]][["GeneSymbol"]]
    lib.size <- colSums(counts)
    MedianLibSize <- median(lib.size)
    if ( min.count=="default" ) {
      CPM.Cutoff <- 10000000/(median(colSums(counts)))
    } else {
      CPM.Cutoff <- min.count/MedianLibSize*1e6
    }
    CPM <- cpm(counts,lib.size=lib.size)
    min.samples <- round(N * ncol(counts))
    test_bools_top <- data.frame(genes=ens)
    for ( j in 1:ncol(CPM) ) {
      cpm_q <- CPM[,j]
      cpm_q <- cpm_q[cpm_q>0]
      q_cutoff_top <- quantile(cpm_q, prob=1-perc_top/100)
      test_bools_top <- cbind(test_bools_top, as.integer(CPM[,j]>q_cutoff_top))
    }

    f1 <- genefilter::kOverA(min.samples, CPM.Cutoff)
    flist <- genefilter::filterfun(f1)
    keep <- genefilter::genefilter(CPM, flist)
    SampMetrics[[i]][["ENSEMBL"]] <- ens[keep]
    SampMetrics[[i]][["Entrez"]] <- ent[keep]
    SampMetrics[[i]][["GeneSizes"]] <- size[keep]
    SampMetrics[[i]][["CountMatrix"]] <- counts[keep,]
    SampMetrics[[i]][["CPM_Matrix"]] <- CPM[keep,]
    SampMetrics[[i]][["GeneSymbol"]] <- symb[keep]

    # top percentile genes
    test_bools_top["genes"] <- NULL
    f1_top <- genefilter::kOverA(1, 0.9)
    flist_top <- genefilter::filterfun(f1_top)
    keep_top <- genefilter::genefilter(test_bools_top, flist_top)
    SampMetrics[[i]][["ENSEMBL_top"]] <- ens[keep_top]
    SampMetrics[[i]][["Entrez_top"]] <- ent[keep_top]
    SampMetrics[[i]][["GeneSymbol_top"]] <- symb[keep_top]
  }
  SampMetrics <- calculateZscore(SampMetrics, cell_types, "CPM")
  return(SampMetrics)
}

TPM_quant_filter <- function(SampMetrics, filt_options, cell_types) {
  N <- filt_options$pos_rep
  perc <- filt_options$percentile
  perc_top <- filt_options$top_percentile
  SampMetrics <- calculateTPM(SampMetrics, cell_types)
  for ( i in 1:length(SampMetrics) ) {

    counts <- SampMetrics[[i]][["CountMatrix"]]
    ens <- SampMetrics[[i]][["ENSEMBL"]]
    ent <- SampMetrics[[i]][["Entrez"]]
    size <- SampMetrics[[i]][["GeneSizes"]]
    symb <- SampMetrics[[i]][["GeneSymbol"]]
    tpm <- SampMetrics[[i]][["TPM_Matrix"]]
    min.samples <- round(N * ncol(tpm))
    test_bools <- data.frame(genes=ens)
    test_bools_top <- test_bools
    for ( j in 1:ncol(tpm) ) {
      tpm_q <- tpm[,j]
      tpm_q <- tpm_q[tpm_q>0]
      q_cutoff <- quantile(tpm_q, prob=1-perc/100)
      q_cutoff_top <- quantile(tpm_q, prob=1-perc_top/100)
      test_bools <- cbind(test_bools, as.integer(tpm[,j]>q_cutoff))
      test_bools_top <- cbind(test_bools_top, as.integer(tpm[,j]>q_cutoff_top))
    }
    test_bools["genes"] <- NULL
    test_bools_top["genes"] <- NULL
    f1 <- genefilter::kOverA(min.samples, 0.9)
    flist <- genefilter::filterfun(f1)
    keep <- genefilter::genefilter(test_bools, flist)
    SampMetrics[[i]][["ENSEMBL"]] <- ens[keep]
    SampMetrics[[i]][["Entrez"]] <- ent[keep]
    SampMetrics[[i]][["GeneSizes"]] <- size[keep]
    SampMetrics[[i]][["CountMatrix"]] <- counts[keep,]
    SampMetrics[[i]][["TPM_Matrix"]] <- tpm[keep,]
    SampMetrics[[i]][["GeneSymbol"]] <- symb[keep]
    f1_top <- genefilter::kOverA(min.samples, 0.9)
    flist_top <- genefilter::filterfun(f1_top)
    keep_top <- genefilter::genefilter(test_bools_top, flist_top)
    SampMetrics[[i]][["ENSEMBL_top"]] <- ens[keep_top]
    SampMetrics[[i]][["Entrez_top"]] <- ent[keep_top]
    SampMetrics[[i]][["GeneSymbol_top"]] <- symb[keep_top]
  }
  SampMetrics <- calculateZscore(SampMetrics, cell_types, "TPM")
  return(SampMetrics)
}

zFPKM_filter <- function(SampMetrics, filt_options, mart, cell_types) {
  N <- filt_options$pos_rep
  perc_top <- filt_options$top_percentile
  for ( i in 1:length(SampMetrics) ) {
    cmat <- SampMetrics[[i]][["CountMatrix"]]
    gnames <- SampMetrics[[i]][["ENSEMBL"]]
    rownames(cmat) <- make.names(gnames, unique = TRUE)
    fmat <- fpkm(cmat, SampMetrics[[i]][["GeneSizes"]],
                 SampMetrics[[i]][["InsertSizes"]])

    ent <- SampMetrics[[i]][["Entrez"]]
    size <- SampMetrics[[i]][["GenesSizes"]]
    symb <- SampMetrics[[i]][["GeneSymbol"]]

    fdf <- data.frame(fmat)
    zmat <- zFPKM(fdf, assayName="FPKM")
    inames <- rownames(zmat)

    min.samples <- round(N * ncol(zmat))

    test_bools_top <- data.frame(genes=ens)
    for ( j in 1:ncol(zmat) ) {
      z_q <- zmat[,j]
      z_q <- z_q[z_q>-3]
      q_cutoff_top <- quantile(z_q, prob=1-perc_top/100)
      test_bools_top <- cbind(test_bools_top, as.integer(zmat[,j]>q_cutoff_top))
    }

    cutoff <- -3
    f1 <- genefilter::kOverA(min.samples, cutoff)
    flist <- genefilter::filterfun(f1)
    keep <- genefilter::genefilter(zmat, flist)
    SampMetrics[[i]][["zFPKM_Matrix"]] <- zmat[keep,]
    SampMetrics[[i]][["ENSEMBL"]] <- inames[keep]
    SampMetrics[[i]][["FPKM_Matrix"]] <- fmat

    gene_inf <- getBM(filters="ensembl_gene_id",
                      attributes=c("ensembl_gene_id", "entrezgene_id",
                                   "start_position", "end_position"),
                      values=SampMetrics[[i]][["ENSEMBL"]], mart=mart)

    gene_inf$size <- (gene_inf$end_position-gene_inf$start_position)
    SampMetrics[[i]][["GeneSizes"]] <- gene_inf$size
    SampMetrics[[i]][["Entrez"]] <- gene_inf$entrezgene_id

     # top percentile genes
    test_bools_top["genes"] <- NULL
    f1_top <- genefilter::kOverA(1, 0.9)
    flist_top <- genefilter::filterfun(f1_top)
    keep_top <- genefilter::genefilter(test_bools_top, flist_top)
    SampMetrics[[i]][["ENSEMBL_top"]] <- gnames[keep_top]
    SampMetrics[[i]][["Entrez_top"]] <- ent[keep_top]
    SampMetrics[[i]][["GeneSymbol_top"]] <- symb[keep_top]
  }
  return(SampMetrics)
}

filterCounts <- function(SampMetrics, technique, filt_options, mart, cell_types) {
  switch(technique,
         cpm = CPM_filter(SampMetrics, filt_options, cell_types),
         zFPKM = zFPKM_filter(SampMetrics, filt_options, mart, cell_types),
         quantile = TPM_quant_filter(SampMetrics, filt_options, cell_types))
}

save_bulk_tests <- function(cmat_file, config_file, out_file, pos_rep=0.5,
                              pos_samp=0.5, technique="quantile", quantile=0.9,
                              min_count=10, top_percentile=25, gene_format="ensembl",
                              species_dataset="hsapiens_gene_ensembl") {

  # condense filter options
  filt_options <- list()
  if ( exists("pos_rep") ) {
    filt_options$pos_rep <- pos_rep
  } else {
    filt_options$pos_rep <- 0.2
  }
  if ( exists("pos_samp") ) {
    filt_options$pos_samp <- pos_samp
  } else {
    filt_options$pos_samp <- 0.5
  }
  if ( exists("quantile") ) {
    filt_options$percentile <- quantile
  } else {
    filt_options$percentile <- 50
  }
  if ( exists("min_count") ) {
    filt_options$min_count <- min_count
  } else {
    filt_options$min_count <- 1
  }
  if ( exists("top_percentile") ) {
    filt_options$top_percentile <- top_percentile
  } else {
    filt_options$top_percentile <- 90
  }

  SampMetrics <- readCountMatrix(cmat_file, config_file)
  #View(SampMetrics[[1]][["CountMatrix"]])
  if ( toupper(species_dataset)=="HUMAN" ) {
    species_dataset <- "hsapiens_gene_ensembl"
  } else if ( toupper(species_dataset)=="MOUSE" ) {
    species_dataset <- "mmusculus_gene_ensembl"
  }
  mart <- useDataset(species_dataset, useMart("ensembl"))
  if ( toupper(gene_format)=="ENSEMBL" ) {
    gene_format <- "ensembl_gene_id"
  } else if (toupper(gene_format)=="SYMBOL") {
    gene_format <- "hgnc_symbol"
  }
  for ( i in 1:length(SampMetrics) ) {
    gene_list <- rownames(SampMetrics[[i]][["CountMatrix"]])
    if ( species_dataset!="hsapiens_gene_ensembl" ) {
          gene_info <- getBM(filters=gene_format,
                       attributes=c("ensembl_gene_id", "entrezgene_id",
                                    "start_position", "end_position"),
                       values=gene_list, mart=mart)
    } else {
           gene_info <- getBM(filters=gene_format,
                       attributes=c("ensembl_gene_id", "hgnc_symbol", "entrezgene_id",
                                    "start_position", "end_position"),
                       values=gene_list, mart=mart)
    }
    gene_info$size <- (gene_info$end_position-gene_info$start_position)
    cdf <- data.frame(SampMetrics[[i]][["CountMatrix"]])
    if ( gene_format=="ensembl_gene_id" ) {
      cdf$ensembl_gene_id <- rownames(SampMetrics[[i]][["CountMatrix"]])
      cdf <- merge(cdf, gene_info, by="ensembl_gene_id")
    } else if ( gene_format=="hgnc_symbol") {
      cdf$hgnc_symbol <- rownames(SampMetrics[[i]][["CountMatrix"]])
      cdf <- merge(cdf, gene_info, by="hgnc_symbol")
    }
    cdf["ensembl_gene_id"] <- NULL
    if ( species_dataset=="hsapiens_gene_ensembl" ) {
      cdf["hgnc_symbol"] <- NULL
    }
    cdf["size"] <- NULL
    cdf["start_position"] <- NULL
    cdf["end_position"] <- NULL
    cdf["entrezgene_id"] <- NULL
    SampMetrics[[i]][["CountMatrix"]] <- cdf
    SampMetrics[[i]][["GeneSizes"]] <- gene_info$size
    SampMetrics[[i]][["ENSEMBL"]] <- gene_info$ensembl_gene_id
    SampMetrics[[i]][["GeneSymbol"]] <- gene_info$hgnc_symbol
    SampMetrics[[i]][["Entrez"]] <- gene_info$entrezgene_id

    rm(gene_info)
  }
  entrez_all <- SampMetrics[[1]][["Entrez"]]
  SampMetrics <- filterCounts(SampMetrics, technique, filt_options, mart, cell_types)
  expressedGenes <- c()
  topGenes <- c()
  for ( i in 1:length(SampMetrics) ) {
    expressedGenes <- c(expressedGenes, SampMetrics[[i]][["Entrez"]])
    topGenes <- c(topGenes, SampMetrics[[i]][["Entrez_top"]])
  }
  expMat <- as.data.frame(table(expressedGenes))
  topMat <- as.data.frame(table(topGenes))
  nc <- length(SampMetrics)
  expMat <- cbind(expMat, "Prop"=expMat$Freq/nc)
  topMat <- cbind(topMat, "Prop"=topMat$Freq/nc)
  SampMetrics[["ExpressionMatrix"]] <- expMat
  SampMetrics[["TopMatrix"]] <- topMat
  SampMetrics[["ExpressedGenes"]] <- as.character(expMat$expressedGenes[expMat$Prop>=pos_samp])
  SampMetrics[["TopGenes"]] <- as.character(topMat$topGenes[topMat$Prop>=pos_samp])
  write_table <- data.frame(entrez_all[!is.na(entrez_all)])
  write_table <- cbind(write_table, rep(0, nrow(write_table)))
  write_table <- cbind(write_table, rep(0, nrow(write_table)))
  for ( i in 1:nrow(write_table) ) {
    if (write_table$entrez_all[i] %in% SampMetrics$ExpressedGenes) {
      write_table[i,2] <- 1
    }
    if (write_table$entrez_all[i] %in% SampMetrics$TopGenes) {
      write_table[i,3] <- 1
    }
  }
  header <- c("ENTREZ_GENE_ID", "expressed", "top")
  write_table <- rbind(header, write_table)
  View(write_table)
  write.table(write_table, out_file, sep=",", row.names=FALSE, col.names=FALSE)
}
"""

bulkio = SignatureTranslatedAnonymousPackage(string, "bulkio")

def load_bulk_supplementary_data(suppfilename):
    if suppfilename=="None":
        return "None"
    suppFullPath = os.path.join(configs.rootdir, 'data', suppfilename)
    supplements = pd.read_csv(suppFullPath, header=0)
    print(supplements)
    supplements = supplements[supplements['SampleName'].str.match("FILENAME")]
    print(supplements)
    Bulk = supplements['InsertSize']
    print(Bulk)
    return Bulk

# Read data from csv files
def load_bulk_tests(Bulk):
    try:
        if not Bulk:
            return None
    except:
        print("bulk exists")

    tests = []
    datas = []
    for test in list(Bulk):
        print(test)
        if test == 'Statistics':
            break
        test = unidecode.unidecode(test)
        # CHANGED BC I CANT USE DOCKER DONT FORGET TO CHANGE BACK
        #output_path = 'G:/GitHub/New Folder/MADRID/docker/pipelines/py/data' + 'Proteomics_{}.csv'.format(test.strip()))
        output_path = os.path.join(configs.rootdir, 'data', 'Bulk_{}.csv'.format(test.strip()))
        testdata = pd.read_csv(output_path, index_col=False)
        #testdat = testdata.applymap(str)
        testdata.drop_duplicates(inplace=True)
        testdata['ENTREZ_GENE_ID'] = testdata['ENTREZ_GENE_ID'].astype(object)
        testdata.set_index('ENTREZ_GENE_ID', inplace=True, drop=True)
        print(testdata.dtypes)
        print('Test Data Load From {}'.format(output_path))
        tests.append(test.strip().replace('ï', 'i'))
        datas.append(testdata)

    bulk_dict = dict(zip(tests, datas))
    return bulk_dict

def main(argv):
    datafile = 'BulkRNAseqDataMatrix.csv'
    suppfile = 'bulk_data_inputs.csv'
    gene_format = 'ensembl'
    species_dataset = 'human'
    pos_rep = 0.5
    pos_samp = 0.5
    top_percentile = 0.25
    technique = 'quantile'
    quantile = 0.9
    min_count = 10
    top_percentile = 25
    gene_format = "ensembl"
    species_dataset = "human"

    try:
        opts, args = getopt.getopt(argv, "hf:c:g:d:r:s:p:t:q:m:",
                ["datafile=", "suppfile=", "gene_format=", "species_dataset=",
                 "expr_prop_rep=", "expr_prop_samp=", "top_percentile=",
                 "technique=", "quantile=", "min_count="])

    except getopt.GetoptError:
        err_str = """
python3 proteomics_gen.py -f <data file> -c <config file>\n
-g <gene format> -d <species_dataset> -r <replicate proportion>\n
-s <sample proportion> -p <top percentile> -t <filtration technique>\n
-q <cutoff quantile (for quantile technique)> -m <min count (for cpm technique)\n
              """
        print(err_str)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            help_str = """
python3 proteomics_gen.py -f <data file> -c <config file>\n
-g <gene format> -d <species_dataset> -r <replicate proportion>\n
-s <sample proportion> -p <top percentile> -t <filtration technique>\n
-q <cutoff quantile (for quantile technique)> -m <min count (for cpm technique)\n
              """
            print(help_str)
            sys.exit()
        elif opt in ("-f", "--datafile"):
            datafilename = arg
        elif opt in ("-c", "--suppfile"):
            suppfilename = arg
        elif opt in ("-g", "--gene_format"):
            gene_format = arg
        elif opt in ("-d", "--species_dataset"):
            species_dataset = arg
        elif opt in ("-r", "--expr_prop_rep"):
            pos_rep = float(arg)
        elif opt in ("-s", "--expr_prop_samp"):
            pos_samp = float(arg)
        elif opt in ("-p", "--top_percentile"):
            top_percentile = int(arg)
        elif opt in ("-t", "--technique"):
            technique = arg
        elif opt in ("-q", "--quantile"):
            quantile = float(arg)
        elif opt in ("-m", "--min_count"):
            min_count = int(arg)
    print('Data file is "{}"'.format(datafilename))
    print('Supplementary Data file is "{}"'.format(suppfilename))

    bulk_config_filepath = os.path.join(configs.rootdir, "data", suppfile)
    bulk_input_filepath = os.path.join(configs.rootdir, "data", datafile)
    #if not os.path.isfile(prote_data_filepath):
    #    proteomics_data.to_csv(prote_data_filepath, index_label='ENTREZ_GENE_ID')

    config_df = pd.read_csv(bulk_config_filepath)
    model_name = "".join(["Bulk_",config_df["InsertSize"][0] ,".csv"])
    bulk_output_filepath = os.path.join(configs.rootdir, "data", model_name)
    print(bulk_output_filepath)

    bulkio.save_bulk_tests(bulk_input_filepath, bulk_config_filepath,
                           bulk_output_filepath, pos_rep=pos_rep,
                           pos_samp=pos_samp, technique=technique,
                           quantile=quantile, min_count=min_count,
                           top_percentile=top_percentile,
                           gene_format=gene_format,
                           species_dataset=species_dataset)
    print("Test data saved to " + bulk_output_filepath)
    # save proteomics data by test
    #proteomics_dict, testdata_dict = save_proteomics_tests(Proteomics, proteomics_data, expr_prop, percentile)
    return True


if __name__ == "__main__":
   main(sys.argv[1:])
