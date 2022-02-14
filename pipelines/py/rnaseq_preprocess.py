#!/usr/bin/python3

from bioservices import BioDBNet
import pandas as pd
from project import configs
import re
import os, time, sys
# import getopt
import argparse
from rpy2.robjects.packages import importr, SignatureTranslatedAnonymousPackage

# limma = importr("limma")
tidyverse = importr("tidyverse")
# edgeR = importr("edgeR")
# genefilter = importr("genefilter")
# biomaRt = importr("biomaRt")
# sjmisc = importr("sjmisc")

# automatically convert ryp2 dataframe to Pandas dataframe
f = open(os.path.join("rscripts", "generate_counts_matrix.R"), "r")
string = f.read()
f.close()

generate_counts_matrix_io = SignatureTranslatedAnonymousPackage(string, 'generate_counts_matrix_io')


def fetch_gene_info(input_values, input_db="Ensembl Gene ID",
                    output_db=["Gene Symbol", "Gene ID", "Chromosomal Location"],
                    delay=15, taxon_id=9606):
    """
    Returns a dataframe with important gene information for future operations in MADRID.

    Fetch gene information from BioDbNet, takes 'input_values' (genes) in format of 'input_db' (default, Ensembl) and
    ond returns dataframe with specified columns as 'output_db' (default is HUGO symbol, Entrez ID, and chromosome
    chromosomal start and end positions).
    """

    s = BioDBNet()
    df_maps = pd.DataFrame([], columns=output_db)
    df_maps.index.name = input_db
    print(input_db)
    i = 0
    batch_len = 300
    print(f"Total Genes to Retrieve: {len(input_values)}")
    while i < len(input_values):
        upper_range = min(i + batch_len, len(input_values))
        print(f"retrieve {i}:{upper_range}")
        df_test = s.db2db(input_db, output_db, input_values[i:upper_range], taxon_id)
        if isinstance(df_test, pd.DataFrame):
            df_maps = pd.concat([df_maps, df_test], sort=False)
        elif df_test == '414':
            print(f"bioDBnet busy, trying again in {delay} seconds")
            time.sleep(delay)
            continue
        i += batch_len
    return df_maps

def create_counts_matrix(tissue_name, technique):
    """
    Create a counts matrix by reading gene counts tables in MADRID_inputs/<tissue name>/<study number>/geneCounts/
    Uses R in backend (generate_counts_matrix.R)
    """
    input_dir = os.path.join(configs.rootdir, 'data', 'MADRID_input', tissue_name)
    print(f"Looking for STAR gene count tables in '{input_dir}'")
    matrix_output_dir = os.path.join(configs.rootdir, 'data', 'data_matrices', tissue_name)
    print(f"Creating Counts Matrix for '{tissue_name}'")
    # call genCountMatrix_old.R to create count matrix from MADRID_input folder
    generate_counts_matrix_io.generate_counts_matrix_main(input_dir, matrix_output_dir)

    return

def create_config_sheet(tissue_name):
    """
    Create configuration sheet at /work/data/config_sheets/rnaseq_data_inputs_auto.xlsx
    based on the gene counts matrix. If using zFPKM normalization technique, fetch mean fragment lengths from
    /work/data/MADRID_inputs/<tissue name>/<study number>/fragmentSizes/
    """

    init_df = pd.DataFrame(columns=['SampleName', 'FragmentSize', 'Layout', 'Strand', 'Group'])

    return

def create_gene_info_file(tissue_name, count_matrix_file, form, taxon_id, gene_output_dir):
    """
    Create gene info file for specified tissue by reading first column in it's count matrix file at
     /work/data/results/<tissue name>/gene_info_<tissue name>.csv
    """

    print(f"Fetching gene info using genes in '{count_matrix_file}'")
    genes = pd.read_csv(count_matrix_file)['genes'].to_list()
    output_db = ['Ensembl Gene ID', 'Gene Symbol', 'Gene ID', 'Chromosomal Location']
    output_db.remove(form)
    gene_info = fetch_gene_info(genes, input_db=form, output_db=output_db, taxon_id=taxon_id)
    gene_info['start_position'] = gene_info['Chromosomal Location'].str.extract("chr_start: (\d+)")
    gene_info['end_position'] = gene_info['Chromosomal Location'].str.extract("chr_end: (\d+)")
    gene_info.index.rename("ensembl_gene_id", inplace=True)
    gene_info.rename(columns={"Gene Symbol": "hgnc_symbol", "Gene ID": "entrezgene_id"}, inplace=True)
    gene_info.drop(['Chromosomal Location'], axis=1, inplace=True)
    gene_info_file = os.path.join(gene_output_dir, ("gene_info_" + tissue_name + ".csv"))
    gene_info.to_csv(gene_info_file)
    print(f"Gene Info file written at '{gene_info_file}'")

    return

def handle_tissue_batch(tissue_names, mode, technique, form, taxon_id):
    """
    Handle iteration through each tissue type and create files according to flag used (config, matrix, info)
    """

    config_file_path = os.path.join(configs.rootdir, "data", "config_sheets", "rnaseq_data_inputs.xlsx")
    for tissue_name in tissue_names:
        tissue_name = tissue_name.strip(" ")
        gene_output_dir = os.path.join(configs.rootdir, "data", "results", tissue_name)
        matrix_output_dir = os.path.join(configs.rootdir, "data", "data_matrices", tissue_name)
        os.makedirs(gene_output_dir, exist_ok=True)
        os.makedirs(matrix_output_dir, exist_ok=True)
        print('Gene info output directory is "{}"'.format(gene_output_dir))
        print('Active gene determination technique is "{}"'.format(technique))

        if mode in ["matrix", "config"]:
            create_counts_matrix(tissue_name, technique)

        count_matrix_file = os.path.join(matrix_output_dir, ("gene_counts_matrix_" + tissue_name + ".csv"))

        if mode == "config":
            init_df = pd.DataFrame(columns=['SampleName', 'FragmentSize', 'Layout', 'Strand', 'Group'])
            init_df.to_excel(config_file_path, sheet_name=tissue_name)
            create_config_sheet('temp')

        create_gene_info_file(tissue_name, count_matrix_file, form, taxon_id, gene_output_dir)

        return


def main(argv):
    """
    Parse arguments to rnaseq_preprocess.py, create a gene info files for each provided tissue at:
    /work/data/results/<tissue name>/gene_info_<tissue name>.csv.

     If using --info-matrix or --info-matrix-config:
    create gene count matrix file at /work/data/data_matrices/<tissue name>/gene_counts_matrix_<tissue name>.csv,

    If using --info-matrix-config:
    create config file at /work/data/config_sheets/rnaseq_data_inputs_auto.xlsx
    """

    parser = argparse.ArgumentParser(
        prog="rnaseq_preprocess.py",
        description="""
            Fetches additional gene information from a provided matrix or gene counts, or optionally creates this
            matrix using gene count files obtained using STAR aligner. Creation of counts matrix from STAR aligner 
            output requires that the 'MADRID_inputs' folder exists and is correctly structured according to the 
            normalization technique being used. A correctly structured folder can be made using our Snakemake-based
            alignment pipeline at:
            https://github.com/HelikarLab/FastqToGeneCounts""",
        epilog="""
            For additional help, please post questions/issues in the MADRID GitHub repo at
            https://github.com/HelikarLab/MADRID or email babessell@gmail.com""",
        usage="python3 $(prog)s [options]"
    )

    parser.add_argument("-n", "--tissue-name",
                        type=str,
                        required=True,
                        dest="tissue_names",
                        help="""Tissue/cell name of models to generate. These names should correspond to the folders
                             in 'MADRID_inputs/' if creating count matrix files, or to
                             'work/data/data_matrices/<tissue name>/gene_counts_matrix_<tissue name>.csv' if supplying
                             the count matrix as an imported .csv file. If making multiple models in a batch, then
                             use the format: \"['tissue1', 'tissue2', ... etc]\". Note the outer double-quotes and the 
                             inner single-quotes are required to be interpreted. This a string, not a python list"""
                        )

    parser.add_argument("-f", "--gene-format",
                        type=str,
                        required=False,
                        default="Ensembl Gene ID",
                        dest="gene_format",
                        help="Format of Genes, accepts 'Ensembl', 'Entrez', or'HGNC symbol'"
                        )

    parser.add_argument("-i", "--taxon-id",
                        required=False,
                        default=9606,
                        dest="taxon_id",
                        help="BioDbNet taxon ID number, also accepts 'human', or 'mouse'"
                        )
    parser.add_argument("-t", "--norm-technique",
                        type=str,
                        required=False,
                        default="tpm-quantile",
                        dest="technique",
                        help="Normalization technique, accepts 'tpm-quantile', 'cpm', or 'zFPKM'"
                        )

    group = parser.add_mutually_exclusive_group(required=True)
    
    group.add_argument("-x", "--info-only",
                       action="store_true",
                       required=False,
                       default=False,
                       dest="info_only",
                       help="Flag if you want to use your own supplied counts matrix at " +
                            "/work/data/data_matrices/<tissue name>/gene_count_matrix_<tissue name>.csv"
                       )

    group.add_argument('-y', "--info-matrix",
                       action="store_true",
                       required=False,
                       default=False,
                       dest="make_matrix",
                       help="""Flag for if you want to make a counts matrix, but not a config file.
                            Requires a correctly structured MADRID_input folder in /work/data/. Can make one using: 
                            https://github.com/HelikarLab/FastqToGeneCounts"""
                       )

    group.add_argument("-z", "--info-matrix-config",
                       action="store_true",
                       required=False,
                       default=False,
                       dest="make_config",
                       help="""Flag if you want to make both a count matrix, and a config file.
                            Requires a correctly structured MADRID_input folder in /work/data/. Can make one using: 
                            https://github.com/HelikarLab/FastqToGeneCounts"""
                       )

    args = parser.parse_args(argv)
    tissue_names = args.tissue_names
    gene_format = args.gene_format
    taxon_id = args.taxon_id
    technique = args.technique
    info_only = args.info_only
    make_matrix = args.make_matrix
    make_config = args.make_config

    tissue_names = tissue_names.strip("[").strip("]").replace("'", "").replace(" ", "").split(",") # convert to py list
    if gene_format.upper() in ["ENSEMBL", "ENSEMBLE", "ENSG", "ENSMUSG", "ENSEMBL ID", "ENSEMBL GENE ID"]:
        form = "Ensembl Gene ID"

    elif gene_format.upper() in ["HGNC SYMBOL", "HUGO", "HUGO SYMBOL", "SYMBOL", "HGNC", "GENE SYMBOL"]:
        form = "Gene Symbol"

    elif gene_format.upper() in ["ENTREZ", "ENTRES", "ENTREZ ID", "ENTREZ NUMBER" "GENE ID"]:
        form = "Gene ID"

    else:  # provided invalid gene format
        print("Gene format (--gene_format) is invalid")
        print("Accepts 'Ensembl', 'Entrez', and 'HGNC symbol'")
        sys.exit()

    # handle species alternative ids
    if type(taxon_id) == str:
        if taxon_id.upper() == "HUMAN" or taxon_id.upper() == "HOMO SAPIENS":
            taxon_id = 9606
        elif taxon_id.upper() == "MOUSE" or taxon_id.upper() == "MUS MUSCULUS":
            taxon_id = 10090
        else:
            print("--taxon-id must be either an integer, or accepted string (\"mouse\", \"human\")")
            sys.exit()
    elif type(taxon_id) != int:
        print("--taxon-id must be either an integer, or accepted string (\"mouse\", \"human\")")
        sys.exit()

    # use mutually exclusive flag to set mode which tells which files to generate
    if info_only:
        mode = "info"
    elif make_matrix:
        mode = "matrix"
    elif make_config:
        mode = "config"

    handle_tissue_batch(tissue_names, mode, technique, form, taxon_id)

    return


if __name__ == "__main__":
    print(sys.argv)
    main(sys.argv[1:])
