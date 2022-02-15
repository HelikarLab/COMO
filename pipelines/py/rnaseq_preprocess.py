#!/usr/bin/python3

from bioservices import BioDBNet
import pandas as pd
from project import configs
import re
import os, time, sys
import argparse
from rpy2.robjects.packages import importr, SignatureTranslatedAnonymousPackage
import glob

# import R libraries
tidyverse = importr("tidyverse")

# read and translate R functions
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


def create_counts_matrix(tissue_name):
    """
    Create a counts matrix by reading gene counts tables in MADRID_inputs/<tissue name>/<study number>/geneCounts/
    Uses R in backend (generate_counts_matrix.R)
    """
    input_dir = os.path.join(configs.rootdir, 'data', 'MADRID_input', tissue_name)
    print(f"Looking for STAR gene count tables in '{input_dir}'")
    matrix_output_dir = os.path.join(configs.rootdir, 'data', 'data_matrices', tissue_name)
    print(f"Creating Counts Matrix for '{tissue_name}'")
    # call generate_counts_matrix.R to create count matrix from MADRID_input folder
    generate_counts_matrix_io.generate_counts_matrix_main(input_dir, matrix_output_dir)

    return


def create_config_sheet(tissue_name, excel_writer):
    """
    Create configuration sheet at /work/data/config_sheets/rnaseq_data_inputs_auto.xlsx
    based on the gene counts matrix. If using zFPKM normalization technique, fetch mean fragment lengths from
    /work/data/MADRID_inputs/<tissue name>/<study number>/fragmentSizes/
    """
    gene_counts_glob = os.path.join(configs.rootdir, "data", "MADRID_input", tissue_name, "geneCounts", "*", "*.tab")
    gene_counts_files = glob.glob(gene_counts_glob, recursive=True)

    out_df = pd.DataFrame(columns=['SampleName', 'FragmentLength', 'Layout', 'Strand', 'Group'])

    for gcfilename in gene_counts_files:
        try:
            label = re.findall(r"S[1-999]R[1-999]r?[1-999]?", gcfilename)[0]

        except IndexError:
            print(f"\nfilename of {gcfilename} is not valid. Should be 'tissueName_SXRYrZ.tab', where X is the "
                  "study/batch number, Y is the replicate number, and Z is the run number. If not a multi-run sample, "
                  "exclude 'rZ' from the filename.")
            sys.exit()

        study_number = re.findall(r"S[1-999]", label)[0]
        run = re.findall(r"r[1-999]", label)
        multi_flag = 0

        if len(run) > 1:
            if run[0] != "r1":
                continue
            else:
                label_glob = study_number + "R*" + "r*"
                runs = [run for run in gene_counts_files if re.search(label_glob, run)]
                multi_flag = 1
                frag_files = []

                for r in runs:
                    r_label = re.findall(r"r?[1-999]?", r)[0]
                    R_label = re.findall(r"R?[1-999]?", r)[0]
                    frag_filename = "".join([tissue_name, "_", study_number, R_label, r_label, "_fragment_size.txt"])
                    frag_files.append(os.path.join(configs.rootdir, "data", "MADRID_input", tissue_name,
                                                   "fragmentSizes", study_number, frag_filename))

        layout_file = tissue_name + "_" + label + "_layout.txt"
        strand_file = tissue_name + "_" + label + "_strandedness.txt"
        frag_file = tissue_name + "_" + label + "_fragment_size.txt"

        tissue_path = os.path.join(configs.rootdir, "data", "MADRID_input", tissue_name)
        layout_path = os.path.join(tissue_path, "layouts", "*", layout_file)
        strand_path = os.path.join(tissue_path, "strandedness", "*", strand_file)
        frag_path = os.path.join(tissue_path, "fragmentSizes", "*", frag_file)

        layout_glob = glob.glob(layout_path, recursive=False)
        strand_glob = glob.glob(strand_path, recursive=False)
        frag_glob = glob.glob(frag_path, recursive=False)

        # Get layout
        if len(layout_glob) < 1:
            print(f"\nNo layout file found for {label}, writing as 'UNKNOWN', this should be defined by user if using "
                  "zFPKM or rnaseq_gen.py will not run")
            layout = "UNKNOWN"
        elif len(layout_glob) > 1:
            print(f"\nMultiple matching layout files for {label}, make sure there is only one copy for each replicate "
                  "in MADRID_input")
            sys.exit()
        else:
            with open(layout_glob[0]) as file:
                layout = file.read().strip()

        # Get strandedness
        if len(strand_glob) < 1:
            print(f"\nNo strandedness file found for {label}, writing as 'UNKNOWN' This will not interfere with the "
                  "analysis since you have already set rnaseq_preprocess.py to infer the strandedness when writing "
                  "the counts matrix")
            strand = "UNKNOWN"
        elif len(strand_glob) > 1:
            print(f"\nMultiple matching strandedness files for {label}, make sure there is only one copy for each "
                  "replicate in MADRID_input")
            sys.exit()
        else:
            with open(strand_glob[0]) as file:
                strand = file.read().strip()

        # Get fragment length
        if len(frag_glob) < 1:
            print(f"\nNo fragment file found for {label}, writing as 'UNKNOWN' This must be defined by the user in "
                  "order to use zFPKM normalization")
            strand = "UNKNOWN"
        elif len(frag_glob) > 1:
            print(f"\nMultiple matching fragment length files for {label}, make sure there is only one copy for each "
                  "replicate in MADRID_input")
            sys.exit()
        else:

            if layout == "single-end":
                mean_fragment_size = 0
            else:
                if not multi_flag:
                    frag_df = pd.read_table(frag_glob[0], low_memory=False)
                    frag_df['meanxcount'] = frag_df['frag_mean'] * frag_df['frag_count']
                    mean_fragment_size = sum(frag_df['meanxcount'] / sum(frag_df['frag_count']))

                else:
                    mean_fragment_sizes = np.array([])
                    library_sizes = np.array([])
                    for ff in frag_files:
                        frag_df = pd.read_table(ff, low_memory=False)
                        frag_df['meanxcount'] = frag_df['frag_mean']*frag_df['frag_count']
                        mean_fragment_size = sum(frag_df['meanxcount']/sum(frag_df['frag_count']))
                        mean_fragment_sizes.append(mean_fragment_size)
                        library_sizes.append(sum(frag_df['frag_count']))

                    mean_fragment_size = sum(mean_fragment_sizes * library_sizes) / sum(library_sizes)

        label = "_".join([tissue_name, re.findall(r"S[1-999]R[1-999]", label)[0]])  # remove run number if there

        new_row = pd.DataFrame({'SampleName': [label],
                                'FragmentLength': [mean_fragment_size],
                                'Layout': [layout],
                                'Strand': [strand],
                                'Group': [study_number]})

        out_df = out_df.append(new_row)

    out_df.to_excel(excel_writer, sheet_name=tissue_name, header=True, index=False)

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


def handle_tissue_batch(tissue_names, mode, form, taxon_id):
    """
    Handle iteration through each tissue type and create files according to flag used (config, matrix, info)
    """
    config_filename = os.path.join(configs.rootdir, "data", "config_sheets", "rnaseq_data_inputs_auto.xlsx")
    if mode == "config":
        writer = pd.ExcelWriter(config_filename)

    for tissue_name in tissue_names:
        tissue_name = tissue_name.strip(" ")
        print(f"Preprocessing {tissue_name}")
        gene_output_dir = os.path.join(configs.rootdir, "data", "results", tissue_name)
        matrix_output_dir = os.path.join(configs.rootdir, "data", "data_matrices", tissue_name)
        os.makedirs(gene_output_dir, exist_ok=True)
        os.makedirs(matrix_output_dir, exist_ok=True)
        print('Gene info output directory is "{}"'.format(gene_output_dir))

        if mode in ["matrix", "config"]:
            create_counts_matrix(tissue_name)

        count_matrix_file = os.path.join(matrix_output_dir, ("gene_counts_matrix_" + tissue_name + ".csv"))

        if mode == "config":
            create_config_sheet(tissue_name, writer)

        create_gene_info_file(tissue_name, count_matrix_file, form, taxon_id, gene_output_dir)

    if mode == "config":
        writer.close()

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

    handle_tissue_batch(tissue_names, mode, form, taxon_id)

    return


if __name__ == "__main__":
    print(sys.argv)
    main(sys.argv[1:])
