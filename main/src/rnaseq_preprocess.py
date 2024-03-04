import argparse
import glob
import os
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from fast_bioservices import BioDBNet, Input, Output, Taxon

import como_utilities
import rpy2_api
from project import Configs

configs = Configs()
r_file_path: Path = Path(configs.root_dir, "src", "rscripts", "generate_counts_matrix.R")


def create_counts_matrix(context_name):
    """
    Create a counts matrix by reading gene counts tables in COMO_input/<context name>/<study number>/geneCounts/
    Uses R in backend (generate_counts_matrix.R)
    """
    input_dir = os.path.join(configs.root_dir, 'data', 'COMO_input', context_name)
    print(f"Looking for STAR gene count tables in '{input_dir}'")
    matrix_output_dir = os.path.join(configs.root_dir, 'data', 'data_matrices', context_name)
    print(f"Creating Counts Matrix for '{context_name}'")
    
    # call generate_counts_matrix.R to create count matrix from COMO_input folder
    rpy2_hook = rpy2_api.Rpy2(r_file_path=r_file_path, data_dir=input_dir, out_dir=matrix_output_dir)
    rpy2_hook.call_function("generate_counts_matrix_main")


def create_config_df(context_name):
    """
    Create configuration sheet at /main/data/config_sheets/rnaseq_data_inputs_auto.xlsx
    based on the gene counts matrix. If using zFPKM normalization technique, fetch mean fragment lengths from
    /work/data/COMO_input/<context name>/<study number>/fragmentSizes/
    """
    gene_counts_glob = os.path.join(configs.root_dir, "data", "COMO_input", context_name, "geneCounts", "*", "*.tab")
    gene_counts_files = glob.glob(gene_counts_glob, recursive=True)
    
    out_df = pd.DataFrame(columns=['SampleName', 'FragmentLength', 'Layout', 'Strand', 'Group'])
    
    for gcfilename in gene_counts_files:
        try:
            # Match S___R___r___
            # \d{1,3} matches 1-3 digits
            # (r\d{1,3})? matches an option "r" followed by three digits
            label = re.findall(r"S\d{1,3}R\d{1,3}(?:r\d{1,3})?", gcfilename)[0]
        
        except IndexError:
            raise IndexError(
                f"\n\nFilename of '{gcfilename}' is not valid. Should be 'contextName_SXRYrZ.tab', where X is the "
                "study/batch number, Y is the replicate number, and Z is the run number."
                "\n\nIf not a multi-run sample, exclude 'rZ' from the filename."
            )
        
        study_number = re.findall(r"S\d{1,3}", label)[0]
        rep_number = re.findall(r"R\d{1,3}", label)[0]
        run = re.findall(r"r\d{1,3}", label)
        
        multi_flag = 0
        if len(run) > 0:
            if run[0] != "r1":
                continue
            else:
                label_glob = study_number + rep_number + "r*"
                runs = [run for run in gene_counts_files if re.search(label_glob, run)]
                multi_flag = 1
                frag_files = []
                
                for r in runs:
                    r_label = re.findall(r"r\d{1,3}", r)[0]
                    R_label = re.findall(r"R\d{1,3}", r)[0]
                    frag_filename = "".join([context_name, "_", study_number, R_label, r_label, "_fragment_size.txt"])
                    frag_files.append(os.path.join(configs.root_dir, "data", "COMO_input", context_name,
                                                   "fragmentSizes", study_number, frag_filename))
        
        layout_file = context_name + "_" + label + "_layout.txt"
        strand_file = context_name + "_" + label + "_strandedness.txt"
        frag_file = context_name + "_" + label + "_fragment_size.txt"
        prep_file = context_name + "_" + label + "_prep_method.txt"
        
        context_path = os.path.join(configs.root_dir, "data", "COMO_input", context_name)
        layout_path = os.path.join(context_path, "layouts", "*", layout_file)
        strand_path = os.path.join(context_path, "strandedness", "*", strand_file)
        frag_path = os.path.join(context_path, "fragmentSizes", "*", frag_file)
        prep_path = os.path.join(context_path, "prepMethods", "*", prep_file)
        
        layout_glob = glob.glob(layout_path, recursive=False)
        strand_glob = glob.glob(strand_path, recursive=False)
        frag_glob = glob.glob(frag_path, recursive=False)
        prep_glob = glob.glob(prep_path, recursive=False)
        
        # Get layout
        if len(layout_glob) < 1:
            print(f"\nNo layout file found for {label}, writing as 'UNKNOWN', this should be defined by user if using "
                  "zFPKM or rnaseq_gen.py will not run")
            layout = "UNKNOWN"
        elif len(layout_glob) > 1:
            print(f"\nMultiple matching layout files for {label}, make sure there is only one copy for each replicate "
                  "in COMO_input")
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
                  "replicate in COMO_input")
            sys.exit()
        else:
            with open(strand_glob[0]) as file:
                strand = file.read().strip()
        
        # Get preparation method
        if len(prep_glob) < 1:
            print(f"No prep file found for {label}, assuming 'total' as in Total RNA library preparation")
            prep = "total"
        elif len(prep_glob) > 1:
            print(f"\nMultiple matching prep files for {label}, make sure there is only one copy for each "
                  "replicate in COMO_input")
            sys.exit()
        else:
            with open(prep_glob[0]) as file:
                prep = file.read().strip().lower()
                if prep not in ["total", "mrna"]:
                    print("prep_method.txt must have either 'mrna' or 'total', or be absent to assume 'total'.")
                    sys.exit()
        
        # Get fragment length
        if len(frag_glob) < 1:
            print(f"\nNo fragment file found for {label}, writing as 'UNKNOWN' This must be defined by the user in "
                  "order to use zFPKM normalization")
            # strand = "UNKNOWN"
            mean_fragment_size = 100
        elif len(frag_glob) > 1:
            print(f"\nMultiple matching fragment length files for {label}, make sure there is only one copy for each "
                  "replicate in COMO_input")
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
                        frag_df = pd.read_table(ff, low_memory=False, sep="\t", on_bad_lines="skip")
                        frag_df['meanxcount'] = frag_df['frag_mean'] * frag_df['frag_count']
                        mean_fragment_size = sum(frag_df['meanxcount'] / sum(frag_df['frag_count']))
                        mean_fragment_sizes = np.append(mean_fragment_sizes, mean_fragment_size)
                        library_sizes = np.append(library_sizes, sum(frag_df['frag_count']))
                    
                    mean_fragment_size = sum(mean_fragment_sizes * library_sizes) / sum(library_sizes)
        
        # label = "_".join([context_name, re.findall(r"S[1-9][0-9]?[0-9]?R[1-9][0-9]?[0-9]?", label)[0]])  # remove run number if there
        label = f"{context_name}_{study_number}{rep_number}"
        
        new_row = pd.DataFrame({'SampleName': [label],
                                'FragmentLength': [mean_fragment_size],
                                'Layout': [layout],
                                'Strand': [strand],
                                'Group': [study_number],
                                'LibraryPrep': [prep]})
        
        out_df = pd.concat([out_df, new_row], sort=True)
        out_df.sort_values('SampleName', inplace=True)
    
    return out_df


def split_config_df(df):
    """
    Split a config dataframe into two seperate ones. One for Total RNA library prep, one for mRNA
    """
    df_t = df[df['LibraryPrep'] == "total"]
    df_m = df[df['LibraryPrep'] == "mrna"]
    
    return df_t, df_m


def split_counts_matrices(count_matrix_all, df_total, df_mrna):
    """
    Split a counts-matrix dataframe into two seperate ones. One for Total RNA library prep, one for mRNA
    """
    matrix_all = pd.read_csv(count_matrix_all)
    matrix_total = matrix_all[["genes"] + [n for n in matrix_all.columns if n in df_total["SampleName"].tolist()]]
    matrix_mrna = matrix_all[["genes"] + [n for n in matrix_all.columns if n in df_mrna["SampleName"].tolist()]]
    
    return matrix_total, matrix_mrna


def create_gene_info_file(matrix_file_list: list[str], input_format: Input, taxon_id):
    """
    Create gene info file for specified context by reading first column in its count matrix file at
     results/<context name>/gene_info_<context name>.csv
    """
    
    print("Fetching gene info")
    gene_info_file = os.path.join(configs.data_dir, "gene_info.csv")
    genes: list[str] = []
    for file in matrix_file_list:
        genes += pd.read_csv(file)["genes"].tolist()
    genes = list(set(genes))
    
    # Create our output database format
    # Do not include values equal to "form"
    # Remove items not equal to `form` because the input database cannot exist as an output database
    output_db: list[Output] = [
        i for i in [
            Output.ENSEMBL_GENE_ID,
            Output.GENE_SYMBOL,
            Output.GENE_ID,
            Output.CHROMOSOMAL_LOCATION
        ]
        if i.value != input_format.value
    ]
    
    biodbnet = BioDBNet()
    gene_info = biodbnet.db2db(
        input_values=genes,
        input_db=input_format,
        output_db=output_db,
        taxon=taxon_id,
    )
    
    gene_info.rename(columns={Output.ENSEMBL_GENE_ID.value: "ensembl_gene_id"}, inplace=True)
    gene_info['start_position'] = gene_info['Chromosomal Location'].str.extract(r"chr_start: (\d+)")
    gene_info['end_position'] = gene_info['Chromosomal Location'].str.extract(r"chr_end: (\d+)")
    gene_info.rename(columns={"Gene Symbol": "hgnc_symbol", "Gene ID": "entrezgene_id"}, inplace=True)
    gene_info.drop(['Chromosomal Location'], axis=1, inplace=True)
    gene_info.to_csv(gene_info_file, index=False)
    print(f"Gene Info file written at '{gene_info_file}'")


def handle_context_batch(context_names, mode, input_format: Input, taxon_id, provided_matrix_file):
    """
    Handle iteration through each context type and create files according to flag used (config, matrix, info)
    """
    trnaseq_config_filename = os.path.join(configs.root_dir, "data", "config_sheets", "trnaseq_data_inputs_auto.xlsx")
    mrnaseq_config_filename = os.path.join(configs.root_dir, "data", "config_sheets", "mrnaseq_data_inputs_auto.xlsx")
    
    tflag = False  # turn on when any total set is found to prevent writer from being init multiple times or empty
    mflag = False  # turn on when any mrna set is found to prevent writer from being init multiple times or empty
    
    print(f"Found {len(context_names)} contexts to process: {', '.join(context_names)}")
    
    tmatrix_files = []
    mmatrix_files = []
    for context_name in context_names:
        context_name = context_name.strip(" ")
        print(f"Preprocessing {context_name}")
        gene_output_dir = os.path.join(configs.root_dir, "data", "results", context_name)
        matrix_output_dir = os.path.join(configs.root_dir, "data", "data_matrices", context_name)
        os.makedirs(gene_output_dir, exist_ok=True)
        os.makedirs(matrix_output_dir, exist_ok=True)
        print('Gene info output directory is "{}"'.format(gene_output_dir))
        
        matrix_path_all = os.path.join(matrix_output_dir, ("gene_counts_matrix_full_" + context_name + ".csv"))
        matrix_path_total = os.path.join(matrix_output_dir, ("gene_counts_matrix_total_" + context_name + ".csv"))
        matrix_path_mrna = os.path.join(matrix_output_dir, ("gene_counts_matrix_mrna_" + context_name + ".csv"))
        matrix_path_prov = os.path.join(matrix_output_dir, provided_matrix_file)
        
        if mode == "make":
            create_counts_matrix(context_name)
            # TODO: warn user or remove samples that are all 0 to prevent density plot error in zFPKM
            df = create_config_df(context_name)
            df_t, df_m = split_config_df(df)
            
            if not df_t.empty:
                if not tflag:
                    tflag = True
                    twriter = pd.ExcelWriter(trnaseq_config_filename)
                
                tmatrix_files.append(matrix_path_total)
                df_t.to_excel(twriter, sheet_name=context_name, header=True, index=False)
            
            if not df_m.empty:
                if not mflag:
                    mflag = True
                    mwriter = pd.ExcelWriter(mrnaseq_config_filename)
                
                mmatrix_files.append(matrix_path_mrna)
                df_m.to_excel(mwriter, sheet_name=context_name, header=True, index=False)
            
            tmatrix, mmatrix = split_counts_matrices(matrix_path_all, df_t, df_m)
            if len(tmatrix.columns) >= 1:
                tmatrix.to_csv(matrix_path_total, header=True, index=False)
            if len(mmatrix.columns) >= 1:
                mmatrix.to_csv(matrix_path_mrna, header=True, index=False)
    
    if mode == "make":
        if tflag:
            twriter.close()
        if mflag:
            mwriter.close()
        
        create_gene_info_file(tmatrix_files + mmatrix_files, input_format, taxon_id)
    
    else:
        matrix_files: list[str] = como_utilities.stringlist_to_list(matrix_path_prov)
        create_gene_info_file(matrix_files, input_format, taxon_id)


def parse_args(argv):
    """
    Parse arguments to rnaseq_preprocess.py, create a gene info files for each provided context at:
    /work/data/results/<context name>/gene_info_<context name>.csv.

     If using --info-matrix or --info-matrix-config:
    create gene count matrix file at /work/data/data_matrices/<context name>/gene_counts_matrix_<context name>.csv,

    If using --info-matrix-config:
    create config file at /work/data/config_sheets/rnaseq_data_inputs_auto.xlsx
    """
    parser = argparse.ArgumentParser(
        prog="rnaseq_preprocess.py",
        description="""
            Fetches additional gene information from a provided matrix or gene counts, or optionally creates this
            matrix using gene count files obtained using STAR aligner. Creation of counts matrix from STAR aligner
            output requires that the 'COMO_input' folder exists and is correctly structured according to the
            normalization technique being used. A correctly structured folder can be made using our Snakemake-based
            alignment pipeline at:
            https://github.com/HelikarLab/FastqToGeneCounts""",
        epilog="""
            For additional help, please post questions/issues in the MADRID GitHub repo at
            https://github.com/HelikarLab/MADRID or email babessell@gmail.com""",
    )
    
    parser.add_argument("-n", "--context-names",
                        type=str,
                        nargs="+",
                        required=True,
                        dest="context_names",
                        help="""Tissue/cell name of models to generate. These names should correspond to the folders
                             in 'COMO_input/' if creating count matrix files, or to
                             'work/data/data_matrices/<context name>/gene_counts_matrix_<context name>.csv' if supplying
                             the count matrix as an imported .csv file. If making multiple models in a batch, then
                             use the format: "context1 context2 context3". """
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
    
    group.add_argument("-p", "--provide-matrix",
                       action="store_true",
                       required=False,
                       default=False,
                       dest="provide_matrix",
                       help="Provide your own count matrix. Requires additional argument '--matrix' which is .csv file "
                            "where colnames are sample names (in contextName_SXRY format) and rownames are genes in "
                            "in format specified by --gene-format"
                       )  # would be nice if this was a directory full matrices in case you want to do in batches
    
    group.add_argument('-y', "--create-matrix",
                       action="store_true",
                       required=False,
                       default=False,
                       dest="make_matrix",
                       help="Flag for if you want to make a counts matrix, but not a config file. "
                            "Requires a correctly structured COMO_input folder in /work/data/. Can make one using: "
                            "https://github.com/HelikarLab/FastqToGeneCounts"
                       )
    
    parser.add_argument("-m", "--matrix",
                        required="--provide-matrix" in argv,  # require if using --provide-matrix flag,
                        dest="provided_matrix_fname",
                        default="SKIP",
                        help="Name of provided counts matrix in "
                             "/work/data/data_matrices/<context name>/<NAME OF FILE>.csv"
                        )
    
    args = parser.parse_args(argv)
    args.context_names = como_utilities.stringlist_to_list(args.context_names)
    
    return args


def main(argv):
    args = parse_args(argv)
    
    if args.gene_format.upper() in ["ENSEMBL", "ENSEMBLE", "ENSG", "ENSMUSG", "ENSEMBL ID", "ENSEMBL GENE ID"]:
        gene_format_database: Input = Input.ENSEMBL_GENE_ID
    
    elif args.gene_format.upper() in ["HGNC SYMBOL", "HUGO", "HUGO SYMBOL", "SYMBOL", "HGNC", "GENE SYMBOL"]:
        gene_format_database: Input = Input.GENE_SYMBOL
    
    elif args.gene_format.upper() in ["ENTREZ", "ENTRES", "ENTREZ ID", "ENTREZ NUMBER" "GENE ID"]:
        gene_format_database: Input = Input.GENE_ID
    
    else:  # provided invalid gene format
        print("Gene format (--gene_format) is invalid")
        print("Accepts 'Ensembl', 'Entrez', and 'HGNC symbol'")
        print(f"You provided: {args.gene_format}")
        sys.exit()
    
    # handle species alternative ids
    if isinstance(args.taxon_id, str):
        if args.taxon_id.upper() == "HUMAN" or args.taxon_id.upper() == "HOMO SAPIENS":
            taxon_id = Taxon.HOMO_SAPIENS
        elif args.taxon_id.upper() == "MOUSE" or args.taxon_id.upper() == "MUS MUSCULUS":
            taxon_id = Taxon.MUS_MUSCULUS
        else:
            print("--taxon-id must be either an integer, or accepted string ('mouse', 'human')")
            sys.exit(1)
    
    elif isinstance(args.taxon_id, int):
        taxon_id = args.taxon_id
    else:
        print("--taxon-id must be either an integer, or accepted string ('mouse', 'human')")
        sys.exit(1)
    
    # use mutually exclusive flag to set mode which tells which files to generate
    if args.provide_matrix:
        mode = "provide"
    elif args.make_matrix:
        mode = "make"
    else:
        print("--provide-matrix or --create-matrix must be set")
        sys.exit(1)
    
    handle_context_batch(
        context_names=args.context_names,
        mode=mode,
        input_format=gene_format_database,
        taxon_id=taxon_id,
        provided_matrix_file=args.provided_matrix_fname
    )


if __name__ == "__main__":
    main(sys.argv[1:])
