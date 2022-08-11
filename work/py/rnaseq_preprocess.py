#!/usr/bin/python3

import asyncio
from bioservices import BioDBNet
import pandas as pd
from project import configs
import re
import os, time, sys, shutil
import argparse
from rpy2.robjects.packages import importr, SignatureTranslatedAnonymousPackage
import glob
import numpy as np

# import R libraries
tidyverse = importr("tidyverse")

# read and translate R functions
f = open(os.path.join(configs.rootdir, "py", "rscripts", "generate_counts_matrix.R"), "r")
string = f.read()
f.close()
generate_counts_matrix_io = SignatureTranslatedAnonymousPackage(string, 'generate_counts_matrix_io')


async def _async_fetch_gene_info(biodbnet: BioDBNet, **kwargs) -> pd.DataFrame:
    input_db: str = kwargs["input_db"]
    output_db: list[str] = kwargs["output_db"]
    input_values: list[str] = kwargs["input_values"]
    taxon_id: int = kwargs["taxon_id"]
    delay: int = kwargs["delay"]
    event_loop = kwargs["event_loop"]

    # print(f"Getting {current_index}:{current_index + len(input_values)}")
    database_convert = await event_loop.run_in_executor(
        None,  # Defaults to ThreadPoolExecutor, uses threads instead of processes. No need to modify
        biodbnet.db2db,  # The function to call
        input_db,       # The following are arguments passed to the function
        output_db,
        input_values,
        taxon_id
    )

    # If the above db2db conversion didn't work, try again until it does
    while not isinstance(database_convert, pd.DataFrame):
        print(f"\nToo many requests to BioDBNet, waiting {delay} seconds and trying again.")
        await asyncio.sleep(delay)
        database_convert = await event_loop.run_in_executor(
            None,  # Defaults to ThreadPoolExecutor, uses threads instead of processes. No need to modify
            _async_fetch_gene_info,  # The function to call
            biodbnet,
            kwargs
        )

    return database_convert


def fetch_gene_info(
        input_values: list[str],
        input_db: str = "Ensembl Gene ID",
        output_db: list[str] = None,
        taxon_id: int = 9606,
        delay: int = 5
) -> pd.DataFrame:
    """
    This function returns a dataframe with important gene information for future operations in MADRID.

    Fetch gene information from BioDBNet
    :param input_values: A list of genes in "input_db" format
    :param input_db: The input database to use (default: "Ensembl Gene ID")
    :param output_db: The output format to use (default: ["Gene Symbol", "Gene ID", "Chromosomal Location"])
    :param delay: The delay in seconds to wait before trying again if bioDBnet is busy (default: 15)
    :param taxon_id: The taxon ID to use (default: 9606)

    :return: A dataframe with specified columns as "output_db" (Default is HUGO symbol, Entrez ID, and chromosome start and end positions)
    """
    if output_db is None:
        output_db: list[str] = ["Gene Symbol", "Gene ID", "Chromosomal Location"]

    biodbnet = BioDBNet()
    dataframe_maps: pd.DataFrame = pd.DataFrame([], columns=output_db)
    dataframe_maps.index.name = input_db

    if taxon_id == 9606:
        batch_len = 500
    else:
        batch_len = 300

    # Create a list of tasks to be awaited
    event_loop = asyncio.new_event_loop()
    asyncio.set_event_loop(event_loop)
    async_tasks = []
    for i in range(0, len(input_values), batch_len):
        # Define an upper range of values to take from input_values
        upper_range = min(i + batch_len, len(input_values))

        # Create a task to await
        async_tasks.append(
            _async_fetch_gene_info(
                biodbnet,
                input_db=input_db,
                output_db=output_db,
                input_values=input_values[i:upper_range],
                taxon_id=taxon_id,
                delay=delay,
                event_loop=event_loop
            )
        )

    # Await all the tasks, must use [0] to get completed results. [1] is a set of pending tasks (i.e., empty).
    database_convert = event_loop.run_until_complete(asyncio.wait(async_tasks))[0]
    event_loop.close()  # Close the event loop to free resources

    # Loop over database_convert to concat them into dataframe_maps
    for df in database_convert:
        dataframe_maps = pd.concat([dataframe_maps, df.result()], sort=False)

    return dataframe_maps


def create_counts_matrix(context_name):
    """
    Create a counts matrix by reading gene counts tables in MADRID_inputs/<context name>/<study number>/geneCounts/
    Uses R in backend (generate_counts_matrix.R)
    """
    input_dir = os.path.join(configs.rootdir, 'data', 'MADRID_input', context_name)
    print(f"Looking for STAR gene count tables in '{input_dir}'")
    matrix_output_dir = os.path.join(configs.rootdir, 'data', 'data_matrices', context_name)
    print(f"Creating Counts Matrix for '{context_name}'")
    # call generate_counts_matrix.R to create count matrix from MADRID_input folder
    generate_counts_matrix_io.generate_counts_matrix_main(input_dir, matrix_output_dir)

    return


def create_config_df(context_name):
    """
    Create configuration sheet at /work/data/config_sheets/rnaseq_data_inputs_auto.xlsx
    based on the gene counts matrix. If using zFPKM normalization technique, fetch mean fragment lengths from
    /work/data/MADRID_inputs/<context name>/<study number>/fragmentSizes/
    """
    gene_counts_glob = os.path.join(configs.rootdir, "data", "MADRID_input", context_name, "geneCounts", "*", "*.tab")
    gene_counts_files = glob.glob(gene_counts_glob, recursive=True)

    out_df = pd.DataFrame(columns=['SampleName', 'FragmentLength', 'Layout', 'Strand', 'Group'])

    for gcfilename in gene_counts_files:
        try:
            label = re.findall(r"S[1-9][0-9]?[0-9]?R[1-9][0-9]?[0-9]?r?[1-9]?[0-9]?[0-9]?", gcfilename)[0]

        except IndexError:
            print(f"\nfilename of {gcfilename} is not valid. Should be 'contextName_SXRYrZ.tab', where X is the "
                  "study/batch number, Y is the replicate number, and Z is the run number. If not a multi-run sample, "
                  "exclude 'rZ' from the filename.")
            sys.exit()

        study_number = re.findall(r"S[1-9][0-9]?[0-9]?", label)[0]
        rep_number = re.findall(r"R[1-9][0-9]?[0-9]?", label)[0]

        run = re.findall(r"r[1-9][0-9]?[0-9]?", label)
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
                    r_label = re.findall(r"r[1-9][0-9]?[0-9]?", r)[0]
                    R_label = re.findall(r"R[1-9][0-9]?[0-9]?", r)[0]
                    frag_filename = "".join([context_name, "_", study_number, R_label, r_label, "_fragment_size.txt"])
                    frag_files.append(os.path.join(configs.rootdir, "data", "MADRID_input", context_name,
                                                   "fragmentSizes", study_number, frag_filename))

        layout_file = context_name + "_" + label + "_layout.txt"
        strand_file = context_name + "_" + label + "_strandedness.txt"
        frag_file = context_name + "_" + label + "_fragment_size.txt"
        prep_file = context_name + "_" + label + "_prep_method.txt"

        context_path = os.path.join(configs.rootdir, "data", "MADRID_input", context_name)
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

        # Get preparation method
        if len(prep_glob) < 1:
            print(f"No prep file found for {label}, assuming 'total' as in Total RNA library preparation")
            prep = "total"
        elif len(prep_glob) > 1:
            print(f"\nMultiple matching prep files for {label}, make sure there is only one copy for each "
                  "replicate in MADRID_input")
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
            #strand = "UNKNOWN"
            mean_fragment_size = 100
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
                        print(ff)
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
    Split a counts matrix dataframe into two seperate ones. One for Total RNA library prep, one for mRNA
    """
    matrix_all = pd.read_csv(count_matrix_all)
    matrix_total = matrix_all[["genes"] + [n for n in matrix_all.columns if n in df_total["SampleName"].tolist()]]
    matrix_mrna = matrix_all[["genes"] + [n for n in matrix_all.columns if n in df_mrna["SampleName"].tolist()]]

    return matrix_total, matrix_mrna



def create_gene_info_file(matrix_file_list, form, taxon_id):
    """
    Create gene info file for specified context by reading first column in it's count matrix file at
     /work/data/results/<context name>/gene_info_<context name>.csv
    """

    print(f"Fetching gene info")
    gene_info_file = os.path.join(configs.datadir, "gene_info.csv")
    if os.path.exists(gene_info_file):
        current_df = pd.read_csv(gene_info_file)
        genes = current_df["ensembl_gene_id"].tolist()
    else:
        genes = []

    for mfile in matrix_file_list:
        add_genes = pd.read_csv(mfile)["genes"].tolist()
        genes = list(set(genes + add_genes))

    output_db = ['Ensembl Gene ID', 'Gene Symbol', 'Gene ID', 'Chromosomal Location']
    output_db.remove(form)

    print(f"Gathering {len(genes)} genes. Please wait...")

    gene_info = fetch_gene_info(input_values=genes, input_db=form, output_db=output_db, taxon_id=taxon_id)
    gene_info['start_position'] = gene_info['Chromosomal Location'].str.extract("chr_start: (\d+)")
    gene_info['end_position'] = gene_info['Chromosomal Location'].str.extract("chr_end: (\d+)")
    gene_info.index.rename("ensembl_gene_id", inplace=True)
    gene_info.rename(columns={"Gene Symbol": "hgnc_symbol", "Gene ID": "entrezgene_id"}, inplace=True)
    gene_info.drop(['Chromosomal Location'], axis=1, inplace=True)
    gene_info.to_csv(gene_info_file)
    print(f"Gene Info file written at '{gene_info_file}'")


def handle_context_batch(context_names, mode, form, taxon_id, provided_matrix_file):
    """
    Handle iteration through each context type and create files according to flag used (config, matrix, info)
    """
    trnaseq_config_filename = os.path.join(configs.rootdir, "data", "config_sheets", "trnaseq_data_inputs_auto.xlsx")
    mrnaseq_config_filename = os.path.join(configs.rootdir, "data", "config_sheets", "mrnaseq_data_inputs_auto.xlsx")

    tflag = False  # turn on when any total set is found to prevent writer from being init multiple times or empty
    mflag = False  # turn on when any mrna set is found to prevent writer from being init multiple times or empty

    tmatrix_files = []
    mmatrix_files = []
    # pmatrix_files = []
    for context_name in context_names:
        context_name = context_name.strip(" ")
        print(f"Preprocessing {context_name}")
        gene_output_dir = os.path.join(configs.rootdir, "data", "results", context_name)
        matrix_output_dir = os.path.join(configs.rootdir, "data", "data_matrices", context_name)
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
            if len(tmatrix.columns) > 1:
                tmatrix.to_csv(matrix_path_total, header=True, index=False)
            if len(mmatrix.columns) > 1:
                mmatrix.to_csv(matrix_path_mrna, header=True, index=False)

    if mode == "make":
        if tflag:
            twriter.close()
        if mflag:
            mwriter.close()

        create_gene_info_file(tmatrix_files + mmatrix_files, form, taxon_id)

    else:
        create_gene_info_file(matrix_path_prov, form, taxon_id)

    return


def main(argv):
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
            output requires that the 'MADRID_inputs' folder exists and is correctly structured according to the 
            normalization technique being used. A correctly structured folder can be made using our Snakemake-based
            alignment pipeline at:
            https://github.com/HelikarLab/FastqToGeneCounts""",
        epilog="""
            For additional help, please post questions/issues in the MADRID GitHub repo at
            https://github.com/HelikarLab/MADRID or email babessell@gmail.com""",
    )

    parser.add_argument("-n", "--context-names",
                        type=str,
                        required=True,
                        dest="context_names",
                        help="""Tissue/cell name of models to generate. These names should correspond to the folders
                             in 'MADRID_inputs/' if creating count matrix files, or to
                             'work/data/data_matrices/<context name>/gene_counts_matrix_<context name>.csv' if supplying
                             the count matrix as an imported .csv file. If making multiple models in a batch, then
                             use the format: \"['context1', 'context2', ... etc]\". Note the outer double-quotes and the 
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
                            "Requires a correctly structured MADRID_input folder in /work/data/. Can make one using: "
                            "https://github.com/HelikarLab/FastqToGeneCounts"
                       )

    parser.add_argument("-m", "--matrix",
                        required="--provide-matrix" in argv,  # require if using --provide-matrix flag,
                        dest="provided_matrix_fname",
                        default="SKIP",
                        help="Name of provided counts matrix in "
                              "/work/data/data_matrices/<context name>/<NAME OF FILE>.csv"
                        )

    args = parser.parse_args()
    context_names = args.context_names
    gene_format = args.gene_format
    taxon_id = args.taxon_id
    provide_matrix = args.provide_matrix
    make_matrix = args.make_matrix
    provided_matrix_fname = args.provided_matrix_fname

    # convert to python list
    context_names = context_names.strip("[").strip("]").replace("'", "").replace(" ", "").split(",")

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
    if provide_matrix:
        mode = "provide"
    elif make_matrix:
        mode = "make"

    handle_context_batch(context_names, mode, form, taxon_id, provided_matrix_fname)


if __name__ == "__main__":
    print(sys.argv)
    main(sys.argv[1:])
