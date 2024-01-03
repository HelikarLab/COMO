import argparse
import sys
import os
import pandas as pd
import numpy as np
import instruments
from pathlib import Path

import project
import rpy2_api

configs = project.Configs()

# read and translate R functions
# f = open(os.path.join(configs.rootdir, "src", "rscripts", "protein_transform.R"), "r")
# string = f.read()
# f.close()
# protein_transform_io = SignatureTranslatedAnonymousPackage(string, "protein_transform_io")
r_file_path = Path(configs.rootdir, "src", "rscripts", "protein_transform.R")


# Load Proteomics
def load_proteomics_data(datafilename, context_name):
    """
    Add description......
    """
    dataFullPath = os.path.join(configs.rootdir, "data", "data_matrices", context_name, datafilename)
    print('Data matrix is at "{}"'.format(dataFullPath))
    
    if os.path.isfile(dataFullPath):
        proteomics_data = pd.read_csv(dataFullPath, header=0)
    
    else:
        print("Error: file not found: {}".format(dataFullPath))
        
        return None
    
    # Preprocess data, drop na, duplicate ';' in symbol,
    proteomics_data["symbol"] = proteomics_data["symbol"].astype(str)
    proteomics_data.dropna(subset=["symbol"], inplace=True)
    pluralnames = proteomics_data[proteomics_data["symbol"].str.contains(";") == True]
    
    for idx, row in pluralnames.iterrows():
        names = row["symbol"].split(";")
        rows = []
        
        for name in names:
            rowcopy = row.copy()
            rowcopy["symbol"] = name
            rows.append(rowcopy)
        proteomics_data.drop(index=idx, inplace=True)
        proteomics_data = pd.concat(
            [proteomics_data, pd.DataFrame(rows)], ignore_index=True
        )
    
    proteomics_data.rename(columns={"symbol": "Gene Symbol"}, inplace=True)
    
    return proteomics_data


# read map to convert to entrez
def load_gene_symbol_map(gene_symbols, filename="proteomics_entrez_map.csv"):
    """
    Add descirption....
    """
    filepath = os.path.join(configs.rootdir, "data", "proteomics_entrez_map.csv")
    if os.path.isfile(filepath):
        sym2id = pd.read_csv(filepath, index_col="Gene Symbol")
    else:
        sym2id = instruments.fetch_entrez_gene_id(gene_symbols, input_db="Gene Symbol")
        sym2id.loc[sym2id["Gene ID"] == "-", ["Gene ID"]] = np.nan
        sym2id.to_csv(filepath, index_label="Gene Symbol")
    
    return sym2id[~sym2id.index.duplicated()]


def abundance_to_bool_group(context_name, group_name, abundance_matrix, rep_ratio, hi_rep_ratio, quantile):
    """
    Descrioption....
    """
    output_dir = os.path.join(configs.rootdir, "data", "results", context_name, "proteomics")
    os.makedirs(output_dir, exist_ok=True)
    
    # write group abundances to individual files
    abundance_filepath = os.path.join(
        configs.datadir,
        "results",
        context_name,
        "proteomics",
        "".join(["protein_abundance_", group_name, ".csv"])
    )
    abundance_matrix.to_csv(abundance_filepath, index_label="ENTREZ_GENE_ID")
    
    # Z-tranform
    rpy2_api.Rpy2(
        r_file_path,
        abundance_filepath,
        output_dir,
        group_name
    ).call_function("protein_transform_main")
    # protein_transform_io.protein_transform_main(abundance_filepath, output_dir, group_name)
    
    # Logical Calculation
    thresholds = abundance_matrix.quantile(quantile, axis=0)
    abundance_matrix_nozero = abundance_matrix.replace(0, np.nan)
    thresholds = abundance_matrix_nozero.quantile(quantile, axis=0)
    testbool = pd.DataFrame(0, columns=list(abundance_matrix), index=abundance_matrix.index)
    
    for col in list(abundance_matrix):
        testbool.loc[abundance_matrix[col] > thresholds[col], [col]] = 1
    
    abundance_matrix["pos"] = (abundance_matrix > 0).sum(axis=1) / abundance_matrix.count(axis=1)
    abundance_matrix["expressed"] = 0
    abundance_matrix.loc[(abundance_matrix["pos"] >= rep_ratio), ["expressed"]] = 1
    abundance_matrix["high"] = 0
    abundance_matrix.loc[(abundance_matrix["pos"] >= hi_rep_ratio), ["high"]] = 1
    
    bool_filepath = os.path.join(output_dir, f"bool_prot_Matrix_{context_name}_{group_name}.csv")
    abundance_matrix.to_csv(bool_filepath, index_label="ENTREZ_GENE_ID")


def to_bool_context(context_name, group_ratio, hi_group_ratio, group_names):
    output_dir = os.path.join(configs.rootdir, "data", "results", context_name, "proteomics")
    merged_df = pd.DataFrame(columns=["ENTREZ_GENE_ID", "expressed", "high"])
    merged_df.set_index(["ENTREZ_GENE_ID"], inplace=True)
    merged_hi_df = merged_df
    
    for group in group_names:
        read_filepath = os.path.join(output_dir, f"bool_prot_Matrix_{context_name}_{group}.csv")
        read_df = pd.read_csv(read_filepath)
        read_df.set_index("ENTREZ_GENE_ID", inplace=True)
        read_df = read_df[['expressed', 'high']]
        
        if not merged_df.empty:
            merged_df = pd.merge(merged_df, read_df["expressed"], right_index=True, left_index=True)
            merged_hi_df = pd.merge(merged_hi_df, read_df["high"], right_index=True, left_index=True)
        
        else:
            merged_df = read_df["expressed"].to_frame()
            merged_hi_df = read_df["high"].to_frame()
    
    if len(merged_df.columns) > 1:
        merged_df.apply(lambda x: sum(x) / len(merged_df.columns) >= group_ratio, axis=1, result_type="reduce")
        merged_hi_df.apply(lambda x: sum(x) / len(merged_hi_df.columns) >= hi_group_ratio, axis=1, result_type="reduce")
    
    out_df = pd.merge(merged_df, merged_hi_df, right_index=True, left_index=True)
    out_filepath = os.path.join(output_dir, f"Proteomics_{context_name}.csv")
    out_df.to_csv(out_filepath, index_label="ENTREZ_GENE_ID")
    print("Test Data Saved to {}".format(out_filepath))


# read data from csv files
def load_proteomics_tests(filename, context_name):
    """
    Description....
    """
    
    def load_empty_dict():
        savepath = os.path.join(
            configs.rootdir,
            "data",
            "data_matrices",
            "placeholder",
            "placeholder_empty_data.csv",
        )
        dat = pd.read_csv(savepath, index_col="ENTREZ_GENE_ID")
        return "dummy", dat
    
    if (not filename or filename == "None"):  # if not using proteomics load empty dummy data matrix
        return load_empty_dict()
    
    inquiry_full_path = os.path.join(configs.rootdir, "data", "config_sheets", filename)
    if not os.path.isfile(inquiry_full_path):  # check that config file exists
        print("Error: file not found {}".format(inquiry_full_path))
        sys.exit()
    
    filename = "Proteomics_{}.csv".format(context_name)
    fullsavepath = os.path.join(
        configs.rootdir, "data", "results", context_name, "proteomics", filename
    )
    if os.path.isfile(fullsavepath):
        data = pd.read_csv(fullsavepath, index_col="ENTREZ_GENE_ID")
        print("Read from {}".format(fullsavepath))
        
        return context_name, data
    
    else:
        print(
            f"Proteomics gene expression file for {context_name} was not found at {fullsavepath}. This may be "
            f"intentional. Contexts where microarray data can be found in /work/data/results/{context_name}/ will "
            f"still be used for other contexts if found."
        )
        
        return load_empty_dict()


def proteomics_gen(TEMP):
    """
    Description
    """
    return "TEMP"


def main(argv):
    """
    Description
    """
    parser = argparse.ArgumentParser(
        prog="proteomics_gen.py",
        description="Description goes here",
        epilog="For additional help, please post questions/issues in the MADRID GitHub repo at "
               "https://github.com/HelikarLab/MADRID or email babessell@gmail.com",
    )
    parser.add_argument(
        "-c",
        "--config-file",
        type=str,
        required=True,
        dest="config_file",
        help="The configuration file for proteomics",
    )
    parser.add_argument(
        "-r",
        "--replicate-ratio",
        type=float,
        required=False,
        default=0.5,
        dest="rep_ratio",
        help="Ratio of replicates required for a gene to be considered active in that group",
    )
    parser.add_argument(
        "-b",
        "--batch-ratio",
        type=float,
        required=False,
        default=0.5,
        dest="group_ratio",
        help="Ratio of groups (batches or studies) required for a gene to be considered active in a context",
    )
    parser.add_argument(
        "-hr",
        "--high-replicate-ratio",
        type=float,
        required=False,
        default=0.5,
        dest="hi_rep_ratio",
        help="Ratio of replicates required for a gene to be considered high-confidence in that group",
    )
    parser.add_argument(
        "-hb",
        "--high-batch-ratio",
        type=float,
        required=False,
        default=0.5,
        dest="hi_group_ratio",
        help="Ratio of groups (batches or studies) required for a gene to be considered high-confidence in a context",
    )
    
    parser.add_argument(
        "-q",
        "--quantile",
        type=int,
        required=False,
        default=25,
        dest="quantile",
        help="The quantile of genes to accept. This should be an integer from 0% (no proteins pass) "
             "to 100% (all proteins pass).",
    )
    args = parser.parse_args()
    
    suppfile = args.config_file
    rep_ratio = args.rep_ratio
    group_ratio = args.group_ratio
    hi_rep_ratio = args.hi_rep_ratio
    hi_group_ratio = args.hi_group_ratio
    quantile = args.quantile / 100
    
    prot_config_filepath = os.path.join(
        configs.rootdir, "data", "config_sheets", suppfile
    )
    print('Config file is at "{}"'.format(prot_config_filepath))
    
    xl = pd.ExcelFile(prot_config_filepath)
    sheet_names = xl.sheet_names
    
    for context_name in sheet_names:
        datafilename = "".join(["protein_abundance_", context_name, ".csv"])
        config_sheet = pd.read_excel(prot_config_filepath, sheet_name=context_name)
        groups = config_sheet["Group"].unique().tolist()
        
        for group in groups:
            group_idx = np.where([True if g == group else False for g in config_sheet["Group"].tolist()])
            cols = np.take(config_sheet["SampleName"].to_numpy(), group_idx).ravel().tolist() + [
                "Gene Symbol",
                "uniprot"
            ]
            cols = np.take(config_sheet["SampleName"].to_numpy(), group_idx).ravel().tolist() + ["Gene Symbol"]
            
            proteomics_data = load_proteomics_data(datafilename, context_name)
            proteomics_data = proteomics_data.loc[:, cols]
            
            sym2id = load_gene_symbol_map(
                gene_symbols=proteomics_data["Gene Symbol"].tolist(),
                filename="proteomics_entrez_map.csv"
            )
            
            # map gene symbol to ENTREZ_GENE_ID
            proteomics_data.dropna(subset=["Gene Symbol"], inplace=True)
            
            try:
                proteomics_data.drop(columns=["uniprot"], inplace=True)
            
            except KeyError:
                pass
            
            proteomics_data = proteomics_data.groupby(["Gene Symbol"]).agg("max")
            proteomics_data["ENTREZ_GENE_ID"] = sym2id["Gene ID"]
            proteomics_data.dropna(subset=["ENTREZ_GENE_ID"], inplace=True)
            proteomics_data.set_index("ENTREZ_GENE_ID", inplace=True)
            
            # save proteomics data by test
            abundance_to_bool_group(context_name, group, proteomics_data, rep_ratio, hi_rep_ratio, quantile)
        
        to_bool_context(context_name, group_ratio, hi_group_ratio, groups)
    
    return True


if __name__ == "__main__":
    main(sys.argv[1:])
