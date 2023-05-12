#!/usr/bin/python3
import os
import sys
import json
import argparse
import pandas as pd
from pathlib import Path
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri


import rpy2_api
import GSEpipelineFast
from project import configs
from async_bioservices import async_bioservices
from async_bioservices.input_database import InputDatabase
from async_bioservices.output_database import OutputDatabase
from async_bioservices.taxon_ids import TaxonIDs
from rpy2.robjects.packages import importr

pandas2ri.activate()

# read and translate R functions
# f = open(os.path.join(configs.rootdir, "py", "rscripts", "DGE.R"), "r")
# string = f.read()
# f.close()
# DGEio = SignatureTranslatedAnonymousPackage(string, "DGEio")
DGEio = rpy2_api.Rpy2(r_file_path=Path(configs.rootdir, "py", "rscripts", "DGE.R"))

# f = open(os.path.join(configs.rootdir, "py", "rscripts", "fitAffy.R"), "r")
# string = f.read()
# f.close()
# affyio = SignatureTranslatedAnonymousPackage(string, "affyio")
affyio = rpy2_api.Rpy2(r_file_path=Path(configs.rootdir, "py", "rscripts", "fitAffy.R"))

# f = open(os.path.join(configs.rootdir, "py", "rscripts", "fitAgilent.R"), "r")
# string = f.read()
# f.close()
# agilentio = SignatureTranslatedAnonymousPackage(string, "agilentio")
agilentio = rpy2_api.Rpy2(r_file_path=Path(configs.rootdir, "py", "rscripts", "fitAgilent.R"))

def pharse_configs(
    config_filepath: str | Path,
    sheet: str
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """
    Read microarray config files
    param: config_filepath - string, path to microarray formatted disease configuration xlsx file
    param: sheet - string, sheet name to read
    return:
    """
    xl = pd.ExcelFile(config_filepath)
    sheet_name = xl.sheet_names
    inqueries = pd.read_excel(config_filepath, sheet_name=sheet_name, header=0)
    inqueries[sheet].fillna(method="ffill", inplace=True)
    df = inqueries[sheet].loc[:, ["GSE ID", "Samples", "GPL ID", "Instrument"]]
    df_target = inqueries[sheet].loc[:, ["Samples", "Experiment"]]
    df_target.rename(
        columns={"Samples": "FileName", "Experiment": "Condition"}, inplace=True
    )
    df_target["FileName"] = df_target["FileName"].astype(str) + ".txt.gz"
    df_target["SampleNumber"] = 1 + df_target.index.values
    df_target = df_target[["SampleNumber", "FileName", "Condition"]]
    df_target["Condition"] = df_target["Condition"].str.lower()
    return df, df_target


def get_microarray_diff_gene_exp(
    config_filepath: str | Path,
    disease_name: str,
    target_path: str | Path,
    taxon_id: TaxonIDs
) -> tuple[pd.DataFrame, str]:
    """
    Get differential gene expression for microarray data

    param: config_filepath - string, path to microarray formatted disease configuration xlsx file
    param: disease_name - string, disease name which should correspond to sheet name in disease config xlsx file
    param: target_path - string, path to save disease targets

    return: dataframe with fold changes, FDR adjusted p-values,
    """
    query_table, df_target = pharse_configs(config_filepath, disease_name)
    df_target.to_csv(target_path, index=False, sep="\t")
    print(query_table)
    sr = query_table["GSE ID"]
    gse_ids = sr[sr.str.match("GSE")].unique()
    gse_id = gse_ids[0]

    inst = query_table["Instrument"]
    inst_name = inst.unique()
    inst_name = inst_name[0]
    print(inst_name)

    gseXXX = GSEpipelineFast.GSEproject(gse_id, query_table, configs.rootdir)
    for key, val in gseXXX.platforms.items():
        raw_dir = os.path.join(gseXXX.gene_dir, key)
        print(f"{key}:{val}, {raw_dir}")
        if inst_name == "affy":
            affy_function = rpy2_api.Rpy2(Path(configs.rootdir, "py", "rscripts", "DGE.R"), raw_dir, target_path)
            affy_function.call_function("fitaffydir")
            diff_exp_df = affyio.call_function("fitaffydir")
        elif inst_name == "agilent":
            agilent_function = rpy2_api.Rpy2(Path(configs.rootdir, "py", "rscripts", "fitAgilent.R"), raw_dir, target_path)
            diff_exp_df = agilent_function.call_function("fitagilent")
        diff_exp_df = ro.conversion.rpy2py(diff_exp_df)
        diff_exp_df.reset_index(inplace=True)
        diff_exp_df = diff_exp_df.rename(columns={'index': 'Affy ID'})
        print(diff_exp_df)
        print(type(diff_exp_df))
        data_top = diff_exp_df.head()
        print(data_top)
        diff_exp_df.rename(columns={"adj.P.Val": "FDR"}, inplace=True)
        print(diff_exp_df)

        input_db: InputDatabase
        if inst_name == "affy":
            input_db = InputDatabase.AFFY_ID
        elif inst_name == "agilent":
            input_db = InputDatabase.AGILENT_ID

        bdnet = async_bioservices.fetch_gene_info(
            input_values=list(map(str, diff_exp_df["Affy ID"].tolist())),
            input_db=input_db,
            output_db=[OutputDatabase.ENSEMBL_GENE_ID, OutputDatabase.GENE_ID, OutputDatabase.GENE_SYMBOL],
            taxon_id=taxon_id
        )

        diff_exp_df["Entrez"] = bdnet["Gene ID"].tolist()
        diff_exp_df["Ensembl"] = bdnet["Ensembl Gene ID"].tolist()
        diff_exp_df["Symbol"] = bdnet["Gene Symbol"].tolist()
        print(diff_exp_df)
        diff_exp_df.rename(columns={"Affy ID": "Affy"}, inplace=True)
    return diff_exp_df, gse_id


def get_rnaseq_diff_gene_exp(
    config_filepath: str | Path,
    disease_name: str,
    context_name: str,
    taxon_id: TaxonIDs
) -> tuple[pd.DataFrame, str]:
    """
    Get differential gene expression for RNA-seq data

    param: config_filepath - string, path to microarray formatted disease configuration xlsx file
    param: disease_name - string, disease name which should correspond to sheet name in disease config xlsx file
    param: context_name - string, context name which should correspond to folder in 'results' folder

    return: dataframe with fold changes, FDR adjusted p-values,
    """
    count_matrix_filename = "".join(["gene_counts_matrix_", disease_name, "_", context_name, ".csv"])
    count_matrix_path = os.path.join(
        configs.datadir,
        "data_matrices",
        context_name,
        "disease",
        count_matrix_filename
    )

    if os.path.exists(count_matrix_path):
        print("Count Matrix File is at ", count_matrix_path)
    else:
        print(f"No count matrix found at {count_matrix_path}. Please make sure file is in the correct location "
              f"with the correct name.")
        sys.exit()

    dge_call = rpy2_api.Rpy2(Path(configs.rootdir, "py", "rscripts", "DGE.R"), count_matrix_path, config_filepath, context_name, disease_name)
    dge_results = dge_call.call_function("DGE_main")
    diff_exp_df: pd.DataFrame = ro.conversion.rpy2py(dge_results)
    gse_id = "rnaseq"

    bdnet = async_bioservices.fetch_gene_info(
        input_values=list(map(str, diff_exp_df["Ensembl"].tolist())),
        input_db=InputDatabase.ENSEMBL_GENE_ID,
        output_db=[OutputDatabase.GENE_ID, OutputDatabase.AFFY_ID, OutputDatabase.GENE_SYMBOL],
        taxon_id=taxon_id
    )

    diff_exp_df["Affy"] = bdnet["Affy ID"].tolist()
    diff_exp_df["Entrez"] = bdnet["Gene ID"].tolist()
    diff_exp_df["Symbol"] = bdnet["Gene Symbol"].tolist()

    return diff_exp_df, gse_id


def write_outputs(
    diff_exp_df: pd.DataFrame,
    gse_id: str,
    context_name: str,
    disease_name: str,
    data_source: str,
    target_path: str
) -> None:
    search_col = "Ensembl" if data_source == "RNASEQ" else "Affy"
    diff_exp_df["logFC"].astype(float)
    diff_exp_df["abs_logFC"] = diff_exp_df["logFC"].abs()
    diff_exp_df["FDR"].astype(float)
    diff_exp_df.sort_values(by="abs_logFC", ascending=False, inplace=True)
    regulated = diff_exp_df[diff_exp_df["FDR"] < 0.05]
    down_regulated = regulated[regulated["logFC"] < 0]
    up_regulated = regulated[regulated["logFC"] > 0]
    diff_exp_df["regulated"] = [
        "unchanged" if gene not in regulated[search_col].tolist()
        else ("upregulated" if gene in up_regulated[search_col].tolist()
              else "downregulated")
        for gene in diff_exp_df[search_col].tolist()
    ]
    up_file = os.path.join(
        configs.rootdir,
        "data",
        "results",
        context_name,
        disease_name,
        f"Disease_UP_{gse_id}.txt"
    )
    os.makedirs(os.path.dirname(up_file), exist_ok=True)

    down_file = os.path.join(
        configs.rootdir,
        "data",
        "results",
        context_name,
        disease_name,
        f"Disease_DOWN_{gse_id}.txt"
    )
    os.makedirs(os.path.dirname(down_file), exist_ok=True)

    up_regulated = up_regulated[up_regulated["Entrez"] != "-"]
    down_regulated = down_regulated[down_regulated["Entrez"] != "-"]

    up_regulated["Entrez"].to_csv(up_file, index=False)
    down_regulated["Entrez"].to_csv(down_file, index=False)
    print(f"Upregulated genes saved to '{up_file}'")
    print(f"Downregulated genes saved to '{down_file}'")

    raw_file = os.path.join(
        configs.rootdir,
        "data",
        "results",
        context_name,
        disease_name,
        f"Raw_Fit_{gse_id}.csv",
    )
    diff_exp_df.drop(columns=["Affy"], inplace=True)  # drop for now bc commas mess up csv parsing, maybe fix later
    diff_exp_df.to_csv(raw_file, index=False)
    print(f"Raw Data saved to '{raw_file}'")

    files_dict = {
        "gse": gse_id,
        "up_regulated": up_file,
        "down_regulated": down_file,
        "raw_data": raw_file,
    }

    files_json = os.path.join(
        configs.rootdir,
        "data",
        "results",
        context_name,
        disease_name,
        "step2_results_files.json",
    )
    os.makedirs(os.path.dirname(files_json), exist_ok=True)

    with open(files_json, "w") as fp:
        json.dump(files_dict, fp)

    if data_source == "MICROARRAY":
        os.remove(target_path)


def main(argv: list[str]) -> None:
    targetfile = "targets.txt"

    parser = argparse.ArgumentParser(
        prog="disease_analysis.py",
        description="Performs differential gene expression analysis to find up and downregulated genes associated "
                    "with a disease. Significant genes are ones that have an FDR adjusted P-value < 0.05 and an "
                    "absolute fold-change greater than the threshold specified, default is 2",
        epilog="For additional help, please post questions/issues in the MADRID GitHub repo at: "
               "https://github.com/HelikarLab/MADRID or email babessell@gmail.com"
    )
    parser.add_argument(
        "-c",
        "--config-file",
        type=str,
        required=True,
        dest="config_file",
        help="The path to the configuration file",
    )
    parser.add_argument(
        "-t",
        "--context-name",
        type=str,
        required=True,
        dest="context_name",
        help="The type of context being used",
    )
    parser.add_argument(
        "-s",
        "--data-source",
        type=str,
        required=True,
        dest="data_source",
        help="Source of data being used, either rnaseq or microarray"
    )
    parser.add_argument(
        "-i",
        "--taxon-id",
        required=False,
        default=9606,
        dest="taxon_id",
        help="BioDbNet taxon ID number, also accepts 'human', or 'mouse'",
    )

    args = parser.parse_args()
    context_name = args.context_name
    config_file = args.config_file
    data_source = args.data_source.upper()
    taxon_id = args.taxon_id

    if data_source == "RNASEQ":
        config_filepath = os.path.join(configs.datadir, "config_sheets", "disease", config_file)
    elif data_source == "MICROARRAY":
        config_filepath = os.path.join(configs.datadir, "config_sheets", "disease", config_file)
    else:
        print(f"{data_source} is not a valid data source, must be either MICROARRAY or RNASEQ.")
        sys.exit()

    try:
        xl = pd.ExcelFile(config_filepath)
    except FileNotFoundError:
        print(f"Config File not found at {config_filepath}")
        sys.exit()
    except ValueError:
        print("Config file must be in xlsx format!")
        sys.exit()
    print("Config file is at ", config_filepath)

    # handle species alternative ids
    set_taxonid: TaxonIDs
    found_valid: bool = True
    if str(taxon_id).isdigit():  # Check if taxon_id is an integer
        if int(taxon_id) == 9606:
            set_taxonid = TaxonIDs.HOMO_SAPIENS
        elif int(taxon_id) == 10090:
            set_taxonid = TaxonIDs.MUS_MUSCULUS
        else:
            found_valid = False
    else:  # Taxon id is a string
        if taxon_id.upper() in ["HUMAN", "HOMO SAPIENS"]:
            set_taxonid = TaxonIDs.HOMO_SAPIENS
        elif taxon_id.upper() in ["MOUSE", "MUS MUSCULUS"]:
            set_taxonid = TaxonIDs.MUS_MUSCULUS
        else:
            found_valid = False
    
    if not found_valid:
        raise ValueError(f'--taxon-id must be either an integer, or accepted string ("mouse", "human"). Received: {taxon_id}')

    sheet_names = xl.sheet_names
    for disease_name in sheet_names:
        target_path = os.path.join(configs.rootdir, "data", targetfile)
        if data_source == "MICROARRAY":
            diff_exp_df, gse_id = get_microarray_diff_gene_exp(config_filepath, disease_name, target_path, set_taxonid)
        elif data_source == "RNASEQ":
            diff_exp_df, gse_id = get_rnaseq_diff_gene_exp(config_filepath, disease_name, context_name, set_taxonid)
        else:
            print("data_source should be either 'microarray' or 'rnaseq'")
            print("Refer to example config file for either type for formatting")
            sys.exit(2)

        write_outputs(diff_exp_df, gse_id, context_name, disease_name, data_source, target_path)


if __name__ == "__main__":
    main(sys.argv[1:])
