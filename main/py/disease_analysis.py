#!/usr/bin/python3
import argparse
import sys
import json
import os
import pandas as pd
import numpy as np
from project import configs
import GSEpipelineFast
import rpy2_api
from pathlib import Path
import GSEpipelineFast

from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr
import rpy2.robjects as ro

from async_bioservices import database_convert
from async_bioservices.input_database import InputDatabase
from async_bioservices.output_database import OutputDatabase
from rpy2.robjects.packages import importr
from pathlib import Path
import os
import pandas as pd

pandas2ri.activate()

# import R libraries
DESeq2 = importr("DESeq2")
edgeR = importr("edgeR")
readxl = importr("readxl")

# read and translate R functions
# f = open(os.path.join(configs.rootdir, "py", "rscripts", "DGE.R"), "r")
# string = f.read()
# f.close()
# DGEio = SignatureTranslatedAnonymousPackage(string, "DGEio")
DGEio_path = Path(configs.rootdir, "py", "rscripts", "DGE.R")
# DGEio = rpy2_api.Rpy2(r_file_path=Path(configs.rootdir, "py", "rscripts", "DGE.R"))
DGEio = rpy2_api.Rpy2(r_file_path=Path(configs.rootdir, "py", "rscripts", "DGE.R"))

# f = open(os.path.join(configs.rootdir, "py", "rscripts", "fitAffy.R"), "r")
# string = f.read()
# f.close()
# affyio = SignatureTranslatedAnonymousPackage(string, "affyio")
affyio_path = Path(configs.rootdir, "py", "rscripts", "fitAffy.R")
# affyio = rpy2_api.Rpy2(r_file_path=Path(configs.rootdir, "py", "rscripts", "fitAffy.R"))
affyio = rpy2_api.Rpy2(r_file_path=Path(configs.rootdir, "py", "rscripts", "fitAffy.R"))

# f = open(os.path.join(configs.rootdir, "py", "rscripts", "fitAgilent.R"), "r")
# string = f.read()
# f.close()
# agilentio = SignatureTranslatedAnonymousPackage(string, "agilentio")
agilentio_path = Path(configs.rootdir, "py", "rscripts", "fitAgilent.R")
# agilentio = rpy2_api.Rpy2(r_file_path=Path(configs.rootdir, "py", "rscripts", "fitAgilent.R"))
agilentio = rpy2_api.Rpy2(r_file_path=Path(configs.rootdir, "py", "rscripts", "fitAgilent.R"))


def breakDownEntrezs(disease):
    """
    Split multi-part entrez ids into multiple rows in the dataframe
    param: disease - df with columns: Gene ID
    return: df with columns:
    """
    disease["Gene ID"] = disease["Gene ID"].str.replace("///", "//")
    single_gene_names = disease[~disease["Gene ID"].str.contains("//")].reset_index(
        drop=True
    )
    multiple_gene_names = disease[
        disease["Gene ID"].str.contains("//")
    ].reset_index(drop=True)
    breaks_gene_names = pd.DataFrame(columns=["Gene ID"])

    for index, row in multiple_gene_names.iterrows():
        for genename in row["Gene ID"].split("//"):
            breaks_gene_names = breaks_gene_names.append(
                {"Gene ID": genename}, ignore_index=True
            )
    gene_expressions = single_gene_names.append(breaks_gene_names, ignore_index=True)

    return gene_expressions


# def get_entrez_id(regulated, output_full_path, in_db, taxon_id, full_flag=False):
#     """
#     Fetch entrez ids using bioDbNet
#     param: regulated -  df, with columns:
#     param: output_full_path - path to
#     param: in_db - bioDBnet input database query name
#     param: taxon_id - bioDBnet taxon id
#     param: full_flag - boolean, True if breaking down multi- entrez ids into seperate rows
#     return: df,
#     """
#     disease = database_convert.fetch_gene_info(
#         input_values=list(regulated.index.values),
#         input_db=in_db,
#         output_db=["Gene Symbol"],
#         taxon_id=taxon_id,
#     )
#     # disease = fetch_gene_info(list(regulated.index.values),
#     #                           input_db=in_db,
#     #                           output_db=["Gene Symbol"],
#     #                           delay=15,
#     #                           taxon_id=taxon_id
#     #                           )
#     disease.reset_index(inplace=True)
#     # disease.drop(columns=["Ensembl Gene ID"], inplace=True)
#     disease.replace(to_replace="-", value=np.nan, inplace=True)
#     if not full_flag:
#         disease.dropna(how="any", subset=disease.columns[[0]], inplace=True)
#         gene_expressions = breakDownEntrezs(disease)
#         gene_expressions.set_index(gene_expressions.columns[0], inplace=True)
#     else:
#         gene_expressions = disease
#
#     # gene_expressions["Gene ID"].to_csv(outputFullPath, index=False)
#     gene_expressions[gene_expressions.columns[0]].to_csv(output_full_path, index=False)
#     return gene_expressions


def pharse_configs(config_filepath, sheet):
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


def get_microarray_diff_gene_exp(config_filepath, disease_name, target_path, taxon_id):
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
            affy_function = affyio.call_function("fitaffydir")
            diff_exp_df = affy_function(raw_dir, target_path)
            # diff_exp_df = affyio.fitaffydir(raw_dir, target_path)
        elif inst_name == "agilent":
            agilent_function = rpy2_api.Rpy2(agilentio_path, raw_dir, target_path)
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

        if inst_name == "affy":
            input_db: InputDatabase = InputDatabase.AFFY_ID
        elif inst_name == "agilent":
            input_db: InputDatabase = InputDatabase.AGILENT_ID

        bdnet = database_convert.fetch_gene_info(
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


def get_rnaseq_diff_gene_exp(config_filepath, disease_name, context_name, taxon_id):
    """
    Get differential gene expression for RNA-seq data

    param: config_filepath - string, path to microarray formatted disease configuration xlsx file
    param: disease_name - string, disease name which should correspond to sheet name in disease config xlsx file
    param: context_name - string, context name which should correspond to folder in 'results' folder

    return: dataframe with fold changes, FDR adjusted p-values,
    """
    count_matrix_filename = "".join(["gene_counts_matrix_", disease_name, "_", context_name, ".csv"])
    count_matrix_path = os.path.join(configs.datadir,
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

    DGEio_function = rpy2_api.Rpy2(DGEio_path, count_matrix_path, config_filepath, context_name, disease_name)
    diff_exp_df = DGEio_function.call_function("DGE_main")
    # diff_exp_df = DGEio.DGE_main(count_matrix_path, config_filepath, context_name, disease_name)
    diff_exp_df = ro.conversion.rpy2py(diff_exp_df)
    gse_id = "rnaseq"

    bdnet = database_convert.fetch_gene_info(
        input_values=list(map(str, diff_exp_df["Ensembl"].tolist())),
        input_db=InputDatabase.ENSEMBL_GENE_ID,
        output_db=[OutputDatabase.GENE_ID, OutputDatabase.AFFY_ID, OutputDatabase.GENE_SYMBOL],
        taxon_id=taxon_id
    )

    diff_exp_df["Affy"] = bdnet["Affy ID"].tolist()
    diff_exp_df["Entrez"] = bdnet["Gene ID"].tolist()
    diff_exp_df["Symbol"] = bdnet["Gene Symbol"].tolist()

    return diff_exp_df, gse_id


def write_outputs(diff_exp_df, gse_id, context_name, disease_name, data_source, target_path):
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
        "GSE": gse_id,
        "UP_Reg": up_file,
        "DN_Reg": down_file,
        "RAW_Data": raw_file,
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


def main(argv):
    targetfile = "targets.txt"
    count_matrix = None

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
    if type(taxon_id) == str:
        if taxon_id.upper() == "HUMAN" or taxon_id.upper() == "HOMO SAPIENS":
            taxon_id = 9606
        elif taxon_id.upper() == "MOUSE" or taxon_id.upper() == "MUS MUSCULUS":
            taxon_id = 10090
        else:
            print('--taxon-id must be either an integer, or accepted string ("mouse", "human")')
            sys.exit()
    elif type(taxon_id) != int:
        print('--taxon-id must be either an integer, or accepted string ("mouse", "human")')
        sys.exit()

    sheet_names = xl.sheet_names
    for disease_name in sheet_names:
        target_path = os.path.join(configs.rootdir, "data", targetfile)
        if data_source == "MICROARRAY":
            diff_exp_df, gse_id = get_microarray_diff_gene_exp(config_filepath, disease_name, target_path, taxon_id)
        elif data_source == "RNASEQ":
            diff_exp_df, gse_id = get_rnaseq_diff_gene_exp(config_filepath, disease_name, context_name, taxon_id)
        else:
            print("data_source should be either 'microarray' or 'rnaseq'")
            print("Refer to example config file for either type for formatting")
            sys.exit(2)

        write_outputs(diff_exp_df, gse_id, context_name, disease_name, data_source, target_path)


if __name__ == "__main__":
    main(sys.argv[1:])
