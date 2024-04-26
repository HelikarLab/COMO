#!/usr/bin/python3
import argparse
import json
import os
import sys
from pathlib import Path

import GSEpipelineFast
import pandas as pd
import rpy2.robjects as ro
import rpy2_api
from fast_bioservices.biodbnet import BioDBNet, Input, Output
from arguments import config_file_arg, context_names_arg, data_source_arg, taxon_id_arg
from project import Configs
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

configs = Configs()

pandas2ri.activate()

# import R libraries
DESeq2 = importr("DESeq2")
edgeR = importr("edgeR")
readxl = importr("readxl")

# read and translate R functions
# f = open(os.path.join(configs.rootdir, "src", "rscripts", "DGE.R"), "r")
# string = f.read()
# f.close()
# DGEio = SignatureTranslatedAnonymousPackage(string, "DGEio")
DGEio = rpy2_api.Rpy2(r_file_path=Path(configs.root_dir, "src", "rscripts", "DGE.R"))

# f = open(os.path.join(configs.rootdir, "src", "rscripts", "fitAffy.R"), "r")
# string = f.read()
# f.close()
# affyio = SignatureTranslatedAnonymousPackage(string, "affyio")
affyio = rpy2_api.Rpy2(
    r_file_path=Path(configs.root_dir, "src", "rscripts", "fitAffy.R")
)

# f = open(os.path.join(configs.rootdir, "src", "rscripts", "fitAgilent.R"), "r")
# string = f.read()
# f.close()
# agilentio = SignatureTranslatedAnonymousPackage(string, "agilentio")
agilentio = rpy2_api.Rpy2(
    r_file_path=Path(configs.root_dir, "src", "rscripts", "fitAgilent.R")
)


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
    multiple_gene_names = disease[disease["Gene ID"].str.contains("//")].reset_index(
        drop=True
    )
    breaks_gene_names = pd.DataFrame(columns=["Gene ID"])

    for index, row in multiple_gene_names.iterrows():
        for genename in row["Gene ID"].split("//"):
            breaks_gene_names = breaks_gene_names.append(  # type: ignore
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
    inqueries[sheet].fillna(method="ffill", inplace=True)  # type: ignore
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
    config_filepath, disease_name, target_path, taxon_id, biodbnet: BioDBNet
):
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

    gseXXX = GSEpipelineFast.GSEproject(gse_id, query_table, configs.root_dir)
    for key, val in gseXXX.platforms.items():
        raw_dir = os.path.join(gseXXX.gene_dir, key)
        print(f"{key}:{val}, {raw_dir}")
        if inst_name == "affy":
            affy_function = affyio.call_function("fitaffydir")
            diff_exp_df = affy_function(raw_dir, target_path)
            # diff_exp_df = affyio.fitaffydir(raw_dir, target_path)
        elif inst_name == "agilent":
            diff_exp_df = agilentio.call_function("fitagilent", raw_dir, target_path)
            # diff_exp_df = agilentio.fitagilent(raw_dir, target_path)
        diff_exp_df = ro.conversion.rpy2py(diff_exp_df)  # type: ignore
        diff_exp_df.reset_index(inplace=True)
        diff_exp_df = diff_exp_df.rename(columns={"index": "Affy ID"})
        print(diff_exp_df)
        print(type(diff_exp_df))
        data_top = diff_exp_df.head()
        print(data_top)
        diff_exp_df.rename(columns={"adj.P.Val": "FDR"}, inplace=True)
        print(diff_exp_df)

        if inst_name == "affy":
            input_db: Input = Input.AFFY_ID
        elif inst_name == "agilent":
            input_db: Input = Input.AGILENT_ID

        conversion = biodbnet.db2db(
            input_values=list(map(str, diff_exp_df["Affy ID"].tolist())),
            input_db=input_db,  # type: ignore
            output_db=[
                Output.ENSEMBL_GENE_ID,
                Output.GENE_ID,
                Output.GENE_SYMBOL,
            ],
            taxon=taxon_id,
        )

        diff_exp_df["Entrez"] = conversion["Gene ID"].tolist()
        diff_exp_df["Ensembl"] = conversion["Ensembl Gene ID"].tolist()
        diff_exp_df["Symbol"] = conversion["Gene Symbol"].tolist()
        print(diff_exp_df)
        diff_exp_df.rename(columns={"Affy ID": "Affy"}, inplace=True)
    return diff_exp_df, gse_id  # type: ignore


def get_rnaseq_diff_gene_exp(
    config_filepath,
    disease_name,
    context_name,
    taxon_id,
    biodbnet: BioDBNet,
):
    """
    Get differential gene expression for RNA-seq data

    param: config_filepath - string, path to microarray formatted disease configuration xlsx file
    param: disease_name - string, disease name which should correspond to sheet name in disease config xlsx file
    param: context_name - string, context name which should correspond to folder in 'results' folder

    return: dataframe with fold changes, FDR adjusted p-values,
    """
    count_matrix_filename = "".join(
        ["gene_counts_matrix_", disease_name, "_", context_name, ".csv"]
    )
    count_matrix_path = os.path.join(
        configs.data_dir,
        "data_matrices",
        context_name,
        "disease",
        count_matrix_filename,
    )

    if os.path.exists(count_matrix_path):
        print("Count Matrix File is at ", count_matrix_path)
    else:
        print(
            f"No count matrix found at {count_matrix_path}. Please make sure file is in the correct location "
            f"with the correct name."
        )
        sys.exit()

    diff_exp_df = DGEio.call_function(
        "DGE_main", count_matrix_path, config_filepath, context_name, disease_name
    )
    diff_exp_df = ro.conversion.rpy2py(diff_exp_df)
    gse_id = "rnaseq"

    conversion = biodbnet.db2db(
        input_values=list(map(str, diff_exp_df["Ensembl"].tolist())),
        input_db=Input.ENSEMBL_GENE_ID,
        output_db=[
            Output.GENE_ID,
            Output.AFFY_ID,
            Output.GENE_SYMBOL,
        ],
        taxon=taxon_id,
    )

    diff_exp_df["Affy"] = conversion["Affy ID"].tolist()
    diff_exp_df["Entrez"] = conversion["Gene ID"].tolist()
    diff_exp_df["Symbol"] = conversion["Gene Symbol"].tolist()

    return diff_exp_df, gse_id


def write_outputs(
    diff_exp_df, gse_id, context_name, disease_name, data_source, target_path
):
    search_col = "Ensembl" if data_source == "RNASEQ" else "Affy"
    diff_exp_df["logFC"].astype(float)
    diff_exp_df["abs_logFC"] = diff_exp_df["logFC"].abs()
    diff_exp_df["FDR"].astype(float)
    diff_exp_df.sort_values(by="abs_logFC", ascending=False, inplace=True)
    regulated = diff_exp_df[diff_exp_df["FDR"] < 0.05]
    down_regulated = regulated[regulated["logFC"] < 0]
    up_regulated = regulated[regulated["logFC"] > 0]
    diff_exp_df["regulated"] = [
        "unchanged"
        if gene not in regulated[search_col].tolist()
        else (
            "upregulated"
            if gene in up_regulated[search_col].tolist()
            else "downregulated"
        )
        for gene in diff_exp_df[search_col].tolist()
    ]
    up_file = os.path.join(
        configs.root_dir,
        "data",
        "results",
        context_name,
        disease_name,
        f"Disease_UP_{gse_id}.txt",
    )
    os.makedirs(os.path.dirname(up_file), exist_ok=True)

    down_file = os.path.join(
        configs.root_dir,
        "data",
        "results",
        context_name,
        disease_name,
        f"Disease_DOWN_{gse_id}.txt",
    )
    os.makedirs(os.path.dirname(down_file), exist_ok=True)

    up_regulated = up_regulated[up_regulated["Entrez"] != "-"]
    down_regulated = down_regulated[down_regulated["Entrez"] != "-"]

    up_regulated["Entrez"].to_csv(up_file, index=False)
    down_regulated["Entrez"].to_csv(down_file, index=False)
    print(f"Upregulated genes saved to '{up_file}'")
    print(f"Downregulated genes saved to '{down_file}'")

    raw_file = os.path.join(
        configs.root_dir,
        "data",
        "results",
        context_name,
        disease_name,
        f"Raw_Fit_{gse_id}.csv",
    )
    diff_exp_df.drop(
        columns=["Affy"], inplace=True
    )  # drop for now bc commas mess up csv parsing, maybe fix later
    diff_exp_df.to_csv(raw_file, index=False)
    print(f"Raw Data saved to '{raw_file}'")

    files_dict = {
        "gse": gse_id,
        "up_regulated": up_file,
        "down_regulated": down_file,
        "raw_data": raw_file,
    }

    files_json = os.path.join(
        configs.root_dir,
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
        "https://github.com/HelikarLab/MADRID or email babessell@gmail.com",
    )
    parser.add_argument(**config_file_arg)
    parser.add_argument(**context_names_arg)
    parser.add_argument(**data_source_arg)
    parser.add_argument(**taxon_id_arg)

    args = parser.parse_args()
    context_name = args.context_name
    config_file = args.config_file
    data_source = args.data_source.upper()
    taxon_id = args.taxon_id
    show_biodbnet_progress = args.show_biodbnet_progress
    use_biodbnet_cache = args.use_biodbnet_cache

    biodbnet = BioDBNet(show_progress=show_biodbnet_progress, cache=use_biodbnet_cache)

    if data_source == "RNASEQ":
        config_filepath = os.path.join(
            configs.data_dir, "config_sheets", "disease", config_file
        )
    elif data_source == "MICROARRAY":
        config_filepath = os.path.join(
            configs.data_dir, "config_sheets", "disease", config_file
        )
    else:
        print(
            f"{data_source} is not a valid data source, must be either MICROARRAY or RNASEQ."
        )
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
    if isinstance(taxon_id, str):
        if taxon_id.upper() == "HUMAN" or taxon_id.upper() == "HOMO SAPIENS":
            taxon_id = 9606
        elif taxon_id.upper() == "MOUSE" or taxon_id.upper() == "MUS MUSCULUS":
            taxon_id = 10090
        else:
            print(
                '--taxon-id must be either an integer, or accepted string ("mouse", "human")'
            )
            sys.exit()
    elif not isinstance(taxon_id, int):
        print(
            '--taxon-id must be either an integer, or accepted string ("mouse", "human")'
        )
        sys.exit()

    sheet_names = xl.sheet_names
    for disease_name in sheet_names:
        target_path = os.path.join(configs.root_dir, "data", targetfile)
        if data_source == "MICROARRAY":
            diff_exp_df, gse_id = get_microarray_diff_gene_exp(
                config_filepath,
                disease_name,
                target_path,
                taxon_id,
                biodbnet=biodbnet,
            )
        elif data_source == "RNASEQ":
            diff_exp_df, gse_id = get_rnaseq_diff_gene_exp(
                config_filepath,
                disease_name,
                context_name,
                taxon_id,
                biodbnet=biodbnet,
            )
        else:
            print("data_source should be either 'microarray' or 'rnaseq'")
            print("Refer to example config file for either type for formatting")
            sys.exit(2)

        write_outputs(
            diff_exp_df, gse_id, context_name, disease_name, data_source, target_path
        )


if __name__ == "__main__":
    main(sys.argv[1:])
