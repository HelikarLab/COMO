# ruff: noqa


import argparse
import json
from pathlib import Path

import pandas as pd
import rpy2.robjects as ro
import rpy2_api
from fast_bioservices import BioDBNet, Input, Output
from project import Config
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import importr

configs = Config()

pandas2ri.activate()

# import R libraries
DESeq2 = importr("DESeq2")
edgeR = importr("edgeR")
readxl = importr("readxl")

DGEio = rpy2_api.Rpy2(r_file_path=Path(configs.code_dir, "rscripts", "DGE.R"))


def get_rnaseq_diff_gene_exp(config_filepath, disease_name, context_name, taxon_id):
    """
    Get differential gene expression for RNA-seq data

    param: config_filepath - path to rna-seq disease configuration xlsx file
    param: disease_name - string, disease name which should correspond to sheet name in disease config xlsx file
    param: context_name - string, context name which should correspond to folder in 'results' folder

    return: dataframe with fold changes, FDR adjusted p-values,
    """
    count_matrix_filename = "".join(["gene_counts_matrix_", disease_name, "_", context_name, ".csv"])
    count_matrix_path: Path = configs.data_dir / "data_matrices" / context_name / "disease" / count_matrix_filename

    if count_matrix_path.exists():
        print("Count Matrix File is at ", count_matrix_path)
    else:
        raise FileNotFoundError(f"Count Matrix File not found at {count_matrix_path}")

    diff_exp_df = DGEio.call_function("DGE_main", count_matrix_path, config_filepath, context_name, disease_name)
    diff_exp_df = ro.conversion.rpy2py(diff_exp_df)
    gse_id = "rnaseq"

    biodbnet = BioDBNet()
    bdnet = biodbnet.db2db(
        input_values=list(map(str, diff_exp_df["Ensembl"].tolist())),
        input_db=Input.ENSEMBL_GENE_ID,
        output_db=[Output.GENE_ID, Output.AFFY_ID, Output.GENE_SYMBOL],
        taxon=taxon_id,
    )

    diff_exp_df["Affy"] = bdnet["Affy ID"].tolist()
    diff_exp_df["Entrez"] = bdnet["Gene ID"].tolist()
    diff_exp_df["Symbol"] = bdnet["Gene Symbol"].tolist()

    return diff_exp_df, gse_id


def write_outputs(diff_exp_df, gse_id, context_name, disease_name, target_path):
    search_col = "Ensembl"
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
        else ("upregulated" if gene in up_regulated[search_col].tolist() else "downregulated")
        for gene in diff_exp_df[search_col].tolist()
    ]
    up_file = configs.result_dir / context_name / disease_name / f"Disease_UP_{gse_id}.txt"
    down_file = configs.result_dir / context_name / disease_name / f"Disease_DOWN_{gse_id}.txt"

    up_file.parent.mkdir(parents=True, exist_ok=True)
    down_file.parent.mkdir(parents=True, exist_ok=True)

    up_regulated = up_regulated[up_regulated["Entrez"] != "-"]
    down_regulated = down_regulated[down_regulated["Entrez"] != "-"]

    up_regulated["Entrez"].to_csv(up_file, index=False)
    down_regulated["Entrez"].to_csv(down_file, index=False)
    print(f"Upregulated genes saved to '{up_file}'")
    print(f"Downregulated genes saved to '{down_file}'")

    raw_file = configs.result_dir / context_name / disease_name / f"Raw_Fit_{gse_id}.csv"
    diff_exp_df.drop(columns=["Affy"], inplace=True)  # drop for now bc commas mess up csv parsing, maybe fix later
    diff_exp_df.to_csv(raw_file, index=False)
    print(f"Raw Data saved to '{raw_file}'")

    files_dict = {
        "gse": gse_id,
        "up_regulated": up_file,
        "down_regulated": down_file,
        "raw_data": raw_file,
    }

    files_json = configs.result_dir / context_name / disease_name / "step2_results_files.json"
    files_json.parent.mkdir(parents=True, exist_ok=True)
    with open(files_json, "w") as fp:
        json.dump(files_dict, fp)


def main():
    target_file = "targets.txt"

    parser = argparse.ArgumentParser(
        prog="disease_analysis.py",
        description="Performs differential gene expression analysis to find up and downregulated genes associated "
        "with a disease. Significant genes are ones that have an FDR adjusted P-value < 0.05 and an "
        "absolute fold-change greater than the threshold specified, default is 2",
        epilog="For additional help, please post questions/issues in the MADRID GitHub repo at: "
        "https://github.com/HelikarLab/MADRID or email babessell@gmail.com",
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
    taxon_id = args.taxon_id
    config_filepath = configs.config_dir / "disease" / config_file

    if not config_filepath.exists():
        raise FileNotFoundError(f"Config file not found at {config_filepath}")
    if not config_filepath.suffix == ".xlsx":
        raise ValueError("Config file must be in xlsx format!")
    print("Config file is at ", config_filepath)
    xl = pd.ExcelFile(config_filepath)

    # handle species alternative ids
    if isinstance(taxon_id, str):
        if taxon_id.upper() == "HUMAN" or taxon_id.upper() == "HOMO SAPIENS":
            taxon_id = 9606
        elif taxon_id.upper() == "MOUSE" or taxon_id.upper() == "MUS MUSCULUS":
            taxon_id = 10090
        else:
            raise ValueError("taxon_id must be either an integer, or accepted string ('mouse', 'human')")
    elif not isinstance(taxon_id, int):
        raise ValueError("taxon_id must be either an integer, or accepted string ('mouse', 'human')")

    sheet_names = xl.sheet_names
    for disease_name in sheet_names:
        target_path = configs.data_dir / target_file
        diff_exp_df, gse_id = get_rnaseq_diff_gene_exp(config_filepath, disease_name, context_name, taxon_id)
        write_outputs(diff_exp_df, gse_id, context_name, disease_name, target_path)


if __name__ == "__main__":
    main()
