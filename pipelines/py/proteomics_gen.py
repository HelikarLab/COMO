#!/usr/bin/python3

import argparse
import sys
import os
import time
import unidecode
import pandas as pd
import numpy as np
import instruments
from project import configs


# Load Proteomics
def load_proteomics_data(datafilename, model_name):
    """
    Add description......
    """
    dataFullPath = os.path.join(
        configs.rootdir, "data", "data_matrices", model_name, datafilename
    )
    print('Data matrix is at "{}"'.format(dataFullPath))
    if os.path.isfile(dataFullPath):
        proteomics_data = pd.read_csv(dataFullPath, header=0)
    else:
        print("Error: file not found: {}".format(dataFullPath))
        return None

    # Preprocess data, drop na, duplicate ';' in Gene names,
    proteomics_data["Gene names"] = proteomics_data["Gene names"].astype(str)
    proteomics_data.dropna(subset=["Gene names"], inplace=True)
    pluralnames = proteomics_data[
        proteomics_data["Gene names"].str.contains(";") == True
    ]
    for idx, row in pluralnames.iterrows():
        # print('{}: {}'.format(idx, row['Gene names']))
        names = row["Gene names"].split(";")
        rows = []
        for name in names:
            rowcopy = row.copy()
            rowcopy["Gene names"] = name
            rows.append(rowcopy)
        proteomics_data.drop(index=idx, inplace=True)
        proteomics_data = pd.concat(
            [proteomics_data, pd.DataFrame(rows)], ignore_index=True
        )

    proteomics_data.rename(columns={"Gene names": "Gene Symbol"}, inplace=True)
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


# determine expression logicals and save results
def save_proteomics_tests(model_name, testdata, expr_prop, top_prop, percentile):
    """
    Descrioption....
    """
    # Logical Calculation
    thresholds = testdata.quantile(1 - percentile, axis=0)
    testbool = pd.DataFrame(0, columns=list(testdata), index=testdata.index)
    for col in list(testdata):
        testbool.loc[testdata[col] > thresholds[col], [col]] = 1

    testdata["pos"] = (testdata > 0).sum(axis=1) / testdata.count(axis=1)
    testdata["expressed"] = 0
    testdata.loc[(testdata["pos"] > expr_prop), ["expressed"]] = 1
    testdata["top"] = 0
    testdata.loc[(testdata["pos"] > top_prop), ["top"]] = 1
    output_dir = os.path.join(configs.rootdir, "data", "results", model_name)
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(
        configs.rootdir,
        "data",
        "results",
        model_name,
        "Proteomics_{}.csv".format(model_name),
    )
    testdata.to_csv(output_path, index_label="ENTREZ_GENE_ID")
    print("Test Data Saved to {}".format(output_path))


# read data from csv files
def load_proteomics_tests(filename, model_name):
    """
    Description....
    """

    def load_empty_dict():
        savepath = os.path.join(
            configs.rootdir,
            "data",
            "data_matrices",
            "dummy",
            "dummy_proteomics_data.csv",
        )
        dat = pd.read_csv(savepath, index_col="ENTREZ_GENE_ID")
        return "dummy", dat

    if (
        not filename or filename == "None"
    ):  # if not using proteomics as a data source load empty dummy data matrix
        return load_empty_dict()

    inquiry_full_path = os.path.join(configs.rootdir, "data", "config_sheets", filename)
    if not os.path.isfile(
        inquiry_full_path
    ):  # check that config file exist (isn't needed to load, but helps user)
        print("Error: file not found {}".format(inquiry_full_path))
        sys.exit()

    filename = "Proteomics_{}.csv".format(model_name)
    fullsavepath = os.path.join(
        configs.rootdir, "data", "results", model_name, filename
    )
    if os.path.isfile(fullsavepath):
        data = pd.read_csv(fullsavepath, index_col="ENTREZ_GENE_ID")
        print("Read from {}".format(fullsavepath))
        return model_name, data
    else:
        print(
            f"Microarray gene expression file for {model_name} was not found at {fullsavepath}. This may be "
            f"intentional. Contexts where microarray data can be found in /work/data/results/{model_name}/ will "
            f"still be used for other contexts if found."
        )
        return load_empty_dict()


def proteomics_gen(temp):
    """
    Description
    """
    return "corn"


def main(argv):
    """
    Description
    """
    suppfile = "proteomics_data_inputs.xlsx"
    expr_prop = 0.5
    top_prop = 0.9
    percentile = 25

    # TODO: Fix the descriptions in these argument parsers
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
        "-e",
        "--expression_proportion",
        type=float,
        required=True,
        dest="expression_proportion",
        help="Ratio of replicates required for a gene to be considered active in that sampele",
    )
    parser.add_argument(
        "-t",
        "--top_proportion",
        type=float,
        required=True,
        dest="top_proportion",
        help="Genes can be considered high confidence if they are expressed in a high proportion of samples. "
             + "High confidence genes will be considered expressed regardless of agreement with other data sources",
    )
    parser.add_argument(
        "-p",
        "--percentile",
        type=float,
        required=True,
        dest="percentile",
        help="The quantile/percentile of genes to accept. This should be a number from 1 to 100",
    )
    args = parser.parse_args()

    suppfile = args.config_file
    expression_proportion = args.expression_proportion
    top_proportion = args.top_proportion
    percentile = args.percentile / 100

    # TODO: Remove this after verifying argparse works
    # try:
    #     opts, args = getopt.getopt(
    #         argv,
    #         "hc:e:t:p:",
    #         ["config_file=", "expr_proportion=", "top_proportion=", "percentile="],
    #     )
    # except getopt.GetoptError:
    #     print(
    #         "python3 proteomics_gen.py -d <config file> -s <supplementary data file> -e <expression_proportion> -t <top_proportion>"
    #     )
    #     sys.exit(2)
    # for opt, arg in opts:
    #     if opt == "-h":
    #         print(
    #             "python3 proteomics_gen.py -d <data file> -s <supplementary data file> -e <expression_proportion> -t <top_proportion>"
    #         )
    #         sys.exit()
    #     elif opt in ("-c", "--config_file"):
    #         suppfile = arg
    #     elif opt in ("-e"):
    #         expr_prop = float(arg)
    #     elif opt in ("-t"):
    #         top_prop = float(arg)
    #     elif opt in ("-p"):
    #         percentile = float(arg) / 100
    # print('Data file is "{}"'.format(datafilename))
    prot_config_filepath = os.path.join(
        configs.rootdir, "data", "config_sheets", suppfile
    )
    print('Config file is at "{}"'.format(prot_config_filepath))

    xl = pd.ExcelFile(prot_config_filepath)
    sheet_names = xl.sheet_names
    for model_name in sheet_names:
        datafilename = "".join(["ProteomicsDataMatrix_", model_name, ".csv"])
        config_sheet = pd.read_excel(
            prot_config_filepath, sheet_name=model_name, header=0
        )
        cols = config_sheet["SampleName"].to_list() + [
            "Gene Symbol",
            "Majority protein IDs",
        ]
        proteomics_data = load_proteomics_data(datafilename, model_name)
        proteomics_data = proteomics_data.loc[:, cols]
        sym2id = load_gene_symbol_map(
            gene_symbols=proteomics_data["Gene Symbol"].tolist(),
            filename="proteomics_entrez_map.csv",
        )

        # map gene symbol to ENTREZ_GENE_ID
        proteomics_data.dropna(subset=["Gene Symbol"], inplace=True)
        proteomics_data.drop(columns=["Majority protein IDs"], inplace=True)
        proteomics_data = proteomics_data.groupby(["Gene Symbol"]).agg("max")
        # proteomics_data.set_index('Gene Symbol', inplace=True)
        proteomics_data["ENTREZ_GENE_ID"] = sym2id["Gene ID"]
        proteomics_data.dropna(subset=["ENTREZ_GENE_ID"], inplace=True)
        proteomics_data.set_index("ENTREZ_GENE_ID", inplace=True)
        # prote_data_filepath = os.path.join(configs.rootdir, 'data', 'proteomics_data.csv')
        # if not os.path.isfile(prote_data_filepath):
        # proteomics_data.to_csv(prote_data_filepath, index_label='ENTREZ_GENE_ID')

        # save proteomics data by test
        save_proteomics_tests(
            model_name, proteomics_data, expr_prop, top_prop, percentile
        )
    return True


if __name__ == "__main__":
    main(sys.argv[1:])
