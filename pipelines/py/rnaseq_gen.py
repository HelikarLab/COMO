#!/usr/bin/python3
import os, time, sys
import pandas as pd
from rpy2.robjects.packages import importr
from rpy2.robjects import r, pandas2ri
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage
from project import configs
import argparse
import re

# enable r to py conversion
pandas2ri.activate()

# import R libraries
limma = importr("limma")
tidyverse = importr("tidyverse")
edgeR = importr("edgeR")
genefilter = importr("genefilter")
biomaRt = importr("biomaRt")
sjmisc = importr("sjmisc")
openxlsx = importr("openxlsx")
zfpkm = importr("zFPKM")
countToFPKM = importr("countToFPKM")

# read and translate R functions
f = open("/home/jupyteruser/work/py/rscripts/rnaseq.R", "r")
string = f.read()
f.close()
rnaseq_io = SignatureTranslatedAnonymousPackage(string, "rnaseq_io")


def load_rnaseq_tests(filename, model_name):
    """
    Load rnaseq results returning a dictionary of test (context, tissue, cell, etc ) names and rnaseq expression data
    """

    if not filename or filename=="None":
        tests = ["dummy"]
        fullsavepath = os.path.join(configs.rootdir, 'data', 'config_sheets', 'dummy_rnaseq_data.csv')
        data = pd.read_csv(fullsavepath, index_col='ENTREZ_GENE_ID')
        datas = [data]
        rnaseq_dict = dict(zip(tests, datas))

        return rnaseq_dict

    inqueryFullPath = os.path.join(configs.rootdir, 'data', 'config_sheets', filename)
    if not os.path.isfile(inqueryFullPath):
        print('Error: file not found {}'.format(inqueryFullPath))

        return None

    tests = []
    datas = []

    filename = 'rnaseq_{}.csv'.format(model_name)
    fullsavepath = os.path.join(configs.rootdir, 'data', 'results', model_name, filename)
    data = pd.read_csv(fullsavepath, index_col='ENTREZ_GENE_ID')
    print('Read from {}'.format(fullsavepath))
    datas.append(data)
    tests.append(model_name)

    rnaseq_dict = dict(zip(tests, datas))

    return rnaseq_dict


def handle_tissue_batch(config_filename, replicate_ratio, sample_ratio, replicate_ratio_high,
                        sample_ratio_high, technique, quantile, min_count):
    """
    Handle iteration through each tissue type and create rnaseq expression file by calling rnaseq.R
    """

    rnaseq_config_filepath = os.path.join(configs.rootdir, "data", "config_sheets", config_filename)
    xl = pd.ExcelFile(rnaseq_config_filepath)
    sheet_names = xl.sheet_names

    for model_name in sheet_names:
        print("model: ", model_name)
        rnaseq_output_file = "".join(["rnaseq_", model_name, ".csv"])
        rnaseq_output_filepath = os.path.join(configs.rootdir, "data", "results",
                                              model_name, rnaseq_output_file)
        rnaseq_input_file = "".join(["gene_counts_matrix_", model_name, ".csv"])
        rnaseq_input_filepath = os.path.join(configs.rootdir, "data", "data_matrices",
                                             model_name, rnaseq_input_file)
        gene_info_file = "".join(["gene_info_", model_name, ".csv"])
        gene_info_filepath = os.path.join(configs.rootdir, "data", "results",
                                          model_name, gene_info_file)

        os.makedirs(os.path.dirname(rnaseq_output_filepath), exist_ok=True)
        print('Input count matrix is at "{}"'.format(rnaseq_input_filepath))
        print('Gene info file is at "{}"'.format(gene_info_filepath))

        rnaseq_io.save_rnaseq_tests(rnaseq_input_filepath, rnaseq_config_filepath,
                                    rnaseq_output_filepath, gene_info_filepath,
                                    replicate_ratio=replicate_ratio, sample_ratio=sample_ratio,
                                    replicate_ratio_high=replicate_ratio_high, sample_ratio_high=sample_ratio_high,
                                    technique=technique, quantile=quantile,
                                    min_count=min_count, model_name=model_name)

        print("Test data saved to " + rnaseq_output_filepath)

    return


def main(argv):
    """
    Generate a list of active and high-confidence genes from a counts matrix using a user defined
    at normalization-technique at /work/data/results/<tissue name>/rnaseq_<tissue_name>.csv

    Currently, can filter raw RNA-seq counts using three normalization techniques. Which are defined in rnaseq.R

    TPM Quantile, where each replicate is normalized with Transcripts-per-million and an upper quantile is taken to
    create a boolean list of active genes for the replicate. Replicates are compared for consensus within the
    study/batch number according to user-defined ratios and then study/batch numbers are checked for consensus
    according to different user defined ratios.   **CITATION NEEDED** **Recomended if user wants more control over the
    size of the model, like a smaller model that allows for only the most expressed reactions, or a larger more
    encompassing one that contains less essential reactions.

    zFPKM method outlined in: https://pubmed.ncbi.nlm.nih.gov/24215113/ can be used. Counts will be normalized using
    zFPKM and genes > -3 will be considered expressed per thier recommendation. Expressed genes will be checked for
    consensus at the replicate and study/batch levels the same as TPM Quantile. **Recommended if user wants to give
    least input over gene essentially determination and use the most standardized method of active gene determination.

    flat cutoff of CPM (counts per million) normalized values, check for consensus the same as other methods.
    """

    parser = argparse.ArgumentParser(
        prog="rnaseq_gen.py",
        description="Generate a list of active and high-confidence genes from a counts matrix using a user defined "
                    "at normalization-technique at /work/data/results/<tissue name>/rnaseq_<tissue_name>.csv: "
                    "https://github.com/HelikarLab/FastqToGeneCounts",
        epilog="For additional help, please post questions/issues in the MADRID GitHub repo at "
               "https://github.com/HelikarLab/MADRID or email babessell@gmail.com",
        usage="python3 $(prog)s [options]"
    )

    parser.add_argument("-c", "--config-file",
                        type=str,
                        required=True,
                        dest="config_filename",
                        help="Name of config .xlsx file in the /work/data/config_files/. Can be generated using "
                             "rnaseq_preprocess.py or manually created and imported into the Juypterlab"
                        )
    
    parser.add_argument("-r", "--replicate-ratio",
                        type=float,
                        required=False,
                        default=0.5,
                        dest="replicate_ratio",
                        help="Ratio of replicates required for a gene to be active within that study/batch group "
                             "Example: 0.7 means that for a gene to be active, at least 70% of replicates in a group "
                             "must pass the cutoff after normalization"
                        )
    
    parser.add_argument("-g", "--group-ratio",
                        type=float,
                        required=False,
                        default=0.5,
                        dest="sample_ratio",
                        help="Ratio of groups (studies or batches) required for a gene to be active "
                             "Example: 0.7 means that for a gene to be active, at least 70% of groups in a study must  "
                             "have passed the replicate ratio test"
                        )

    parser.add_argument("-rh", "--high-replicate-ratio",
                        type=float,
                        required=False,
                        default=1.,
                        dest="replicate_ratio_high",
                        help="Ratio of replicates required for a gene to be considered high-confidence. "
                             "High-confidence genes ignore consensus with other data-sources like proteomics or "
                             "microarray. Example: 0.9 means that for a gene to be high-confidence, at least 90% of "
                             "replicates in a group must pass the cutoff after normalization"
                        )

    parser.add_argument("-gh", "--high-group-ratio",
                        type=float,
                        required=False,
                        default=1.,
                        dest="sample_ratio_high",
                        help="Ratio of groups (studies/batches) required for a gene to be considered high-confidence " 
                             "within that group. HHigh-confidence genes ignore consensus with other data-sources like "
                             "proteomics or microarray. Example: 0.9 means that for a gene to be high-confidence, at "
                             "least 90% of groups in a study must have passed the replicate ratio test"
                        )

    parser.add_argument("-t", "--filt-technique",
                        type=str,
                        required=False,
                        default="quantile-tpm",
                        dest="technique",
                        help="Technique to normalize and filter counts with. Either 'zfpkm', 'quantile-tpm' or "  
                             "'flat-cpm'. More info about each method is discussed in pipeline.ipynb."
                        )

    parser.add_argument("-q", "--quantile",
                        type=int,
                        required=False,
                        default=25,
                        dest="quantile",
                        help="Cutoff used for quantile-tpm normalization and filtration technique. Example: 25 means "
                             "that genes with TPM > 75% percentile wik=ll be considered active for that replicate."
                        )

    parser.add_argument("-m", "--min-count",
                        required=False,
                        default="default",
                        dest="min_count",
                        help="Cutoff used for cpm. Minimum number of counts to be considered expressed, alternatively "
                             "use 'default' to use method outlined in CITATION NEEDED"
                        )

    args = parser.parse_args(argv)

    config_filename = args.config_filename
    replicate_ratio = args.replicate_ratio
    sample_ratio = args.sample_ratio
    replicate_ratio_high = args.replicate_ratio_high
    sample_ratio_high = args.sample_ratio_high
    technique = args.technique
    quantile = args.quantile
    min_count = args.min_count

    if re.search("tpm", technique.lower()) or re.search("quantile", technique.lower()):
        technique = "quantile"
    elif re.search("cpm", technique.lower()):
        technique = "cpm"
    elif re.search("zfpkm", technique.lower()):
        technique = "zfpkm"
    else:
        print("Normalization-filtration technique not recognized. Must be 'tpm-quantile', 'cpm', or 'zfpkm'.")
        sys.exit()

    if int(quantile) > 100 or int(quantile) < 1:
        print("Quantile must be between 1 - 100")

    print('Config file is "{}"'.format(config_filename))

    handle_tissue_batch(config_filename, replicate_ratio, sample_ratio, replicate_ratio_high,
                        sample_ratio_high, technique, quantile, min_count)

    return


if __name__ == "__main__":
   main(sys.argv[1:])
