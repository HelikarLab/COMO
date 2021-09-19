#!/usr/bin/python3
import os, time, sys
#sys.stdout = open("/home/jupyteruser/work/py/rlogs/sys.out", 'w')
import pandas as pd
import getopt
from rpy2.robjects.packages import importr
from rpy2.robjects import r, pandas2ri
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage
from project import configs
from contextlib import redirect_stderr
import unidecode

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

# read and translate R functions
f = open("/home/jupyteruser/work/py/rscripts/bulk.R", "r")
string = f.read()
f.close()
bulkio = SignatureTranslatedAnonymousPackage(string, "bulkio")

# read data from csv files
def load_bulk_tests(filename):
    if not filename or filename=="None":
        tests = ["dummy"]
        fullsavepath = os.path.join(configs.rootdir, 'data', 'config_sheets', 'dummy_bulk_data.csv')
        data = pd.read_csv(fullsavepath, index_col='ENTREZ_GENE_ID')
        datas = [data]
        bulk_dict = dict(zip(tests, datas))
        return bulk_dict

    inqueryFullPath = os.path.join(configs.rootdir, 'data', 'config_sheets', filename)
    if not os.path.isfile(inqueryFullPath):
        print('Error: file not found {}'.format(inqueryFullPath))
        return None
    xl = pd.ExcelFile(inqueryFullPath)
    sheet_names = xl.sheet_names

    tests = []
    datas = []
    for model_name in sheet_names:
        # print(list(inqueries[i]))
        filename = 'Bulk_{}.csv'.format(model_name)
        fullsavepath = os.path.join(configs.rootdir, 'data', 'results', model_name, filename)
        data = pd.read_csv(fullsavepath, index_col='ENTREZ_GENE_ID')
        print('Read from {}'.format(fullsavepath))
        datas.append(data)
        tests.append(model_name)
    bulk_dict = dict(zip(tests, datas))
    return bulk_dict

def main(argv):
    #default vals
    pos_rep = 0.5
    pos_samp = 0.5
    top_rep = 0.5
    top_samp = 0.5
    quantile = 0.9
    min_count = 10

    try:
        opts, args = getopt.getopt(argv, "hc:r:s:x:y:t:q:m:",
                ["suppfile=", "expr_prop_rep=", "expr_prop_samp=", 
                 "top_prop_rep", "top_prop_samp", "technique=", 
                 "quantile=", "min_count="])

    except getopt.GetoptError:
        err_str = """
        python3 bulk_gen.py -c <config file> -r <replicate expression ratio> \n
        -s <sample expression ratio> -x <replicate high-confidence ratio> \n
        -y < sample high-confidence ratio> -t <filter technique> \n
        -q <cutoff quantile> -m <cutoff count>
              """
        print(err_str)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            help_str = """
            python3 bulk_gen.py -c <config file> -r <replicate expression ratio> \n
            -s <sample expression ratio> -x <replicate high-confidence ratio> \n
            -y < sample high-confidence ratio> -t <filter technique> \n
            -q <cutoff quantile> -m <cutoff count>
              """
            print(help_str)
            sys.exit()
        elif opt in ("-c", "--suppfile"):
            suppfile = arg
        elif opt in ("-r", "--expr_prop_rep"):
            pos_rep = float(arg)
        elif opt in ("-s", "--expr_prop_samp"):
            pos_samp = float(arg)
        elif opt in ("-x", "--top_prop_rep"):
            top_rep = float(arg)
        elif opt in ("-y", "--top_prop_samp"):
            top_samp = float(arg)
        elif opt in ("-t", "--technique"):
            technique = arg
        elif opt in ("-q", "--quantile"):
            quantile = int(arg)
        elif opt in ("-m", "--min_count"):
            min_count = int(arg)
    
    print('Config file is "{}"'.format(suppfile))

    bulk_config_filepath = os.path.join(configs.rootdir, "data", "config_sheets", suppfile)
    xl = pd.ExcelFile(bulk_config_filepath)
    sheet_names = xl.sheet_names

    for model_name in sheet_names:
        bulk_output_file = "".join(["Bulk_", model_name, ".csv"])
        bulk_output_filepath = os.path.join(configs.rootdir, "data", "results",
                                            model_name, bulk_output_file)
        bulk_input_file = "".join(["BulkRNAseqDataMatrix_", model_name, ".csv"])
        bulk_input_filepath = os.path.join(configs.rootdir, "data", "data_matrices", 
                                           model_name, bulk_input_file)
        gene_info_file = "".join(["GeneInfo_", model_name, ".csv"])
        gene_info_filepath = os.path.join(configs.rootdir, "data", "results",
                                          model_name, gene_info_file)
        
        os.makedirs(os.path.dirname(bulk_output_filepath), exist_ok=True)
        print('Input count matrix is at "{}"'.format(bulk_input_filepath))
        print('Gene info file is at "{}"'.format(gene_info_filepath))

        bulkio.save_bulk_tests(bulk_input_filepath, bulk_config_filepath,
                               bulk_output_filepath, gene_info_filepath, 
                               pos_rep=pos_rep, pos_samp=pos_samp,
                               top_rep=top_rep, top_samp=top_samp,
                               technique=technique, quantile=quantile,
                               min_count=min_count, model_name=model_name)
    
        print("Test data saved to " + bulk_output_filepath)

    return True


if __name__ == "__main__":
   main(sys.argv[1:])
