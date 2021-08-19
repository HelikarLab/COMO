#!/usr/bin/python3

import os, time, sys
import pandas as pd
import getopt
from rpy2.robjects.packages import importr
from rpy2.robjects import r, pandas2ri
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage
from project import configs
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

# read config file
def load_bulk_supplementary_data(suppfilename):
    if not suppfilename or suppfilename=="None":
        return "None"
    suppFullPath = os.path.join(configs.rootdir, 'data', suppfilename)
    supplements = pd.read_csv(suppFullPath, header=0)
    supplements = supplements[supplements['SampleName'].str.match("FILENAME")]
    Bulk = supplements['InsertSize']
    return Bulk

# read data from csv files
def load_bulk_tests(Bulk):
    try:
        if not Bulk or Bulk=="None":
            tests = ["dummy"]
            fullsavepath = os.path.join(configs.rootdir, 'data', "dummy_data.csv")
            datas = ["dummy_data"]
            bulk_dict = dict(zip(tests, datas))
            return bulk_dict
    except:
        print("bulk exists")

    tests = []
    datas = []
    for test in list(Bulk):
        if test == 'Statistics':
            break
        test = unidecode.unidecode(test)
        output_path = os.path.join(configs.rootdir, 'data', 'Bulk_{}.csv'.format(test.strip()))
        testdata = pd.read_csv(output_path, index_col=False)
        #testdat = testdata.applymap(str)
        testdata.drop_duplicates(inplace=True)
        testdata['ENTREZ_GENE_ID'] = testdata['ENTREZ_GENE_ID'].astype(object)
        testdata.set_index('ENTREZ_GENE_ID', inplace=True, drop=True)
        print('Test Data Load From {}'.format(output_path))
        tests.append(test.strip().replace('ï', 'i'))
        datas.append(testdata)

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
        opts, args = getopt.getopt(argv, "hi:c:g:r:s:x:y:t:q:m:",
                ["datafile=", "suppfile=", "gene_format=", "species_dataset=",
                 "expr_prop_rep=", "expr_prop_samp=", "top_percentile=",
                 "technique=", "quantile=", "min_count="])

    except getopt.GetoptError:
        err_str = """
python3 proteomics_gen.py -f <data file> -c <config file>\n
-g <gene format> -d <species_dataset> -r <replicate proportion>\n
-s <sample proportion> -p <top percentile> -t <filtration technique>\n
-q <cutoff quantile (for quantile technique)> -m <min count (for cpm technique)\n
              """
        print(err_str)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            help_str = """
python3 proteomics_gen.py -f <data file> -c <config file>\n
-g <gene format> -d <species_dataset> -r <replicate proportion>\n
-s <sample proportion> -p <top percentile> -t <filtration technique>\n
-q <cutoff quantile (for quantile technique)> -m <min count (for cpm technique)\n
              """
            print(help_str)
            sys.exit()
        elif opt in ("-i", "--datafile"):
            datafile = arg
        elif opt in ("-c", "--suppfile"):
            suppfile = arg
        elif opt in ('-g', '--gene_info_file'):
            gene_info_file = arg
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
    
    print('Data file is "{}"'.format(datafile))
    print('Supplementary Data file is "{}"'.format(suppfile))
    print('Gene info file is "{}"'.format(gene_info_file))

    bulk_config_filepath = os.path.join(configs.rootdir, "data", suppfile)
    bulk_input_filepath = os.path.join(configs.rootdir, "data", datafile)
    gene_info_filepath = os.path.join(configs.rootdir, "data", gene_info_file)
    #if not os.path.isfile(prote_data_filepath):
    #    proteomics_data.to_csv(prote_data_filepath, index_label='ENTREZ_GENE_ID')

    #config_df = pd.read_csv(bulk_config_filepath)
    xl = pd.ExcelFile(bulk_config_filepath)
    sheet_names = xl.sheet_names
    #model_name = "".join(["Bulk_",config_df["InsertSize"][0] ,".csv"])
    for model_name in sheet_names:
        bulk_filename = "".join(["Bulk_", model_name, ".csv"])
        bulk_output_filepath = os.path.join(configs.rootdir, "data", bulk_filename)
        print('Output File is "{}"'.format(bulk_output_filepath))

        bulkio.save_bulk_tests(bulk_input_filepath, bulk_config_filepath,
                               bulk_output_filepath, gene_info_filepath, 
                               pos_rep=pos_rep, pos_samp=pos_samp,
                               top_rep=top_rep, top_samp=top_samp,
                               technique=technique, quantile=quantile,
                               min_count=min_count)
    
        print("Test data saved to " + bulk_output_filepath)
    # save proteomics data by test
    #proteomics_dict, testdata_dict = save_proteomics_tests(Proteomics, proteomics_data, expr_prop, percentile)
    return True


if __name__ == "__main__":
   main(sys.argv[1:])