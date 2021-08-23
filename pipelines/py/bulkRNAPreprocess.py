#!/usr/bin/python3
from bioservices import BioDBNet
import pandas as pd
from project import configs
import re
import os, time, sys
import getopt
from rpy2.robjects.packages import importr, SignatureTranslatedAnonymousPackage
#from rpy2.robjects import r, pandas2ri
#import rpy2.robjects as ro
#from rpy2.robjects.conversion import localconverter
#from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage

#pandas2ri.activate()

#limma = importr("limma")
tidyverse = importr("tidyverse")
#edgeR = importr("edgeR")
#genefilter = importr("genefilter")
#biomaRt = importr("biomaRt")
#sjmisc = importr("sjmisc")

# automatically convert ryp2 dataframe to Pandas dataframe
f = open("/home/jupyteruser/work/py/rscripts/genCountMatrix.R", "r")
string = f.read()
f.close()

genCountMatrixio = SignatureTranslatedAnonymousPackage(string, "genCountMatrixio")

def fetch_gene_info(input_values, input_db='Ensembl Gene ID',
                    output_db=['Gene Symbol','Gene ID','Chromosomal Location'],
                    delay=15):
 
    s = BioDBNet()   
    # input_db = 'Agilent ID'
    # input_values = df_results.ProbeName.tolist()

    df_maps = pd.DataFrame([],columns=output_db)
    df_maps.index.name=input_db
    i = 0
    # for i in range(0,len(input_values),500):
    while i < len(input_values):
        print('retrieve {}:{}'.format(i,min(i+500,len(input_values))))
        df_test = s.db2db(input_db, output_db, input_values[i:min(i+500,len(input_values))], 9606)
        if isinstance(df_test, pd.DataFrame):
            df_maps = pd.concat([df_maps, df_test], sort=False)
        elif df_test == '414':
            print("bioDBnet busy, try again in {} seconds".format(delay))
            time.sleep(delay)
            continue
        i += 500
    return df_maps

def main(argv):
    try:
        opts, args = getopt.getopt(argv, "hn:c:f:t:",
                                   ["tissue_name=",
                                    "create_counts_matrix=",
                                    "gene_format=",
                                    "technique="])
    except getopt.GetoptError:
        errstr = """
        python3 bulkRNAPreprocess.py -n <tissue_names> -c <create_counts_matrix> \n
        -f <gene_format> -t <technique>
        """
        print(errstr)
        sys.exit(2)
    for opt, arg in opts:
        if opt == '-h':
            helpstr = """
            python3 bulkRNAPreprocess.py -n <tissue_names> -c <create_counts_matrix> \n
            -f <gene_format> -t <technique>
            """
            print(helpstr)
            sys.exit()
        elif opt in ("-n", "--tissue_names"):
            tissue_names = arg
        elif opt in ("-c", "--create_counts_matrix"):
            make_matrix = arg
        elif opt in ("-f", "--gene_format"):
            gene_format = arg  
        elif opt in ("-t", "--technique"):
            technique = arg
            
    tissue_names = tissue_names.strip("[").strip("]").replace("'", "").split(",")
    #tissue_names = tissue_names.split(",")
    for tissue_name in tissue_names:
        print(tissue_name)
        input_dir = os.path.join(configs.rootdir, 'data', 'STAR_output', tissue_name)
        gene_output_dir = os.path.join(configs.rootdir, 'data', 'results', tissue_name)
        matrix_output_dir = os.path.join(configs.rootdir, 'data', 'data_matrices', tissue_name)
        os.makedirs(gene_output_dir, exist_ok=True)
        os.makedirs(matrix_output_dir, exist_ok=True) 
        print('Input directory is "{}"'.format(input_dir))
        print('Gene info output directory is "{}"'.format(gene_output_dir))
        print('Active gene determination technique is "{}"'.format(technique))
        
        if make_matrix:
            matrix_output_dir = os.path.join(configs.rootdir, 'data', 'data_matrices', tissue_name)  
            print("Creating Counts Matrix")
            genCountMatrixio.genCountMatrix_main(input_dir, matrix_output_dir, technique)
            
        geneCountFile = os.path.join(matrix_output_dir, ("BulkRNAseqDataMatrix_"+tissue_name+".csv"))
        print('Fetching gene info using genes in "{}"'.format(geneCountFile))
        genes = pd.read_csv(geneCountFile)['genes'].to_list()
        output_db=['Ensembl Gene ID', 'Gene Symbol', 'Gene ID', 'Chromosomal Location']               
        if gene_format.upper()=="ENSEMBL":
             form = "Ensembl Gene ID"                     
        elif gene_format.upper()=="ENTREZ":
             form = "Gene ID"
        elif gene_format.upper()=="SYMBOL":
              form = "Gene Symbol"                       
        output_db.remove(form)  
        gene_info = fetch_gene_info(genes, input_db=form, output_db=output_db)
        gene_info['start_position'] = gene_info['Chromosomal Location'].str.extract("chr_start: (\d+)")
        gene_info['end_position'] = gene_info['Chromosomal Location'].str.extract("chr_end: (\d+)")
        gene_info.index.rename("ensembl_gene_id", inplace=True)
        gene_info.rename(columns={"Gene Symbol": "hgnc_symbol", "Gene ID": "entrezgene_id"}, inplace=True)
        gene_info.drop(['Chromosomal Location'], axis=1, inplace=True)
        gene_info_file = os.path.join(gene_output_dir, ("GeneInfo_"+tissue_name+".csv"))
        gene_info.to_csv(gene_info_file)
        print('Gene Info file written at "{}"'.format(gene_info_file))
        
    return True


if __name__ == "__main__":
    print(sys.argv)
    main(sys.argv[1:])
