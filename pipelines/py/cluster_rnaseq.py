# !/usr/bin/python3

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
ggplot2 = importr("ggplot2")
ggrepel = importr("ggrepel")
tidyverse = importr("tidyverse")
#factoMineR = importr("FactoMineR")
uwot = importr("uwot")

# read and translate R functions
f = open("/home/jupyteruser/work/py/rscripts/cluster_samples.R", "r")
string = f.read()
f.close()
cluster_io = SignatureTranslatedAnonymousPackage(string, "cluster_io")


def main(argv):
    """
    Cluster RNA-seq Data
    """

    parser = argparse.ArgumentParser(
        prog="cluster_rnaseq.py",
        description="Cluster RNA-seq Data using Multiple Correspondence Analysis or UMAP. Clusters at the replicate, "
                    "batch/study, and context levels.",
        epilog="For additional help, please post questions/issues in the MADRID GitHub repo at "
               "https://github.com/HelikarLab/MADRID or email babessell@gmail.com",
    )

    parser.add_argument("-c", "--context-names",
                        type=str,
                        required=True,
                        dest="context_names",
                        help="""Tissue/cell name of models to generate. If making multiple models in a batch, then
                             use the format: \"['context1', 'context2', ... etc]\". Note the outer double-quotes and the 
                             inner single-quotes are required to be interpreted. This a string, not a python list"""
                        )

    parser.add_argument("-t", "--filt-technique",
                        type=str,
                        required=True,
                        dest="technique",
                        help="'zfpkm', 'quantile', or 'cpm'"
                        )

    parser.add_argument("-a", "--cluster-algorithm",
                        type=str,
                        required=False,
                        default="umap",
                        dest="clust_algo",
                        help="""Clustering algorithm to use. 'mca' or 'umap'."""
                        )

    parser.add_argument("-r", "--replicate-ratio",
                        type=str,
                        required=False,
                        default=0.9,
                        dest="rep_ratio",
                        help="""Ratio of genes active in replicates for a batch/study to be active"""
                        )

    parser.add_argument("-b", "--batch-ratio",
                        type=str or float,
                        required=False,
                        default=0.9,
                        dest="batch_ratio",
                        help="""Ratio of genes active in a batch/study to be active in the context"""
                        )

    parser.add_argument("-m", "--min-count",
                        type=str or int,
                        required=False,
                        default="default",
                        dest="min_count",
                        help="""Ratio of genes active in a batch/study to be active in the context"""
                        )

    parser.add_argument("-q", "--quantile",
                        type=str or int,
                        required=False,
                        default=0.5,
                        dest="quantile",
                        help="""Ratio of genes active in a batch/study to be active in the context"""
                        )

    args = parser.parse_args(argv)

    wd = os.path.join(configs.datadir, "results")
    context_names = args.context_names.strip("[").strip("]").replace("'", "").replace(" ", "").split(",")
    technique = args.technique.lower()
    clust_algo = args.clust_algo.lower()
    rep_ratio = args.rep_ratio
    batch_ratio = args.batch_ratio
    min_count = args.min_count
    quantile = args.quantile


    if type(min_count)==str and not min_count.lower() == "default":
        try:
            min_count = int(min_count)
        except ValueError:
            print(f"--min-count must be either 'default' or an integer > 0")
            sys.exit()

    if (type(quantile)==str and not quantile.lower() == "default") or 0 < quantile > 100:
        try:
            quantile = int(quantile)
        except ValueError:

            print(f"--quantile must be either 'default' or an integer between 0 and 100")
            sys.exit()

    if (type(rep_ratio)==str and not rep_ratio.lower() == "default") or 0 < rep_ratio > 1.0:
        try:
            rep_ratio = float(rep_ratio)
        except ValueError:
            print("--rep-ratio must be 'default' or a float between 0 and 1")
            sys.exit()

    if (type(batch_ratio) == str and not batch_ratio.lower() == "default") or 0 < rep_ratio > 1.0:
        try:
            batch_ratio = float(batch_ratio)
        except ValueError:
            print("--batch-ratio must be 'default' or a float between 0 and 1")
            sys.exit()

    if technique.lower() not in ["quantile", "tpm", "cpm", "zfpkm"]:
        print("--technique must be either 'quantile', 'cpm', 'zfpkm'")
        sys.exit()

    if technique.lower() == "tpm":
        technique = "quantile"

    if clust_algo.lower() not in ["mca", "umap"]:
        print("--technique must be either 'umap', 'mca'")
        sys.exit()

    cluster_io.cluster_samples_main(wd, context_names, technique, clust_algo,  n_neigh_rep="default",
                                    n_neigh_batch="default", n_neigh_context="default", rep_ratio=rep_ratio,
                                    batch_ratio=batch_ratio, quantile=quantile, min_count=min_count)
    return


if __name__ == "__main__":
    main(sys.argv[1:])
