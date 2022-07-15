# !/usr/bin/python3

import os
import sys
from rpy2.robjects.packages import importr
from rpy2.robjects import r
from rpy2.robjects import pandas2ri
from rpy2.robjects.packages import SignatureTranslatedAnonymousPackage
from project import configs
import argparse

# enable r to py conversion
pandas2ri.activate()

# import R libraries
ggplot2 = importr("ggplot2")
ggrepel = importr("ggrepel")
tidyverse = importr("tidyverse")
# factoMineR = importr("FactoMineR")
uwot = importr("uwot")

# read and translate R functions
f = open(os.path.join(configs.rootdir, "py", "rscripts", "cluster_samples.R"), "r")
string = f.read()
f.close()
cluster_io = SignatureTranslatedAnonymousPackage(string, "cluster_io")


def parse_args(argv) -> argparse.Namespace:
    """
    This function is responsible for parsing arguments as they are received from the command line
    
    :param argv: The list of arguments directly from the command line
    :return: The parsed arguments in the form of an argparse.Namespace object
    """
    parser = argparse.ArgumentParser(
        prog="cluster_rnaseq.py",
        description="Cluster RNA-seq Data using Multiple Correspondence Analysis or UMAP. Clusters at the replicate, "
        "batch/study, and context levels.",
        epilog="For additional help, please post questions/issues in the MADRID GitHub repo at "
        "https://github.com/HelikarLab/MADRID or email babessell@gmail.com",
    )
    parser.add_argument(
        "-n",
        "--context-names",
        type=str,
        required=True,
        dest="context_names",
        help="""Tissue/cell name of models to generate. If making multiple models in a batch, then
                             use the format: \"['context1', 'context2', ... etc]\". Note the outer double-quotes and the
                             inner single-quotes are required to be interpreted. This a string, not a python list""",
    )
    parser.add_argument(
        "-t",
        "--filt-technique",
        type=str,
        required=True,
        dest="technique",
        help="'zfpkm', 'quantile', or 'cpm'",
    )
    parser.add_argument(
        "-a",
        "--cluster-algorithm",
        type=str,
        required=False,
        default="umap",
        dest="clust_algo",
        help="""Clustering algorithm to use. 'mca' or 'umap'.""",
    )
    parser.add_argument(
        "-l",
        "--label",
        type=str,
        required=False,
        default=True,
        dest="label",
        help="""True to label replicate/batch/context names on the plots. May be ugly for large sets""",
    )
    parser.add_argument(
        "-d",
        "--min-dist",
        type=float,
        required=False,
        default=0.01,
        dest="min_dist",
        help="""Minimum distance for UMAP clustering. Must be between 0 and 1""",
    )
    parser.add_argument(
        "-r",
        "--replicate-ratio",
        type=str,
        required=False,
        default=0.9,
        dest="rep_ratio",
        help="""Ratio of genes active in replicates for a batch/study to be active""",
    )
    parser.add_argument(
        "-b",
        "--batch-ratio",
        type=str or float,
        required=False,
        default=0.9,
        dest="batch_ratio",
        help="""Ratio of genes active in a batch/study to be active in the context""",
    )
    parser.add_argument(
        "-nr",
        "--n-neighbors-rep",
        type=str or float,
        required=False,
        default="default",
        dest="n_neigh_rep",
        help="""N nearest neighbors for replicate clustering, 'default' is total number of replicates""",
    )
    parser.add_argument(
        "-nb",
        "--n-neighbors-batch",
        type=str or float,
        required=False,
        default="default",
        dest="n_neigh_batch",
        help="""N nearest neighbors for batch clustering, 'default' is total number of batches""",
    )
    parser.add_argument(
        "-nc",
        "--n-neighbors-context",
        type=str or float,
        required=False,
        default="default",
        dest="n_neigh_cont",
        help="""N nearest neighbors for context clustering, 'default' is total number of contexts""",
    )
    parser.add_argument(
        "-c",
        "--min-count",
        type=str or int,
        required=False,
        default="default",
        dest="min_count",
        help="""Ratio of genes active in a batch/study to be active in the context""",
    )
    parser.add_argument(
        "-q",
        "--quantile",
        type=str or int,
        required=False,
        default=0.5,
        dest="quantile",
        help="""Ratio of genes active in a batch/study to be active in the context""",
    )
    parser.add_argument(
        "-s",
        "--seed",
        type=int,
        required=False,
        default=12345,
        dest="seed",
        help="""Random seed for clustering algorithm initialization""",
    )
    args = parser.parse_args(argv)
    return args


def validate_args(args: argparse.Namespace) -> argparse.Namespace:
    """
    This function is responsible for ensuring the parsed arguments are valid.
    
    :param args: The arguments parsed from the parse_args() function
    :return: The validated arguments in the form of an argparse.Namespace object
    """
    
    args.valid_arguments = True
    
    if type(args.min_count) == str and str(args.min_count).lower() != "default":
        try:
            args.min_count = int(args.min_count)
        except ValueError:
            print("--min-count must be either 'default' or an integer > 0")
            args.valid_arguments = False
    elif type(args.min_count) != str and args.min_count < 0:
        print("--min-count must be either 'default' or an integer > 0")
        args.valid_arguments = False

    if type(args.quantile) == str and args.quantile.lower() != "default":
        try:
            args.quantile = int(args.quantile)
        except ValueError:
            print("--quantile must be either 'default' or an integer between 0 and 100")
            args.valid_arguments = False
    elif type(args.quantile) != str and 0 > args.quantile > 1:
        print("--quantile must be either 'default' or an integer between 0 and 100")
        args.valid_arguments = False

    if type(args.rep_ratio) == str and str(args.rep_ratio).lower() != "default":
        try:
            args.rep_ratio = float(args.rep_ratio)
        except ValueError:
            print("--rep-ratio must be 'default' or a float between 0 and 1")
            args.valid_arguments = False
    elif type(args.rep_ratio) != str and 0 > args.rep_ratio > 1.0:
        print("--rep-ratio must be 'default' or a float between 0 and 1")

    if type(args.batch_ratio) == str and str(args.batch_ratio).lower() != "default":
        try:
            args.batch_ratio = float(args.batch_ratio)
        except ValueError:
            print("--batch-ratio must be 'default' or a float between 0 and 1")
            args.valid_arguments = False
    elif type(args.batch_ratio) != str and 0 > args.batch_ratio > 1.0:
        print("--batch-ratio must be 'default' or a float between 0 and 1")

    if args.technique.lower() not in ["quantile", "tpm", "cpm", "zfpkm"]:
        print("--technique must be either 'quantile', 'cpm', 'zfpkm'")
        args.valid_arguments = False
    elif args.technique.lower() == "tpm":
        args.technique = "quantile"

    if args.clust_algo.lower() not in ["mca", "umap"]:
        print("--technique must be either 'umap', 'mca'")
        args.valid_arguments = False

    if type(args.min_dist) != str and 0 > args.min_dist > 1.0:
        print("--min_dist must be a float between 0 and 1")
        args.valid_arguments = False

    if type(args.n_neigh_rep) == str and str(args.n_neigh_rep).lower() != "default":
        try:
            args.n_neigh_rep = int(args.n_neigh_rep)
        except ValueError:
            print(f"--n_neigh_rep must be either 'default' or an integer greater than 1 and less than or equal to "
                  f"the total number of replicates being clustered across all contexts.")
            args.valid_arguments = False
    elif type(args.n_neigh_rep) != str and args.n_neigh_rep < 2:
        print(f"--n_neigh_rep must be either 'default' or an integer greater than 1 and less than or equal to "
              f"the total number of replicates being clustered across all contexts.")
        args.valid_arguments = False

    if type(args.n_neigh_batch) == str and str(args.n_neigh_batch).lower() != "default":
        try:
            args.n_neigh_batch = int(args.n_neigh_batch)
        except ValueError:
            print(f"--n_neigh_batch must be either 'default' or an integer greater than 1 and less than or equal to "
                  f"the total number of batches being clustered across all contexts.")
            args.valid_arguments = False
    elif type(args.n_neigh_batch) != str and args.n_neigh_batch < 2:
        print(f"--n_neigh_batch must be either 'default' or an integer greater than 1 and less than or equal to "
              f"the total number of batches being clustered across all contexts.")
        args.valid_arguments = False

    if type(args.n_neigh_cont) == str and str(args.n_neigh_cont).lower() != "default":
        try:
            args.n_neigh_cont = int(args.n_neigh_cont)
        except ValueError:
            print(f"--n_neigh_batch must be either 'default' or an integer greater than 1 and less than or equal to "
                  f"the total number of batches being clustered across all contexts.")
            args.valid_arguments = False
            
    if type(args.n_neigh_batch) != str and args.n_neigh_cont < 2:
        print(f"--n_neigh_context must be either 'default' or an integer greater than 1 and less than or equal to "
              f"the total number of contexts being clustered.")
        args.valid_arguments = False
        
    return args


def main(argv):
    """
    Cluster RNA-seq Data
    """
    args = parse_args(argv)
    args = validate_args(args)
    
    if not args.valid_arguments:
        sys.exit(1)

    wd = os.path.join(configs.datadir, "results")
    context_names = (
        args.context_names.strip("[")
        .strip("]")
        .replace("'", "")
        .replace(" ", "")
        .split(",")
    )

    cluster_io.cluster_samples_main(
        wd,
        context_names,
        args.technique.lower(),
        args.clust_algo.lower(),
        args.label.upper(),
        min_dist=args.min_dist,
        n_neigh_rep=args.n_neigh_rep,
        n_neigh_batch=args.n_neigh_batch,
        n_neigh_cont=args.n_neigh_cont,
        rep_ratio=args.rep_ratio,
        batch_ratio=args.batch_ratio,
        quantile=args.quantile,
        min_count=args.min_count,
        seed=args.seed,
    )


if __name__ == "__main__":
    main(sys.argv[1:])
