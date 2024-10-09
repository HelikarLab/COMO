import argparse
from pathlib import Path

import numpy as np
import rpy2_api
from como_utilities import stringlist_to_list
from project import Config

# read and translate R functions
configs = Config()
r_file_path = Path(configs.code_dir, "rscripts", "cluster_samples.R")


def main() -> None:
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
        help="""Ratio of active genes in a batch/study to be active in the context""",
    )
    parser.add_argument(
        "-q",
        "--quantile",
        type=str or int,
        required=False,
        default=0.5,
        dest="quantile",
        help="""Ratio of active genes in a batch/study to be active in the context""",
    )
    parser.add_argument(
        "-s",
        "--seed",
        type=int,
        required=False,
        default=-1,
        dest="seed",
        help="""Random seed for clustering algorithm initialization""",
    )
    args = parser.parse_args()

    context_names = stringlist_to_list(args.context_names)
    technique = args.technique.lower()
    clust_algo = args.clust_algo.lower()
    label = args.label
    rep_ratio = args.rep_ratio
    batch_ratio = args.batch_ratio
    min_count = args.min_count
    quantile = args.quantile
    min_dist = args.min_dist
    n_neigh_rep = args.n_neigh_rep
    n_neigh_batch = args.n_neigh_batch
    n_neigh_cont = args.n_neigh_cont

    # Set a random seed if none provided
    if int(args.seed) == -1:
        seed = np.random.randint(0, 100000)
    else:
        seed = args.seed

    if isinstance(min_count, str) and min_count.lower() == "default":
        try:
            min_count = int(min_count)
        except ValueError:
            raise ValueError("--min-count must be either 'default' or an integer > 0")
    if not isinstance(min_count, str) and min_count < 0:
        raise ValueError("--min-count must be either 'default' or an integer > 0")

    if isinstance(quantile, str) and not quantile.lower() == "default":
        try:
            quantile = int(quantile)
        except ValueError:
            raise ValueError("--quantile must be either 'default' or an integer between 0 and 100")
    if not isinstance(quantile, str) and 0 > quantile > 100:
        raise ValueError("--quantile must be either 'default' or an integer between 0 and 100")

    if isinstance(rep_ratio, str) and not rep_ratio.lower() == "default":
        try:
            rep_ratio = float(rep_ratio)
        except ValueError:
            raise ValueError("--rep-ratio must be 'default' or a float between 0 and 1")
    if not isinstance(rep_ratio, str) and 0 > rep_ratio > 1.0:
        raise ValueError("--rep-ratio must be 'default' or a float between 0 and 1")

    if isinstance(batch_ratio, str) and not batch_ratio.lower() == "default":
        try:
            batch_ratio = float(batch_ratio)
        except ValueError:
            raise ValueError("--batch-ratio must be 'default' or a float between 0 and 1")
    if not isinstance(batch_ratio, str) and 0 > batch_ratio > 1.0:
        raise ValueError("--batch-ratio must be 'default' or a float between 0 and 1")

    if technique.lower() not in ["quantile", "tpm", "cpm", "zfpkm"]:
        raise ValueError("--technique must be either 'quantile', 'tpm', 'cpm', 'zfpkm'")

    if technique.lower() == "tpm":
        technique = "quantile"

    if clust_algo.lower() not in ["mca", "umap"]:
        raise ValueError("--clust_algo must be either 'mca', 'umap'")

    if not isinstance(min_dist, str) and 0 > min_dist > 1.0:
        raise ValueError("--min_dist must be a float between 0 and 1")

    if isinstance(n_neigh_rep, str) and not n_neigh_rep.lower() == "default":
        try:
            n_neigh_rep = int(n_neigh_rep)
        except ValueError:
            raise ValueError(
                "--n_neigh_rep must be either 'default' or an integer greater than 1 and less than or equal to "
                "the total number of replicates being clustered across all contexts."
            )
    if not isinstance(n_neigh_rep, str) and n_neigh_rep < 2:
        raise ValueError(
            "--n_neigh_rep must be either 'default' or an integer greater than 1 and less than or equal to "
            "the total number of replicates being clustered across all contexts."
        )

    if isinstance(n_neigh_batch, str) and not n_neigh_batch.lower() == "default":
        try:
            n_neigh_batch = int(n_neigh_batch)
        except ValueError:
            raise ValueError(
                "--n_neigh_batch must be either 'default' or an integer greater than 1 and less than or equal to "
                "the total number of batches being clustered across all contexts."
            )
    if not isinstance(n_neigh_batch, str) and n_neigh_batch < 2:
        raise ValueError(
            "--n_neigh_batch must be either 'default' or an integer greater than 1 and less than or equal to "
            "the total number of batches being clustered across all contexts."
        )

    if isinstance(n_neigh_cont, str) and not n_neigh_cont.lower() == "default":
        try:
            n_neigh_cont = int(n_neigh_cont)
        except ValueError:
            raise ValueError(
                "--n_neigh_batch must be either 'default' or an integer greater than 1 and less than or equal to "
                "the total number of batches being clustered across all contexts."
            )
    if not isinstance(n_neigh_cont, str) and n_neigh_cont < 2:
        raise ValueError(
            "--n_neigh_context must be either 'default' or an integer greater than 1 and less than or equal to "
            "the total number of contexts being clustered."
        )

    cluster_samples = rpy2_api.Rpy2(
        r_file_path=r_file_path,
        wd=configs.result_dir,
        context_names=context_names,
        technique=technique,
        clust_algo=clust_algo,
        label=label,
        min_dist=min_dist,
        n_neigh_rep=n_neigh_rep,
        n_neigh_batch=n_neigh_batch,
        n_neigh_cont=n_neigh_cont,
        rep_ratio=rep_ratio,
        batch_ratio=batch_ratio,
        quantile=quantile,
        min_count=min_count,
        seed=seed,
    )
    cluster_samples.call_function("cluster_samples_main")


if __name__ == "__main__":
    main()
