# !/usr/bin/python3

# from rpy2.robjects.packages import importr
# from rpy2.robjects import pandas2ri
import argparse
import os
import sys
from pathlib import Path

import numpy as np
import rpy2_api
from arguments import (
    batch_ratio_arg,
    cluster_algorithm_arg,
    context_names_arg,
    filtering_technique_arg,
    label_arg,
    min_count_arg,
    min_dist_arg,
    num_neighbors_batch_arg,
    num_neighbors_context_arg,
    num_neighbors_replicate_arg,
    quantile_arg,
    random_seed_arg,
    replicate_ratio_arg,
)
from como_utilities import stringlist_to_list
from project import Configs

# enable r to py conversion
# pandas2ri.activate()

# import R libraries
# ggplot2 = importr("ggplot2")
# ggrepel = importr("ggrepel")
# tidyverse = importr("tidyverse")
# uwot = importr("uwot")

# read and translate R functions
configs = Configs()
r_file_path = Path(configs.root_dir, "src", "rscripts", "cluster_samples.R")


# f = open(os.path.join(configs.rootdir, "src", "rscripts", "cluster_samples.R"), "r")
# string = f.read()
# f.close()
# cluster_io = SignatureTranslatedAnonymousPackage(string, "cluster_io")


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

    parser.add_argument(**context_names_arg)
    parser.add_argument(**filtering_technique_arg)
    parser.add_argument(**cluster_algorithm_arg)
    parser.add_argument(**label_arg)
    parser.add_argument(**min_dist_arg)
    parser.add_argument(**replicate_ratio_arg)
    parser.add_argument(**batch_ratio_arg)
    parser.add_argument(**num_neighbors_replicate_arg)
    parser.add_argument(**num_neighbors_batch_arg)
    parser.add_argument(**num_neighbors_context_arg)
    parser.add_argument(**min_count_arg)
    parser.add_argument(**quantile_arg)
    parser.add_argument(**random_seed_arg)
    args = parser.parse_args()

    wd = os.path.join(configs.data_dir, "results")
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
            raise ValueError(
                "--quantile must be either 'default' or an integer between 0 and 100"
            )
    if not isinstance(quantile, str) and 0 > quantile > 100:
        raise ValueError(
            "--quantile must be either 'default' or an integer between 0 and 100"
        )

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
            raise ValueError(
                "--batch-ratio must be 'default' or a float between 0 and 1"
            )
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
                f"--n_neigh_rep must be either 'default' or an integer greater than 1 and less than or equal to "
                f"the total number of replicates being clustered across all contexts."
            )
    if not isinstance(n_neigh_rep, str) and n_neigh_rep < 2:
        raise ValueError(
            f"--n_neigh_rep must be either 'default' or an integer greater than 1 and less than or equal to "
            f"the total number of replicates being clustered across all contexts."
        )

    if isinstance(n_neigh_batch, str) and not n_neigh_batch.lower() == "default":
        try:
            n_neigh_batch = int(n_neigh_batch)
        except ValueError:
            raise ValueError(
                f"--n_neigh_batch must be either 'default' or an integer greater than 1 and less than or equal to "
                f"the total number of batches being clustered across all contexts."
            )
    if not isinstance(n_neigh_batch, str) and n_neigh_batch < 2:
        raise ValueError(
            f"--n_neigh_batch must be either 'default' or an integer greater than 1 and less than or equal to "
            f"the total number of batches being clustered across all contexts."
        )

    if isinstance(n_neigh_cont, str) and not n_neigh_cont.lower() == "default":
        try:
            n_neigh_cont = int(n_neigh_cont)
        except ValueError:
            raise ValueError(
                f"--n_neigh_batch must be either 'default' or an integer greater than 1 and less than or equal to "
                f"the total number of batches being clustered across all contexts."
            )
    if not isinstance(n_neigh_cont, str) and n_neigh_cont < 2:
        raise ValueError(
            f"--n_neigh_context must be either 'default' or an integer greater than 1 and less than or equal to "
            f"the total number of contexts being clustered."
        )

    cluster_samples = rpy2_api.Rpy2(
        r_file_path=r_file_path,
        wd=wd,
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
    # cluster_io = rpy2_api.Rpy2(r_file_path=r_file_path)
    # cluster_io_function = cluster_io.call_function("cluster_samples_main")
    # cluster_io_function(
    #     wd,
    #     context_names,
    #     technique,
    #     clust_algo,
    #     label,
    #     min_dist=min_dist,
    #     n_neigh_rep=n_neigh_rep,
    #     n_neigh_batch=n_neigh_batch,
    #     n_neigh_cont=n_neigh_cont,
    #     rep_ratio=rep_ratio,
    #     batch_ratio=batch_ratio,
    #     quantile=quantile,
    #     min_count=min_count,
    #     seed=seed,
    # )


if __name__ == "__main__":
    main(sys.argv[1:])
