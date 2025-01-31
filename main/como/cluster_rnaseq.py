from __future__ import annotations

import argparse
from curses.ascii import isdigit
from dataclasses import dataclass
from typing import Any

import numpy as np

from como.data_types import LogLevel
from como.utils import _log_and_raise_error, stringlist_to_list


@dataclass
class _Arguments:
    context_names: list[str]
    filtering_technique: str
    cluster_algorithm: str
    label_plot: str
    min_distance: float
    replicate_ratio: Any
    batch_ratio: Any
    num_replicate_neighbors: Any
    num_batch_neighbors: Any
    num_context_neighbors: Any
    min_active_count: int | str
    quantile: Any
    seed: int

    def __post_init__(self):  # noqa: C901, ignore too complex
        self.filtering_technique = self.filtering_technique.lower()
        self.cluster_algorithm = self.cluster_algorithm.lower()

        if self.seed == -1:
            self.seed = np.random.randint(0, 100_000)

        if (isdigit(self.min_active_count) and int(self.min_active_count) < 0) or self.min_active_count != "default":
            _log_and_raise_error(
                "min_active_count must be either 'default' or an integer > 0",
                error=ValueError,
                level=LogLevel.ERROR,
            )

        if (isdigit(self.quantile) and 0 > int(self.quantile) > 100) or self.quantile != "default":
            _log_and_raise_error(
                "quantile must be either 'default' or an integer between 0 and 100",
                error=ValueError,
                level=LogLevel.ERROR,
            )

        if (isdigit(self.replicate_ratio) and 0 > self.replicate_ratio > 1.0) or self.replicate_ratio != "default":
            _log_and_raise_error(
                "--rep-ratio must be either 'default' or a float between 0 and 1",
                error=ValueError,
                level=LogLevel.ERROR,
            )

        if (isdigit(self.batch_ratio) and 0 > self.batch_ratio > 1.0) or self.batch_ratio != "default":
            _log_and_raise_error(
                "--batch-ratio must be either 'default' or a float between 0 and 1",
                error=ValueError,
                level=LogLevel.ERROR,
            )

        if self.filtering_technique.lower() not in {"quantile", "tpm", "cpm", "zfpkm"}:
            _log_and_raise_error(
                "--technique must be either 'quantile', 'tpm', 'cpm', 'zfpkm'",
                error=ValueError,
                level=LogLevel.ERROR,
            )

        if self.filtering_technique.lower() == "tpm":
            self.filtering_technique = "quantile"

        if self.cluster_algorithm.lower() not in {"mca", "umap"}:
            _log_and_raise_error(
                "--clust_algo must be either 'mca', 'umap'",
                error=ValueError,
                level=LogLevel.ERROR,
            )

        if 0 > self.min_distance > 1.0:
            _log_and_raise_error(
                "--min_dist must be a float between 0 and 1",
                error=ValueError,
                level=LogLevel.ERROR,
            )

        if (
            isdigit(self.num_replicate_neighbors) and self.num_replicate_neighbors < 1
        ) or self.num_replicate_neighbors != "default":
            _log_and_raise_error(
                "--n-neighbors-rep must be either 'default' or an integer > 1",
                error=ValueError,
                level=LogLevel.ERROR,
            )

        if (
            isdigit(self.num_batch_neighbors) and self.num_batch_neighbors < 1
        ) or self.num_batch_neighbors != "default":
            _log_and_raise_error(
                "--n-neighbors-batch must be either 'default' or an integer > 1",
                error=ValueError,
                level=LogLevel.ERROR,
            )

        if (
            isdigit(self.num_context_neighbors) and self.num_context_neighbors < 1
        ) or self.num_context_neighbors != "default":
            _log_and_raise_error(
                "--n-neighbors-context must be either 'default' or an integer > 1",
                error=ValueError,
                level=LogLevel.ERROR,
            )


def _parse_args() -> _Arguments:
    parser = argparse.ArgumentParser(
        prog="cluster_rnaseq.py",
        description="Cluster RNA-seq Data using Multiple Correspondence Analysis or UMAP. Clusters at the replicate, "
        "batch/study, and context levels.",
        epilog="For additional help, please post questions/issues in the MADRID GitHub repo at "
        "https://github.com/HelikarLab/MADRID or email babessell@gmail.com",
    )
    parser.add_argument(
        "--context-names",
        type=str,
        required=True,
        dest="context_names",
        help="Tissue/cell name of models to generate",
    )
    parser.add_argument(
        "--filtering-technique",
        type=str,
        required=True,
        dest="filtering_technique",
        help="'zfpkm', 'quantile', or 'cpm'",
    )
    parser.add_argument(
        "--cluster-algorithm",
        type=str,
        required=False,
        default="umap",
        dest="cluster_algorithm",
        help="Clustering algorithm to use. 'mca' or 'umap'.",
    )
    parser.add_argument(
        "--label-plot",
        type=str,
        required=False,
        default=True,
        dest="label_plot",
        help="Set to True to label replicate/batch/context names on the plots. May be ugly for large sets",
    )
    parser.add_argument(
        "--min-distance",
        type=float,
        required=False,
        default=0.01,
        dest="min_distance",
        help="Minimum distance for UMAP clustering. Must be between 0 and 1",
    )
    parser.add_argument(
        "-r",
        "--replicate-ratio",
        type=str,
        required=False,
        default=0.9,
        dest="replicate_ratio",
        help="Ratio of genes active in replicates for a batch/study to be active",
    )
    parser.add_argument(
        "-b",
        "--batch-ratio",
        type=str or float,
        required=False,
        default=0.9,
        dest="batch_ratio",
        help="Ratio of genes active in a batch/study to be active in the context",
    )
    parser.add_argument(
        "--num-replicate-neighbors",
        type=str or float,
        required=False,
        default="default",
        dest="num_replicate_neighbors",
        help="N nearest neighbors for replicate clustering, 'default' is total number of replicates",
    )
    parser.add_argument(
        "-nb",
        "--num-batch-neighbors",
        type=str or float,
        required=False,
        default="default",
        dest="num_batch_neighbors",
        help="N nearest neighbors for batch clustering, 'default' is total number of batches",
    )
    parser.add_argument(
        "--num-context-neighbors",
        type=str or float,
        required=False,
        default="default",
        dest="num_context_neighbors",
        help="N nearest neighbors for context clustering, 'default' is total number of contexts",
    )
    parser.add_argument(
        "--min-active-count",
        type=str or int,
        required=False,
        default="default",
        dest="min_active_count",
        help="Ratio of active genes in a batch/study to be active in the context",
    )
    parser.add_argument(
        "--quantile",
        type=str or int,
        required=False,
        default=0.5,
        dest="quantile",
        help="Ratio of active genes in a batch/study to be active in the context",
    )
    parser.add_argument(
        "--seed",
        type=int,
        required=False,
        default=-1,
        dest="seed",
        help="Random seed for clustering algorithm initialization",
    )
    args = parser.parse_args()
    args.context_names = stringlist_to_list(args.context_names)
    return _Arguments(**vars(args))
