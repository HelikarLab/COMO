"""
The purpose of this file is to create a uniform interface for the arguments that are passed to python files
"""

import argparse
import typing


def ranged_type(value_type, min_value, max_value) -> typing.Callable:
    """
    Return function handle of an argument type function for ArgumentParser checking a range

    From: stackoverflow.com/a/71112312

    Args:
        value_type (float | int): Value-type to convert the arg to
        min_value (float | int): Minimum acceptable value
        max_value (float | int): Maximum acceptable argument

    Returns:
        float | int: Function handle of an argument type function for ArgumentParser
    """

    def range_checker(arg: str):
        try:
            f = value_type(arg)
        except ValueError:
            raise argparse.ArgumentTypeError(
                f"Must be a valid {value_type}, received {arg}"
            )
        if f < min_value or f > max_value:
            raise argparse.ArgumentTypeError(
                f"Must be in range [{min_value}, {max_value}], received {f}"
            )
        return f

    return range_checker


taxon_id_arg = {
    "name_or_flags": ["--taxon_id"],
    "type": int | str,
    "required": True,
    "dest": "taxon_id",
    "help": "The taxon ID of the organism, such as 9096 or 'Mus Musculus'",
}

context_names_arg = {
    "name_or_flags": ["--context-names"],
    "type": str,
    "required": True,
    "dest": "context_names",
    "help": """Tissue/cell name of models to generate. If making multiple models in a batch, then use the format: "context1 context2 context3" """,
}

filtering_technique_arg = {
    "name_or_flags": ["--filtering-technique"],
    "type": str,
    "choices": ["zfpkm", "quantile", "cpm"],
    "required": True,
    "dest": "filtering_technique",
    "help": "The filtering technique to use. Options are: 'zfpkm', 'quantile', or 'cpm'",
}

cluster_algorithm_arg = {
    "name_or_flags": ["--cluster-algorithm"],
    "type": str,
    "choices": ["mca", "umap"],
    "required": False,
    "default": "umap",
    "dest": "clust_algo",
    "help": "The clustering algorithm to use. Options are: 'mca' or 'umap'.",
}

label_arg = {
    "name_or_flags": ["--no-label"],
    "type": bool,
    "action": "store_false",
    "default": "True",
    "required": False,
    "dest": "label",
    "help": "Do not label replicates/batches/context names on the plot. May be ugly for large sets.",
}

min_dist_arg = {
    "name_or_flags": ["--min-dist"],
    "type": ranged_type(float, 0.0, 1.0),
    "required": False,
    "default": 0.01,
    "dest": "min_dist",
    "help": "The minimum distance for UMAP clustering. Must be between 0 and 1.",
}

replicate_ratio_arg = {
    "name_or_flags": ["--replicate-ratio"],
    "type": ranged_type(float, 0.0, 1.0),
    "required": False,
    "default": 0.9,
    "dest": "replicate_ratio",
    "help": "The ratio of genes that must be active in a replicate for the associated batch to be considered active. Must be between 0 and 1 (inclusive).",
}

batch_ratio_arg = {
    "name_or_flags": ["--batch-ratio"],
    "type": ranged_type(float, 0.0, 1.0),
    "required": False,
    "default": 0.9,
    "dest": "batch_ratio",
    "help": "The ratio of genes that must be active in a batch for the associated context to be considered active. Must be between 0 and 1 (inclusive).",
}

num_neighbors_replicate_arg = {
    "name_or_flags": ["--num-neighbors-replicate"],
    "type": str | int,
    "required": False,
    "default": "default",
    "dest": "n_neigh_rep",
    "help": "N-nearest neighbors for replicate clustering. 'default' uses the total number of replicates. Must be 'default' or a number less than the total number of replicates.",
}

num_neighbors_batch_arg = {
    "name_or_flags": ["--num-neighbors-batch"],
    "type": str | int,
    "required": False,
    "default": "default",
    "dest": "n_neigh_batch",
    "help": "N-nearest neighbors for batch clustering. 'default' uses the total number of batches. Must be 'default' or a number less than the total number of batches.",
}

num_neighbors_context_arg = {
    "name_or_flags": ["--num-neighbors-context"],
    "type": str | int,
    "required": False,
    "default": "default",
    "dest": "n_neigh_cont",
    "help": "N-nearest neighbors for context clustering. 'default' uses the total number of contexts. Must be 'default' or a number less than the total number of contexts.",
}

min_count_arg = {
    "name_or_flags": ["--min-count"],
    "type": int | str,
    "required": False,
    "default": "default",
    "dest": "min_count",
    "help": "The minimum number of cells that must express a gene for it to be considered active. If 'default' is used, the minimum count is set to 10e6.",
}

quantile_arg = {
    "name_or_flags": ["--quantile"],
    "type": ranged_type(float, 0.0, 1.0),
    "required": False,
    "default": 0.5,
    "dest": "quantile",
    "help": "The quantile to use for filtering. Must be between 0 and 1 (inclusive).",
}

random_seed_arg = {
    "name_or_flags": ["--seed"],
    "type": int,
    "required": False,
    "default": -1,
    "dest": "random_seed",
    "help": "The random seed to use for clustering initialization. If -1 is used, the seed is randomly generated.",
}
