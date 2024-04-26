"""
The purpose of this file is to create a uniform interface for the arguments that are passed to python files
"""

import argparse
import math
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

reference_model_filepath_arg = {
    "name_or_flags": ["--reference-model-filepath"],
    "type": str,
    "required": True,
    "dest": "modelfile",
    "help": "Name of Genome-scale metabolic model to seed the context model to. For example, the GeneralModelUpdatedV2.mat, is a modified Recon3D model. We also provide iMM_madrid for mouse. The file type can be in .mat, .xml, or .json format.",
}

active_genes_filepath_arg = {
    "name_or_flags": ["--active-genes-filepath"],
    "type": str,
    "required": True,
    "dest": "genefile",
    "help": "Path to logical table of active genes output from merge_xomics.py called ActiveGenes_contextName_Merged.csv. Should be in the corresponding context/context folder inside /main/data/results/contextName/. The json file output from the function using the context of interest as the key can be used here.",
}

objective_function_arg = {
    "name_or_flags": ["--objective"],
    "type": str,
    "required": False,
    "default": "biomass_reeaction",
    "dest": "objective",
    "help": "The reaction ID of the objective function to use. Generally a biomass function (biomass_maintenance, biomass_maintenance_noTrTr, etc.).",
}

boundary_reactions_filepath_arg = {
    "name_or_flags": ["--boundary-reactions-filepath"],
    "type": str,
    "required": False,
    "default": None,
    "dest": "boundary_rxns_filepath",
    "help": "Path to file contains the exchange (media), sink, and demand reactions which the model should use to fulfill the reactions governed by transcriptomic and proteomics data inputs. It must be a csv or xlsx with three columns: Rxn, Lowerbound, Upperbound. If not specified, COMO will allow ALL BOUNDARY REACTIONS THAT ARE OPEN IN THE REFERENCE MODEL TO BE USED!",
}

exclude_reactions_filepath_arg = {
    "name_or_flags": ["--exclude-reactions-filepath"],
    "type": str,
    "required": False,
    "default": None,
    "dest": "exclude_rxns_filepath",
    "help": "Filepath to file that contains reactions which will be removed from active reactions the model to use when seeding, even if considered active from transcriptomic and proteomics data inputs. It must be a csv or xlsx with one column of reaction IDs consistent with the reference model",
}

force_reactions_filepath_arg = {
    "name_or_flags": ["--force-reactions-filepath"],
    "type": str,
    "required": False,
    "default": None,
    "dest": "force_rxns_filepath",
    "help": "Filepath to file that contains reactions which will be added to active reactions for the model to use when seeding (unless it causes the model to be unsolvable), regardless of results of transcriptomic and proteomics data inputs. It must be a csv or xlsx with one column of reaction IDs consistent with the reference model",
}

reconstruction_algorithm_arg = {
    "name_or_flags": ["--algorithm"],
    "type": str,
    "required": False,
    "default": "GIMME",
    "choices": ["GIMME", "FASTCORE", "iMAT", "tINIT"],
    "dest": "algorithm",
    "help": "Algorithm used to seed context specific model to the Genome-scale model. Can be either GIMME, FASTCORE, iMAT, or tINIT.",
}

imat_low_threshold_arg = {
    "name_or_flags": ["--low-threshold"],
    "type": ranged_type(int, -math.inf, math.inf),
    "required": False,
    "default": -5,
    "dest": "low_threshold",
    "help": "Lower threshold for iMAT algorithm",
}

imat_high_threshold_arg = {
    "name_or_flags": ["--high-threshold"],
    "type": ranged_type(int, -math.inf, math.inf),
    "required": False,
    "default": -3,
    "dest": "high_threshold",
    "help": "Middle to high bin cutoff for iMAT algorithm.",
}

reconstruction_solver_arg = {
    "name_or_flags": ["--solver"],
    "type": str,
    "required": False,
    "default": "glpk",
    "dest": "solver",
    "help": "Solver used to seed model and attempt to solve objective. Default is GLPK, also takes GUROBI but you must mount a container license to the Docker to use. An academic license can be obtained for free. See the README on the Github or Dockerhub for information on mounting this license.",
}

output_filetypes_arg = {
    "name_or_flags": ["--output-filetypes"],
    "type": str,
    "required": False,
    "default": "mat",
    "dest": "output_filetypes",
    "help": "Filetypes to save seeded model type. Can be either a string with one filetype such as 'xml' or multiple in the format 'extension1 extension2 extension3'. If you want to output in all 3 accepted formats,  would be: 'mat xml json'.",
}
