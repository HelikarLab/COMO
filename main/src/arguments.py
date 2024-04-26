"""
The purpose of this file is to create a uniform interface for the arguments that are passed to python files
"""

import argparse
import math
from typing import Callable, Union


def ranged_type(value_type, min_value, max_value) -> Callable:
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
    "flag": "--taxon-id",
    "type": str,
    "required": True,
    "dest": "taxon_id",
    "help": "The taxon ID of the organism, such as 9096 or 'Mus Musculus'",
}

context_names_arg = {
    "flag": "--context-names",
    "type": str,
    "required": True,
    "dest": "context_names",
    "help": "Tissue/cell name of models to generate. If making multiple models in a batch, then use the format: 'context1 context2 context3' ",
}

filtering_technique_arg = {
    "flag": "--filtering-technique",
    "type": str,
    "choices": ["zfpkm", "quantile", "cpm"],
    "required": True,
    "dest": "filtering_technique",
    "help": "The filtering technique to use. Options are: 'zfpkm', 'quantile', or 'cpm'",
}

cluster_algorithm_arg = {
    "flag": "--cluster-algorithm",
    "type": str,
    "choices": ["mca", "umap"],
    "required": False,
    "default": "umap",
    "dest": "clust_algo",
    "help": "The clustering algorithm to use. Options are: 'mca' or 'umap'.",
}

label_arg = {
    "flag": "--no-label",
    "type": bool,
    "action": "store_false",
    "default": "True",
    "required": False,
    "dest": "label",
    "help": "Do not label replicates/batches/context names on the plot. May be ugly for large sets.",
}

min_dist_arg = {
    "flag": "--min-dist",
    "type": ranged_type(float, 0.0, 1.0),
    "required": False,
    "default": 0.01,
    "dest": "min_dist",
    "help": "The minimum distance for UMAP clustering. Must be between 0 and 1.",
}

replicate_ratio_arg = {
    "flag": "--replicate-ratio",
    "type": ranged_type(float, 0.0, 1.0),
    "required": False,
    "default": 0.9,
    "dest": "replicate_ratio",
    "help": "The ratio of genes that must be active in a replicate for the associated batch to be considered active. Must be between 0 and 1 (inclusive).",
}

batch_ratio_arg = {
    "flag": "--batch-ratio",
    "type": ranged_type(float, 0.0, 1.0),
    "required": False,
    "default": 0.9,
    "dest": "batch_ratio",
    "help": "The ratio of genes that must be active in a batch for the associated context to be considered active. Must be between 0 and 1 (inclusive).",
}

num_neighbors_replicate_arg = {
    "flag": "--num-neighbors-replicate",
    "type": str | int,
    "required": False,
    "default": "default",
    "dest": "n_neigh_rep",
    "help": "N-nearest neighbors for replicate clustering. 'default' uses the total number of replicates. Must be 'default' or a number less than the total number of replicates.",
}

num_neighbors_batch_arg = {
    "flag": "--num-neighbors-batch",
    "type": str | int,
    "required": False,
    "default": "default",
    "dest": "n_neigh_batch",
    "help": "N-nearest neighbors for batch clustering. 'default' uses the total number of batches. Must be 'default' or a number less than the total number of batches.",
}

num_neighbors_context_arg = {
    "flag": "--num-neighbors-context",
    "type": str | int,
    "required": False,
    "default": "default",
    "dest": "n_neigh_cont",
    "help": "N-nearest neighbors for context clustering. 'default' uses the total number of contexts. Must be 'default' or a number less than the total number of contexts.",
}

min_count_arg = {
    "flag": "--min-count",
    "type": int | str,
    "required": False,
    "default": "default",
    "dest": "min_count",
    "help": "The minimum number of cells that must express a gene for it to be considered active. If 'default' is used, the minimum count is set to 10e6.",
}

quantile_arg = {
    "flag": "--quantile",
    "type": ranged_type(float, 0.0, 1.0),
    "required": False,
    "default": 0.5,
    "dest": "quantile",
    "help": "The quantile to use for filtering. Must be between 0 and 1 (inclusive).",
}

random_seed_arg = {
    "flag": "--seed",
    "type": int,
    "required": False,
    "default": -1,
    "dest": "random_seed",
    "help": "The random seed to use for clustering initialization. If -1 is used, the seed is randomly generated.",
}

reference_model_filepath_arg = {
    "flag": "--reference-model-filepath",
    "type": str,
    "required": True,
    "dest": "modelfile",
    "help": "Name of Genome-scale metabolic model to seed the context model to. For example, the GeneralModelUpdatedV2.mat, is a modified Recon3D model. We also provide iMM_madrid for mouse. The file type can be in .mat, .xml, or .json format.",
}

active_genes_filepath_arg = {
    "flag": "--active-genes-filepath",
    "type": str,
    "required": True,
    "dest": "genefile",
    "help": "Path to logical table of active genes output from merge_xomics.py called ActiveGenes_contextName_Merged.csv. Should be in the corresponding context/context folder inside /main/data/results/contextName/. The json file output from the function using the context of interest as the key can be used here.",
}

objective_function_arg = {
    "flag": "--objective",
    "type": str,
    "required": False,
    "default": "biomass_reeaction",
    "dest": "objective",
    "help": "The reaction ID of the objective function to use. Generally a biomass function (biomass_maintenance, biomass_maintenance_noTrTr, etc.).",
}

boundary_reactions_filepath_arg = {
    "flag": "--boundary-reactions-filepath",
    "type": str,
    "required": False,
    "default": None,
    "dest": "boundary_rxns_filepath",
    "help": "Path to file contains the exchange (media), sink, and demand reactions which the model should use to fulfill the reactions governed by transcriptomic and proteomics data inputs. It must be a csv or xlsx with three columns: Rxn, Lowerbound, Upperbound. If not specified, COMO will allow ALL BOUNDARY REACTIONS THAT ARE OPEN IN THE REFERENCE MODEL TO BE USED!",
}

exclude_reactions_filepath_arg = {
    "flag": "--exclude-reactions-filepath",
    "type": str,
    "required": False,
    "default": None,
    "dest": "exclude_rxns_filepath",
    "help": "Filepath to file that contains reactions which will be removed from active reactions the model to use when seeding, even if considered active from transcriptomic and proteomics data inputs. It must be a csv or xlsx with one column of reaction IDs consistent with the reference model",
}

force_reactions_filepath_arg = {
    "flag": "--force-reactions-filepath",
    "type": str,
    "required": False,
    "default": None,
    "dest": "force_rxns_filepath",
    "help": "Filepath to file that contains reactions which will be added to active reactions for the model to use when seeding (unless it causes the model to be unsolvable), regardless of results of transcriptomic and proteomics data inputs. It must be a csv or xlsx with one column of reaction IDs consistent with the reference model",
}

reconstruction_algorithm_arg = {
    "flag": "--algorithm",
    "type": str,
    "required": False,
    "default": "GIMME",
    "choices": ["GIMME", "FASTCORE", "iMAT", "tINIT"],
    "dest": "algorithm",
    "help": "Algorithm used to seed context specific model to the Genome-scale model. Can be either GIMME, FASTCORE, iMAT, or tINIT.",
}

imat_low_threshold_arg = {
    "flag": "--low-threshold",
    "type": ranged_type(int, -math.inf, math.inf),
    "required": False,
    "default": -5,
    "dest": "low_threshold",
    "help": "Lower threshold for iMAT algorithm",
}

imat_high_threshold_arg = {
    "flag": "--high-threshold",
    "type": ranged_type(int, -math.inf, math.inf),
    "required": False,
    "default": -3,
    "dest": "high_threshold",
    "help": "Middle to high bin cutoff for iMAT algorithm.",
}

reconstruction_solver_arg = {
    "flag": "--solver",
    "type": str,
    "required": False,
    "default": "glpk",
    "dest": "solver",
    "help": "Solver used to seed model and attempt to solve objective. Default is GLPK, also takes GUROBI but you must mount a container license to the Docker to use. An academic license can be obtained for free. See the README on the Github or Dockerhub for information on mounting this license.",
}

output_filetypes_arg = {
    "flag": "--output-filetypes",
    "type": str,
    "required": False,
    "default": "mat",
    "dest": "output_filetypes",
    "help": "Filetypes to save seeded model type. Can be either a string with one filetype such as 'xml' or multiple in the format 'extension1 extension2 extension3'. If you want to output in all 3 accepted formats,  would be: 'mat xml json'.",
}

config_file_arg = {
    "flag": "--config-file",
    "type": str,
    "required": True,
    "dest": "config_file",
    "help": "The path to the configuration file",
}

data_source_arg = {
    "flag": "--data-source",
    "type": str,
    "choices": ["rnaseq", "microarray"],
    "required": True,
    "dest": "data_source",
    "help": "Source of data being used, either rnaseq or microarray",
}

context_model_filepath_arg = {
    "flag": "--context-model",
    "type": str,
    "required": True,
    "dest": "model",
    "help": "The context-specific model file. Must end in .mat, .xml, or .json.",
}

disease_names_arg = {
    "flag": "--disease-names",
    "type": str,
    "required": True,
    "dest": "disease",
    "help": "The disease names of the organism, such as 'cancer' or 'diabetes'. Only accepts a single disease at a time.",
}

disease_up_filepath_arg = {
    "flag": "--disease-up",
    "type": str,
    "required": True,
    "dest": "disease_up",
    "help": "The upregulated gene filepath for the disease.",
}

disease_down_filepath_arg = {
    "flag": "--disease-down",
    "type": str,
    "required": True,
    "dest": "disease_down",
    "help": "The downregulated gene filepath for the disease.",
}

raw_drug_filepath_arg = {
    "flag": "--raw-drug-file",
    "type": str,
    "required": True,
    "dest": "raw_drug_file",
    "help": "The raw drug file path.",
}

reference_flux_filepath_arg = {
    "flag": "--reference-flux-file",
    "type": str,
    "required": False,
    "dest": "ref_flux_file",
    "help": "The reference flux file path.",
}

test_all_genes_arg = {
    "flag": "--test-all",
    "action": "store_true",
    "required": False,
    "default": False,
    "dest": "test_all",
    "help": "Test all genes in the model, even ones predicted to have little-to-no effect.",
}

parsimonious_fba_arg = {
    "flag": "--parsimonious",
    "action": "store_true",
    "required": False,
    "default": False,
    "dest": "pars_flag",
    "help": "Use parsimonious FBA to calculate the optimal reference solution. Only valid if a reference flux file is not being provided.",
}

merge_distribution_arg = {
    "flag": "--merge-distribution",
    "action": "store_true",
    "required": False,
    "default": False,
    "dest": "merge_distro",
    "help": "Flag to merge zFPKM distributions. Required if using iMAT reconstruction algorithm in create_context_specific_model.py. Must have run rnaseq_gen.py with 'zFPKM' as '--technique'. If --proteomics-config-file is given will merge proteomics distributions with zFPKM distributions using a weighted scheme.",
}

keep_gene_score_arg = {
    "flag": "--keep-gene-scores",
    "action": "store_true",
    "required": False,
    "default": True,
    "dest": "keep_gene_score",
    "help": "When merging z-score distributions of expression, if using both protein abundance and transcipt zFPKM flag true if you wish to keep z-score of genes with no protein data, flag false if you wish to discard and treat as no expression",
}

microarray_config_filename_arg = {
    "flag": "--microarray-config-file",
    "type": str,
    "required": False,
    "dest": "microarray_file",
    "help": "The name of the microarray configuration file",
}

total_rnaseq_filename_arg = {
    "flag": "--total-rnaseq-config-file",
    "type": str,
    "required": False,
    "dest": "trnaseq_file",
    "help": "The name of the total RNA-seq file",
}

mrnaseq_filename_arg = {
    "flag": "--mrnaseq-config-file",
    "type": str,
    "required": False,
    "dest": "mrnaseq_file",
    "help": "The name of the mRNA RNA-seq file",
}

scrnaseq_filename_arg = {
    "flag": "--scrnaseq-config-file",
    "type": str,
    "required": False,
    "dest": "scrnaseq_file",
    "help": "The name of the single-cell RNA-seq file",
}

proteomics_config_filename_arg = {
    "flag": "--proteomics-config-file",
    "type": str,
    "required": False,
    "dest": "proteomics_file",
    "help": "The name of the proteomics configuration file",
}

expression_requirement_arg = {
    "flag": "--expression-requirement",
    "type": int | str,
    "required": False,
    "default": "default",
    "dest": "expression_requirement",
    "help": "Number of sources with active gene for it to be considered active even if it is not a high confidence-gene. Set to 'default' to set expression requirement to the total number of data sources provided.",
}

requirement_adjustment_arg = {
    "flag": "--requirement-adjust",
    "type": str,
    "choices": ["progressive", "regressive", "flat", "custom"],
    "required": False,
    "default": "flat",
    "dest": "requirement_adjust",
    "help": "Technique to adjust expression requirement based on differences in number of provided data source types.",
}

custom_requirement_file_arg = {
    "flag": "--custom-requirement-file",
    "type": str,
    "required": False,
    "default": "SKIP",
    "dest": "custom_file",
    "help": "Name of .xlsx file where the first column is a context name and the second column is the expression requirement for that context",
}

no_high_confidence_genes_arg = {
    "flag": "--no-hc",
    "action": "store_true",
    "required": False,
    "default": False,
    "dest": "no_hc",
    "help": "Use this to prevent high-confidence genes forcing a gene to be used in final model irrespective of other other data sources",
}

no_na_adjustment_arg = {
    "flag": "--no-na-adjust",
    "action": "store_true",
    "required": False,
    "default": False,
    "dest": "no_na",
    "help": "Use this to prevent genes missing in a data source library, but present in others from subtracting 1 from the expression requirement per data source that gene is missing in",
}

total_rnaseq_weight_arg = {
    "flag": "--total-rnaseq-weight",
    "type": ranged_type(int, 0, math.inf),
    "required": False,
    "default": 1,
    "dest": "tweight",
    "help": "Total RNA-seq weight for merging zFPKM distribution",
}

mrnaseq_weight_arg = {
    "flag": "--mrnaseq-weight",
    "type": ranged_type(int, 0, math.inf),
    "required": False,
    "default": 1,
    "dest": "mweight",
    "help": "mRNA RNA-seq weight for merging zFPKM distribution",
}

scrnaseq_weight_arg = {
    "flag": "--single-cell-rnaseq-weight",
    "type": ranged_type(int, 0, math.inf),
    "required": False,
    "default": 1,
    "dest": "sweight",
    "help": "Single-cell RNA-seq weight for merging zFPKM distribution",
}

proteomics_weight_arg = {
    "flag": "--protein-weight",
    "type": ranged_type(int, 0, math.inf),
    "required": False,
    "default": 2,
    "dest": "pweight",
    "help": "Proteomics weight for merging zFPKM distribution",
}

high_replicate_ratio_arg = {
    "flag": "--high-replicate-ratio",
    "type": ranged_type(float, 0.0, 1.0),
    "required": False,
    "default": 0.5,
    "dest": "hi_replicate_ratio",
    "help": "Ratio of replicates required for a gene to be considered high-confidence in that group",
}

high_batch_ratio_arg = {
    "flag": "--high-batch-ratio",
    "type": ranged_type(float, 0.0, 1.0),
    "required": False,
    "default": 0.5,
    "dest": "hi_batch_ratio",
    "help": "Ratio of groups (batches or studies) required for a gene to be considered high-confidence in a context",
}

min_zfpkm_arg = {
    "flag": "--min-zfpkm",
    "type": int,
    "required": False,
    "default": -3,
    "dest": "min_zfpkm",
    "help": "Cutoff used for zfpkm. Minimum zfpkm to be considered expressed. According to PMC3870982, should be between -3 and -2.",
}

library_prep_arg = {
    "flag": "--library-prep",
    "type": str,
    "choices": ["total", "mrna", "scrna"],
    "required": True,
    "dest": "prep",
    "help": "Library preparation method used for the data.",
}

# fmt: off
gene_format_arg = {
    "flag": "--gene-format",
    "type": str,
    "required": False,
    "default": "Ensembl Gene ID",
    # "choices": [
    #     "ENSEMBL", "ENSEMBLE", "ENSG", "ENSMUSG", "ENSEMBL ID", "ENSEMBL GENE ID",
    #     "HGNC SYMBOL", "HUGO", "HUGO SYMBOL", "SYMBOL", "HGNC",
    #     "GENE SYMBOL", "ENTREZ", "ENTRES", "ENTREZ ID", "ENTREZ NUMBER" "GENE ID",
    # ],
    "dest": "gene_format",
    "help": "The format of the gene IDs in the data. Must be one of 'entrez', 'symbol', or 'ensembl'.",
}
# fmt: on

provide_count_matrix_arg = {
    "flag": "--provide-matrix",
    "action": "store_true",
    "required": False,
    "default": False,
    "dest": "provide_matrix",
    "help": "Provide your own count matrix. Requires additional argument '--matrix' which is .csv file where colnames are sample names (in contextName_SXRY format) and rownames are genes in in format specified by --gene-format",
}

create_matrix_arg = {
    "flag": "--create-matrix",
    "action": "store_true",
    "required": False,
    "default": False,
    "dest": "make_matrix",
    "help": "Flag for if you want to make a counts matrix, but not a config file. Requires a correctly structured COMO_input folder in /work/data/. Can make one using: https://github.com/HelikarLab/FastqToGeneCounts",
}

provided_matrix_filename_arg = {
    "flag": "--matrix",
    "type": str,
    "required": False,
    "default": "SKIP",
    "dest": "provided_matrix_fname",
    "help": "The file name of the provided count matrix",
}
