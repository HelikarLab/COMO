# This is the Python version of removeLowScoreGenes (used in ftinit) originally written in MATLAB.
# Source code: https://github.com/SysBioChalmers/RAVEN/blob/cd82b5c6ba26a88bed2a9ca179c682a6434fd541/INIT/removeLowScoreGenes.m

# remove_low_score_genes: remove low-scoring genes from model
# This function removes genes from a model based on their scores, a step used by the tINIT package. The function recognizes and differentiates
# between isozymes and subunits of an enzyme complex. Genes are removed from each grRule, subject to the following conditions:
# 1) At least one gene must remain associated with the reaction
# 2) Genes involved in a complex (joined by ANDs) are not removed

## INPUT
# model: Model structure from which genes are to be removed
# gene_score: A vector of scores associated with the model genes. Genes with a positive score will remain in the model, whereas
#             genes with a negative score will try to be removed.
#             If all genes associated with a reaction have a negative score, then the least-negative gene will remain; if there is a tie,
#             one will be selected at random.
#             If a negative-scoring gene is a subunit in a complex, it will not be removed; however, the entire complex may be removed.
#             See the following example cases:
#              Original: G1 or (G2 and G3 and G4)
#              Negative: G1, G2
#              New: G2 and G3 and G4
#
#              Original: G1 or (G2 and G3) or (G4 and G5)
#              Negative: G1, G2
#              New: G4 and G5 [using default complexScoring]
#
#              Original: (G1 and (G2 or G3) and G4)
#              Negative: G2
#              New: (G1 and G3 and G4)
#
# isoenzyme_scoring: Method used to combine the scores of multiple isozymes; 'min','max','median', or 'average'. (optional, default 'max')
# complex_scoring: Method used to combine the scores of enzyme complex subunits: 'min','max', 'median', or 'average'. (optional, default 'min')

## OUTPUT
# new_method: Model with updated gene, gr_rules, and rxn_gene_mat fields after attempting to remove negative-score genes
# rem_genes: A list of negative-score genes that were fully removed from the model. Negative-score genes that were removed from some gr_rules
# but not all will not be included in this list.

import re
from typing import List, Tuple

import numpy as np
from cobra import Model


def remove_low_score_genes(
    model: Model, gene_score: np.ndarray, isoenzyme_scoring: str = "max", complex_scoring: str = "min"
) -> Tuple[Model, List[str]]:
    # Convert logical operators to symbols if needed
    gr_rules = model.gr_rules.copy()
    symbolic_rules = any("&" in rule or "|" in rule for rule in gr_rules)

    if not symbolic_rules:
        gr_rules = [re.sub(r"and", "&", rule) for rule in gr_rules]
        gr_rules = [re.sub(r" or ", "|", rule) for rule in gr_rules]

    # Get unique rules and their mapping
    u_rules, rule_ind = np.unique(gr_rules, return_inverse=True)

    # Process each unique rule
    new_rules = u_rules.copy()
    for i, rule in enumerate(u_rules):
        if not rule or "|" not in rule:
            continue  # Skip empty rules or rules without ORs
        elif "&" in rule:
            new_rules[i] = process_complex_rule(rule, model.genes, gene_score, isoenzyme_scoring, complex_scoring)
        else:
            new_rules[i], _ = process_simple_rule(rule, model.genes, gene_score, isoenzyme_scoring, complex_scoring)

    # Create new model with updated rules
    new_model = model.copy()
    new_model.gr_rules = new_rules[rule_ind]

    # Restore original operator format if needed
    if not symbolic_rules:
        new_model.gr_rules = [re.sub(r"\|", "or", rule) for rule in new_model.gr_rules]
        new_model.gr_rules = [re.sub(r"&", "and", rule) for rule in new_model.gr_rules]

    # Regenerate gene list and rxn_gene_mat
    genes, rxn_gene_mat = get_genes_format_rules(new_model.gr_rules)
    new_model.genes = genes
    new_model.rxn_gene_mat = rxn_gene_mat

    # Find and remove genes no longer present
    rem_ind = ~np.isin(model.genes, new_model.genes)
    rem_genes = model.genes[rem_ind].tolist()

    # Cleanup other gene_related fields
    for field in ["gene_short_names", "proteins", "gene_miriams", "gene_from", "gene_comps"]:
        if hasattr(new_model, field):
            current_field = getattr(new_model, field)
            if isinstance(current_field, list):
                setattr(new_model, field, [val for i, val in enumerate(current_field) if not rem_ind[i]])
    return new_model, rem_genes


def process_simple_rule(rule: str, genes: List[str], g_scores: np.ndarray, isozyme_scoring: str, complex_scoring: str) -> Tuple[str, float]:
    # Helper function to process rules with only ANDs or ORs
    rule_genes = list(set(re.findall(r"[^&|\(\)]+", rule)))
    gene_ind = [genes.index(g) for g in rule_genes]

    if len(rule_genes) < 2:
        return rule, g_scores[gene_ind[0]] if rule_genes else None

    if "&" not in rule:  # Isozyme case
        neg_ind = g_scores[gene_ind] < 0
        if np.any(~neg_ind):
            kept_genes = [g for g, keep in zip(rule_genes, ~neg_ind) if keep]
            updated_rule = "|".join(kept_genes)
            if rule.startswith("("):
                updated_rule = f"({updated_rule})"

            # Recalculate indices for scoring
            rule_genes = kept_genes
            gene_ind = [genes.index(g) for f in rule_genes]
        else:
            # Keep most positive gene
            best_gene = rule_genes[np.argmax(g_scores[gene_ind])]
            updated_rule = best_gene
            gene_ind = [genes.index(best_gene)]
    elif "|" not in rule:  # Complex case
        updated_rule = rule
    else:
        raise ValueError("Mixed AND/OR rules should use process_complex_rule")
    # Score the rule
    valid_scores = g_scores[gene_ind][~np.isnan(g_scores[gene_ind])]
    if len(valid_scores) == 0:
        return updated_rule, None

    if isozyme_scoring.lower() == "max":
        score = np.max(valid_scores)
    elif isozyme_scoring.lower() == "min":
        score = np.min(valid_scores)
    elif isozyme_scoring.lower() == "median":
        score = np.median(valid_scores)
    elif isozyme_scoring.lower() == "average":
        score = np.mean(valid_scores)

    return updated_rule, score


def process_complex_rule(rule: str, genes: List[str], g_scores: np.ndarray, isozyme_scoring: str, complex_scoring: str) -> str:
    """Helper function to process rules containing both AND and OR expressions"""
    # Implement similar logic as MATLAB version
    # This would involve recursive processing of nested rules
    raise NotImplementedError("Complex rule processing not implemented")


def get_genes_from_gr_rules(gr_rules: List[str]) -> Tuple[List[str], np.ndarray]:
    """Extract unique genes and create reaction-gene matrix"""
    # Extract all unique genes from rules
    all_genes = set()
    for rule in gr_rules:
        if rule:
            all_genes.update(re.findall(r"[^&|\(\)]+", rule))

    # Create binary matrix
    genes = sorted(all_genes)
    rxn_gene_mat = np.zeros((len(gr_rules), len(genes)), dtype=bool)

    for i, rule in enumerate(gr_rules):
        if rule:
            rule_genes = re.findall(r"[^&|\(\)]+", rule)
            for g in rule_genes:
                rxn_gene_mat[i, genes.index(g)] = True

    return genes, rxn_gene_mat
