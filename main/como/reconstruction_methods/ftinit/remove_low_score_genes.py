# This is the Python version of removeLowScoreGenes (used in ftinit) originally written in MATLAB.
# Source code: https://github.com/SysBioChalmers/RAVEN/blob/cd82b5c6ba26a88bed2a9ca179c682a6434fd541/INIT/removeLowScoreGenes.m
from copy import deepcopy


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
import random
import numpy as np
from copy import deepcopy
from cobra import Model

def remove_low_score_genes(model, gene_scores, isoenzyme_scoring='max', complex_scoring='min'):
    if len(model.genes) != len(gene_scores):
        raise ValueError("The dimentions of the model genes and gene_scores do not match")

    gene_scores = np.array(gene_scores, dtype=float)
    new_model = deepcopy(model)

    # Convert logical operators to symbolic form
    gr_rules = [gpr.rule_str if hasattr(gpr, 'rule_str') else "" for gpr in new_model.reactions]
    symbolic_rules = any('&' in rule or '|' in rule for rule in gr_rules)
    gr_rules = [re.sub(r' and ', ' & ', rule) for rule in gr_rules]
    gr_rules = [re.sub(r' or ', ' | ', rule) for rule in gr_rules]

    # Get the unique rules
    unique_rules = list(dict.fromkeys(gr_rules)) # preserves order
    rule_map = {rule: i for i, rule in enumerate(unique_rules)}

    # Process each unique rule
    new_rules = unique_rules.copy()
    for i, rule in enumerate(unique_rules):
        if not rule or '|' not in rule:
            continue
        elif '&' in rule:
            new_rules[i] = _process_complex_rule(rule, model.genes, gene_scores, isoenzyme_scoring, complex_scoring)
        else:
            new_rules[i] = _process_simple_rule(rule, model.genes, gene_scores, isoenzyme_scoring,complex_scoring)[0]

    # Re-map updated rules back to model reactions
    updated_gr_rules = [new_rules[rule_map[rule]] for rule in gr_rules]

    # Restore original logical formatting if needed
    if not symbolic_rules:
        updated_gr_rules = [re.sub(r' \| ', ' or ', rule) for rule in updated_gr_rules]
        updated_gr_rules = [re.sub(r' & ', ' and ', rule) for rule in updated_gr_rules]

        # Updated gene-reactions associations
        for rxn, new_rule in zip(new_model.reactions, updated_gr_rules):
            rxn.gene_reaction_rule = new_rule

        # Determine removed genes
        current_genes = [g.id for g in new_model.genes]
        removed_genes = [g.id for g in model.genes if g.id not in current_genes]
    else:
        removed_genes = []
    return new_model, removed_genes

def _process_simple_rule(rule, genes, g_scores, isozyme_scoring, complex_scoring):
    """Handle rules that are all ORs (isozymes) or all ANDs (complexes)"""
    rule_genes = list(set(re.findall(r'[^&|\(\) ]+', rule)))
    gene_indices = [genes.index(g) for g in rule_genes if g in genes]

    if len(rule_genes) < 2:
        r_score = g_scores[gene_indices[0]] if gene_indices else np.nan
        return rule. r_score
    if '&' not in rule:
        neg_idx = g_scores[gene_indices] < 0
        if all(neg_idx):
            # Keep least negative, break ties randomly
            vals = g_scores[gene_indices] + np.random.rand(len(gene_indices)) * 1e-8
            keep_idx = np.argmax(vals)
            updated_rule = rule_genes[keep_idx]
        elif np.sum(~neg_idx) == 1:
            updated_rule = rule_genes[np.where(~neg_idx)[0][0]]
        else:
            updated_rule = ' | '.join([rule_genes[i] for i in range(len(rule_genes)) if not neg_idx[i]])
            if rule.strip().startswith("("):
                updated_rule = f"({updated_rule})"
        scored_method = isozyme_scoring

    elif '|' not in rule: # complex case
        updated_rule = rule
        score_method = complex_scoring
    else:
        raise ValueError('Mixed AND/OR rules must be handled by _process_complex_rule.')

    # Score updated rule
    updated_genes = list(set(re.findall(r'[^&|\(\) ]+', updated_rule)))
    updated_indices = [genes.index(g) for g in updated_genes if g in genes]
    if score_method.lower() == 'min':
        r_score = np.nanmin(g_scores[updated_indices])
    elif scored_method.lower() == 'max':
        r_score = np.nanmax(g_scores[updated_indices])
    elif scored_method.lower() == 'median':
        r_score = np.nanmedian(g_scores[updated_indices])
    elif scored_method.lower() == 'average':
        r_score = np.nanmean(g_scores[updated_indices])
    else:
        raise ValueError(f"Invalid scoring method: {scored_method}")
    return updated_rule, r_score

def _process_complex_rule(rule, genes, g_scores, isozyme_scoring, complex_scoring):
    """Handle rules with both ANDs and ORs by breaking them into subsets."""
    search_patterns = [
        r'|([^&|\(\) ]+( & [^|\(\) ]+)+\)', # all ANDs
        r'|([^&|\(\) ]+( \| [^|\(\) ]+)+\)' # all ORs
    ]

    subsets = []
    c = 1
    r_orig = None
    for _ in range(100): # arbitrary high number
        if r_orig == rule:
            break
        r_orig = rule
        for pat in search_patterns:
            matches = re.findall(pat, rule)
            if matches:
                subsets.extend(matches)
                for  m in matches:
                    rule = re.sub(pat, f"#{c}#", rule, count=1)
                    c+=1

    subsets.append(rule)

    # Process subsets recursively
    for i in range(len(subsets)):
        new_rule, subset_score = _process_simple_rule(subsets[i], genes, g_scores, isozyme_scoring, complex_scoring)
        subsets[i] = new_rule
        g_scores = np.append(g_scores, subset_score)
        genes = genes + [f"#{i+1}#"]

    # Reconstruction updated rule
    updated_rule = subsets[-1]
    for i in range(c - 1, 0, -1):
        updated_rule = updated_rule.replace(f"#{i}#", subsets[i-1])

    return updated_rule

