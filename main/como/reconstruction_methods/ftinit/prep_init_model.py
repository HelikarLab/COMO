# This is the Python version of prepINITModel (used in ftinit) originally written in MATLAB.
# Source code: https://github.com/SysBioChalmers/RAVEN/blob/cd82b5c6ba26a88bed2a9ca179c682a6434fd541/INIT/prepINITModel.m

# prep_init_model: The purpose of this function is to run time-consuming calculation steps that are not dependent on the RNA-Seq data

## INPUT
# orig_ref_model: The model o use. Expected to be something such as Human_GEM, Mouse_GEM, etc.
# task_structure: The essential tasks. Can be loaded with for example
#                 task_structure = pars_task_list(../data/metabolicTasks_Essential.txt)
# spont_rxn_names: The spontaneous rxns. (optional, default {})
# convert_genes: If true the gees are converted to gene names (from ENSEMBL) (optional, default false)
#                custom_rxns_to_ignore: These reactions can be ignored in the ignore mask (specifying b7=1) (optional, default = {})
#                ext_comp: Name of the external compartment, typically 's' or 'e'. This is used for identifying exch and import rxns
#                (optional, default = 'e')
# skip_scaling: If true the scaling step is not run on the minimal model. The scaling is there to remove large differences between the stoichiometric
#                coefficients within a reaction, since such differences creates numerical issues in the ftINIT algorithm due to limitations
#                 in solver resolution. However, the current scaling step is also risky and may lead to that some reactions cannot carry flux, since
#                 it changes the stoichiometry of the reaction with large differences. If you experience problems where the solution is infeasible,
#                 it may be worth trying to turn off the scaling. Note that it is only the min_model that is scaled, the scaling will not be present
#                 in the final model. Default: (optional, default = false)

## OUTPUT
# prep_data: The resulting prep_data structure which is used as input to ftINIT.

import cobra
import numpy as np
import sympy
from cobra import Model, Reaction
from typing import List, Tuple

from reconstruction_methods.ftinit.reverse_rxns import reverse_rxns

###-----------------------------------------------Helper Functions-----------------------------------------------###
def standardize_gr_rules(model: Model, verbose: bool = False) -> Tuple[List[str], np.ndarray]:
    """
    Standardize gene rules in the model to consistent format

    Args:
        model: CobraPy model
        verbose: Print progress messages

    Returns:
        Tuple of (standardized rules, reaction-gene matrix)
    """
    if verbose:
        print("Standardizing gene rules...")

    # Initialize empty rules if none exist
    if not hasattr(model, 'grRules') or not model.grRules:
        model.grRules = ['' for _ in model.reactions]

    # Create reaction-gene matrix
    rxn_gene_mat = np.zeros((len(model.reactions), len(model.genes)), dtype=bool)

    for i, rxn in enumerate(model.reactions):
        if hasattr(rxn, 'gene_reaction_rule'):
            # Standardize AND/OR formatting
            rule = rxn.gene_reaction_rule
            rule = rule.replace(' and ', ' & ').replace(' or ', ' | ')
            model.grRules[i] = rule

            # Update reaction-gene matrix
            for gene in rxn.genes:
                gene_idx = model.genes.index(gene.id)
                rxn_gene_mat[i, gene_idx] = True

    return model.grRules, rxn_gene_mat

def simplify_model(model: Model, remove_dead_ends: bool = True, remove_rev: bool = False,
                   remove_orphans: bool = True, remove_blocked: bool = True,
                   remove_irrelevant: bool = True, verbose: bool = False) -> Model:
    """
    Simplify model by removing various reaction types

    Args:
        model: CobraPy model
        remove_dead_ends: Remove dead-end metabolites
        remove_rev: Remove reversible reactions
        remove_orphans: Remove orphan metabolites
        remove_blocked: Remove blocked reactions
        remove_irrelevant: Remove irrelevant reactions
        verbose: Print progress messages

    Returns:
        Simplified model
    """
    simplified = model.copy()

    if remove_dead_ends:
        if verbose: print("Removing dead-end metabolites...")
        dead_ends = [met for met in simplified.metabolites
                     if len(met.reactions) == 0]
        simplified.remove_metabolites(dead_ends)

    if remove_rev:
        if verbose: print("Removing reversible reactions...")
        rev_rxns = [rxn for rxn in simplified.reactions if rxn.reversibility]
        simplified.remove_reactions(rev_rxns)

    if remove_blocked:
        if verbose: print("Removing blocked reactions...")
        blocked = cobra.flux_analysis.variability.find_blocked_reactions(simplified)
        simplified.remove_reactions(blocked)

    return simplified
def remove_reactions(model: Model, rxn_ids: List[str], remove_metabolites: bool = False,
                     update_gprs: bool = True) -> Model:
    """
    Remove reactions from model

    Args:
        model: CobraPy model
        rxn_ids: List of reaction IDs to remove
        remove_metabolites: Remove unused metabolites
        update_gprs: Update gene rules after removal

    Returns:
        Model with reactions removed
    """
    new_model = model.copy()
    rxns_to_remove = [rxn for rxn in new_model.reactions if rxn.id in rxn_ids]
    new_model.remove_reactions(rxns_to_remove,
                               remove_metabolites=remove_metabolites)

    if update_gprs:
        standardize_gr_rules(new_model)

    return new_model

def check_tasks(model: Model, tasks: List[dict] = None, verbose: bool = False,
                check_essential: bool = True, check_all: bool = False) -> Tuple:
    """
    Check if model can perform metabolic tasks

    Args:
        model: CobraPy model
        tasks: List of task dictionaries
        verbose: Print progress
        check_essential: Check essential reactions
        check_all: Check all reactions

    Returns:
        Tuple of (task report, essential reaction matrix, essential metabolites, essential fluxes)
    """
    if tasks is None:
        return [], np.zeros((len(model.reactions), 0)), [], np.zeros((len(model.reactions), 0))

    task_report = []
    essential_rxn_mat = np.zeros((len(model.reactions), len(tasks)), dtype=bool)
    essential_fluxes = np.zeros((len(model.reactions), len(tasks)))

    for i, task in enumerate(tasks):
        if verbose: print(f"Checking task {i+1}/{len(tasks)}: {task['id']}")

        # Convert task to cobra constraints
        with model as temp_model:
            # Apply task constraints
            for rxn_id, bounds in task['constraints'].items():
                if rxn_id in temp_model.reactions:
                    temp_model.reactions.get_by_id(rxn_id).bounds = bounds

            # Solve and check feasibility
            solution = temp_model.optimize()
            task_report.append({
                'id': task['id'],
                'ok': solution.status == 'optimal',
                'flux': solution.objective_value
            })

            if check_essential and solution.status == 'optimal':
                # Find essential reactions for this task
                for j, rxn in enumerate(temp_model.reactions):
                    with temp_model:
                        rxn.knock_out()
                        sol = temp_model.optimize()
                        if sol.status != 'optimal' or sol.objective_value < 0.01 * solution.objective_value:
                            essential_rxn_mat[j, i] = True
                            essential_fluxes[j, i] = solution.fluxes[rxn.id]

    return task_report, essential_rxn_mat, [], essential_fluxes

def merge_linear(model: Model, rxns_to_merge: List[str] = None) -> Tuple[Model, List[str], List[int]]:
    """
    Merge linearly dependent reactions

    Args:
        model: CobraPy model
        rxns_to_merge: Specific reactions to merge (None for all)

    Returns:
        Tuple of (merged model, original reaction IDs, group IDs)
    """
    merged = model.copy()
    S = cobra.util.array.create_stoichiometric_matrix(merged)

    if rxns_to_merge is None:
        rxns_to_merge = [rxn.id for rxn in merged.reactions]

    # Find linearly dependent reactions
    _, pivots = sympy.Matrix(S.T).rref()
    independent = [i for i in range(S.shape[1]) if i in pivots]
    dependent = [i for i in range(S.shape[1]) if i not in pivots]

    # Create mapping of merged reactions
    orig_rxn_ids = []
    group_ids = []
    next_group = 1

    for i in independent:
        orig_rxn_ids.append(merged.reactions[i].id)
        group_ids.append(next_group)
        next_group += 1

    for i in dependent:
        # Find which independent reaction this depends on
        for j in independent:
            if np.allclose(S[:,i], S[:,j]):
                orig_rxn_ids.append(merged.reactions[i].id)
                group_ids.append(group_ids[j])
                break

    # Create new model with merged reactions
    new_rxns = []
    for group in np.unique(group_ids):
        if group > 0:  # 0 means not merged
            members = [i for i, g in enumerate(group_ids) if g == group]
            if len(members) > 1:
                # Create new merged reaction
                base_rxn = merged.reactions[members[0]]
                new_rxn = Reaction(f"merged_{group}")
                new_rxn.add_metabolites(base_rxn.metabolites)
                new_rxn.gene_reaction_rule = " or ".join([merged.reactions[i].gene_reaction_rule for i in members])
                new_rxns.append(new_rxn)

    # Build new model
    new_model = Model()
    new_model.add_metabolites([met.copy() for met in merged.metabolites])
    new_model.add_reactions([rxn.copy() for rxn in merged.reactions if rxn.id not in rxns_to_merge])
    new_model.add_reactions(new_rxns)

    return new_model, orig_rxn_ids, group_ids

def rescale_model_for_init(model: Model) -> Model:
    """
    Rescale model stoichiometry to improve numerical stability

    Args:
        model: CobraPy model

    Returns:
        Rescaled model
    """
    scaled = model.copy()
    for rxn in scaled.reactions:
        # Find greatest common divisor of coefficients
        coeffs = [abs(v) for v in rxn.metabolites.values()]
        if len(coeffs) > 0:
            gcd_val = np.gcd.reduce([int(round(c*100)) for c in coeffs])/100
            if gcd_val > 0:
                # Scale reaction
                for met in rxn.metabolites:
                    rxn.add_metabolites({met: -rxn.metabolites[met]/gcd_val}, combine=False)
    return scaled

###-----------------------------------------------Main Function-----------------------------------------------###

def prep_init_model(
    orig_ref_model: Model, task_structure, spont_rxn_names=None, convert_genes=False, custom_rxns_to_ignore=None, ext_comp="e", skip_scaling=False
):
    if spont_rxn_names is None:
        spont_rxn_names = []
    if custom_rxns_to_ignore is None:
        custom_rxns_to_ignore = []

    print("Step 1: Gene rules")
    orig_ref_model.gr_rules, orig_ref_model.rxn_gene_mat = standardize_gr_rules(orig_ref_model, True)

    if convert_genes:
        orig_ref_model.gr_rules, orig_ref_model.genes, orig_ref_model.rxn_gene_mat = translate_gr_rules(orig_ref_model.gr_rules, "None")

    print("Step 2: First Simplification")
    _, deleted_dead_end_rxns = simplify_model(orig_ref_model, True, False, True, True, True)

    # Get reduced model
    c_model = remove_reactions(orig_ref_model, deleted_dead_end_rxns, False, True)

    print("Step 3: Check tasks (~10min)")
    if task_structure:
        b_model = close_model(c_model)
        task_report, essential_rxn_mat, _, essential_fluxes = check_tasks(b_model, None, True, False, True, task_structure)

        # Extract essential reactions
        sel = np.sum(essential_rxn_mat, axis=1) > 0
        sel_ind = np.where(sel)[0]
        essential_rxns = b_model.rxns[sel_ind]

        # Find metabolites present in task_structure
        task_mets = np.union1d(np.concatenate(task_structure.inputs), np.concatenate(task_structure.outputs))
        task_mets = np.union1d(task_mets, parse_rxn_equ(np.concatenate(task_structure.equations)))
        model_mets = [f"{name}[{comp}]" for name, comp in zip(c_model.met_Names, c_mocel.comps[c_model.met_comps])]
        in_model, met_ind = np.isin(task_mets, model_mets, return_indices=True)
        essential_mets_for_tasks = c_model.mets[met_ind[in_model]]

        # Remove tasks that cannot be performed
        task_structure = [task for task, ok in zip(task_structure, task_report.ok) if ok]
    else:
        essential_rxns = []
        essential_mets_for_tasks = []
        task_report = []

    # Remove metabolites separately to avoid removing those needed for tasks
    unused_mets = c_model.mets[np.all(c_model.S == 0, axis=1)]
    c_model = remove_mets(c_model, np.setdiff1d(unused_mets, essential_mets_for_tasks))

    if task_structure:
        # Determine direction for essential reactions
        essential_rev_dir = np.zeros(len(essential_rxns), dtype=bool)
        pp = np.zeros(len(essential_rxns))
        nn = np.zeros(len(essential_rxns))

        for i in range(len(sel_ind)):
            pos = np.sum(essential_fluxes[sel_ind[i], essential_rxn_mat[sel_ind[i], :]] > 0)
            neg = np.sum(essential_fluxes[sel_ind[i], essential_rxn_mat[sel_ind[i], :] < 0])
            essential_rev_dir[i] = pos, neg
            pp[i] = pos
            nn[i] = neg

        # Create a minimal model for the MILP
        min_model1 = c_model.copy()
        min_model1 = reverse_rxns(min_model1, min_model1.rxns[sel_ind[essential_rev_dir]])
        min_model1.rev[sel_ind] = 0
        min_model1.lb[sel_ind] = 0

    else:
        min_model1 = c_model

    print('Step 4: Second simplification (~1 hour)')
    min_model2 = simplify_model(min_model1, False, False, False, False, False, False, True)

    print('Step 5: Final Work')
    # Merge all linear dependednt reactions
    min_model3, orig_rxn_ids, group_ids = merge_linear(min_model2, [])

    if task_structure:
        # Find all potential reactions
        rxn_ind_orig = np.isin(orig_rxn_ids, essential_rxns)
        ids = np.unique(group_ids[rxn_ind_orig])
        ids = ids[ids > 0 ]
        rxn_candidates = np.unique(np.concatenate([essential_rxns, orig_rxn_ids[np.isin(group_ids, ids)]]))
        new_essential_rxns = min_model3.rxns[np.isin(min_model3.rxns, rxn_candidates)]
    else:
        new_essential_rxns = []

    # Identify reactions to ignore in the minModel2
    print('Identifying reactions to ignore...')
    exch_rxn_ind = get_exchange_rxns(min_model2)
    to_ignore_exch = np.zeros(len(min_model2.rxns), dtype=bool)
    to_ignore_exch[exch_rxn_ind] = True
    to_ignore_import_rxns = np.zeros(len(min_model2.rxns), dtype=bool)
    to_ignore_simple_transp = np.zeros(len(min_model2.rxns), dtype=bool)
    to_ignore_adv_transp = np.zeros(len(min_model2.rxns), dtype=bool)
    to_ignore_s = np.zeros(len(min_model2.rxns), dtype=bool)

    # Identify simple transport reactions with no GPRs
    num_mets = np.sum(min_model2.S != 0, axis=1)
    scomp = np.where(min_model2.comps == ext_comp)[0]][0]

    for i in range(len(min_model2.rxns)):
        if min_model2.gr_rules[i] == ' ' not num_mets[i] == 2:
            met_comps = min_model2.met_comps[min_model2.S[:,i] != 0]
            met_names = min_model2.met_names[min_model2.S[:, i] != 0]
            if mets_comps[0] != mets_comps[1] and (mets_comps[0] == scomp or mets_comps[1] == scomp):
                if met_names[0] == met_names[1]:
                    to_ignore_import_rxns[i] = True
            elif mets_comps[0] != mets_comps[1]:
                if met_names[0] == met_names[1]:
                    to_ignore_simple_transp[i] = True
        else:
            mets_ind = np.where(min_model2.S[:,i] != 0)[0]
            if len(mets_ind) % 2 == 0 and len(mets_ind) > 2 and min_model2.gr_rules[i] == ' ':
                sv_vals = min_model2.S[mets_ind, i].A.flatten()
                met_names = min_model2.met_names[mets_ind]
                comps = min_model2.met_comps[mets_ind]
                success = True
                while len(met_names) > 0:
                    if len(met_names) < 2:
                        raise ValueError("We should never arrive here")
                    met_match = np.where(met_names[1:] == met_names[0])[0] + 1
                    if len(met_match) != 1:
                        success = False
                        break
                    if (sv_vals[0] + sv_vals[met_match]) != 0:
                        success = False
                        break
                    if comps[0] == comps[met_match]:
                        success = False
                        break
                    sv_vals = np.delete(sv_vals, [0, met_match])
                    mat_names = np.delete(met_names, [0, met_match])
                    comps = np.delete(comps, [0, met_match])
                to_ignore_adv_transp[i] = success
    # Spontaneous reactions
    to_ignore_spont = np.isin(min_model2.rxns, spont_rxn_names)

    # Reactions without GPRs in the external compartment
    s_comp_ind = np.where(min_model2.comps == ext_comp)[0]
    for i in range(len(min_model2.rxns)):
        mets_ind = np.where(min_model2.S[:,i] != 0)[0]
        if len(mets_ind) > 1 and min_model2.gr_rules[i] == ' ':
            if np.all(min_model2.met_comps[mets_ind] == s_comp_ind):
                to_ignore_s[i] = True

    to_ignore_all_without_gprs = np.array([rule == '' for rule in min_model2.gr_rules])

    # Scaling the model
    if skip_scaling:
        scaled_min_model = min_model3
    else:
        scaled_min_model = rescale_model_for_init(min_model3)
        scaled_min_model.ub[scaled_min_model.ub > 0] = 1000
        scaled_min_model.lb[scaled_min_model.lb < 0] = -1000

    # Create data structure with preprocessing results
    prep_data = {
        'task_report': task_report,
        'essential_rxns': new_essential_rxns,
        'task_structure': task_structure,
        'ref_model': c_model,
        'min_model': scaled_min_model,
        'ref_model_with_bm': closed_model(c_model),
        'group_ids': group_ids,
        'essential_mets_for_tasks': essential_mets_for_tasks,
        'to_ignore_exch': to_ignore_exch,
        'to_ignore_import_rxns': to_ignore_import_rxns,
        'to_ignore_simple_transp': to_ignore_simple_transp,
        'to_ignore_adv_transp': to_ignore_adv_transp,
        'to_ignore_spont': to_ignore_spont,
        'to_ignore_s': to_ignore_s,
        'to_ignore_custom_rxns': np.isin(min_model2.rxns, custom_rxns_to_ignore),
        'to_ignore_all_without_gprs': to_ignore_all_without_gprs
    }
    return prep_data



