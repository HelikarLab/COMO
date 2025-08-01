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

import copy

from cobra import Model
from cobra.flux_analysis import find_blocked_reactions, pfba


def prep_init_model(
    orig_ref_model: Model,
    task_struct: list = [],
    spont_rxn_names: list = [],
    convert_genes: bool = False,
    custom_rxns_to_ignore: list = [],
    ext_comp: str = "e",
    skip_scaling: bool = False,
):
    print("Step 1: Gene rules")
    ref_model = copy.deepcopy(orig_ref_model)

    print("Step 2: First simplification - removing blocked reactions")
    blocked_rxns = find_blocked_reactions(ref_model)
    c_model = ref_model.copy()
    c_model.remove_reactions(blocked_rxns, remove_orphans=True)

    print("Step 3: Check tasks (~10 min)")
    essential_rxns = set()
    essential_mets = set()
    task_report = []

    for task in task_struct:
        task_success, flux_dist = check_task(c_model, task)
        task_report.append({"ok": task_success, "task": task})
        if task_success:
            for rxn_id, flux in flux_dist.items():
                if abs(flux) > 1e-6:
                    essential_rxns.add(rxn_id)
            essential_mets.update(task["inputs"] + task["outputs"])

    unused_mets = [met.id for met in c_model.metabolites if all(abs(v) < 1e-9 for v in met.reactions)]
    unused_mets = list(set(unused_mets) - essential_mets)
    c_model.remove_metabolites(unused_mets)

    print("Step 4: Second simplification")
    min_model = c_model.copy()

    if not skip_scaling:
        min_model = scale_model_for_init(min_model)

    print("Step 5: Identify reactions to ignore")
    to_ignore_exch = set(r.id for r in min_model.exchanges)
    to_ignore_spont = set(spont_rxn_names)
    to_ignore_import = set()
    to_ignore_transp = set()

    for rxn in min_model.reactions:
        mets = list(rxn.metabolites)
        comps = [met.compartment for met in mets]

        if len(mets) == 2 and rxn.gene_reaction_rule == "":
            names = [met.name for met in mets]
            if names[0] == names[1] and comps[0] != comps[1]:
                to_ignore_import.add(rxn.id)

        if rxn.gene_reaction_rule == "" and all(met.compartment == ext_comp for met in mets):
            to_ignore_transp.add(rxn.id)

    to_ignore_custom = set(custom_rxns_to_ignore)

    prep_data = {
        "taskReport": task_report,
        "essentialRxns": list(essential_rxns),
        "taskStruct": task_struct,
        "refModel": c_model,
        "minModel": min_model,
        "refModelWithBM": c_model.copy(),
        "essentialMetsForTasks": list(essential_mets),
        "toIgnoreExch": list(to_ignore_exch),
        "toIgnoreImportRxns": list(to_ignore_import),
        "toIgnoreSimpleTransp": list(to_ignore_transp),
        "toIgnoreSpont": list(to_ignore_spont),
        "toIgnoreCustomRxns": list(to_ignore_custom),
    }

    return prep_data


def check_task(model: Model, task: dict):
    """Evaluates whether a metabolic task can be performed.

    task = {
        'inputs': [metabolite IDs],
        'outputs': [metabolite IDs],
        'equation': 'A[e] + B[e] => C[c] + D[c]'
    }

    Returns:
        success (bool), flux_distribution (dict)

    """
    task_model = model.copy()

    try:
        from cobra import Reaction

        # Create task reaction
        reaction = Reaction("TASK_RXN")
        reaction.build_reaction_from_string(task["equation"], task_model)
        task_model.add_reactions([reaction])

        # Set the task reaction as the objective
        task_model.objective = reaction
        solution = pfba(task_model)

        if solution.status != "optimal" or solution.fluxes[reaction.id] < 1e-6:
            return False, {}

        flux_dist = solution.fluxes.to_dict()
        return True, flux_dist

    except Exception as e:
        print(f"Error parsing task equation: {task['equation']}")
        print(f"Exception: {e}")
        return False, {}


def scale_model_for_init(model):
    for rxn in model.reactions:
        max_coeff = max(abs(val) for val in rxn.metabolites.values())
        if max_coeff > 0 and max_coeff != 1:
            scaled_coeffs = {met: coeff / max_coeff for met, coeff in rxn.metabolites.items()}
            rxn.add_metabolites(scaled_coeffs, combine=False)
    return model
