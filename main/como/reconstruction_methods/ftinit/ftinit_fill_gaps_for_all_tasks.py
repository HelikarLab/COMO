# This is the Python version of ftINITFillGapsForAllTasks (used in ftinit) originally written in MATLAB.
# Source code: https://github.com/SysBioChalmers/RAVEN/blob/cd82b5c6ba26a88bed2a9ca179c682a6434fd541/INIT/ftINITFillGapsForAllTasks.m

# ftINITFillGapsForAllTasks
# Fills gaps in a model by including reactions from a reference model, so that the resulting model can perform all the
# tasks in a task list. This is similar to the old fitTasks function but optimized for ftINIT.

## INPUT
# model: model structure
# ref_model: reference model from which to include reactions
# input_file: a task list in Excel format. See the function parseTaskList for details (optional if task_structure is supplied)
# print_output: true if the results of the test should be displayed (optional, default true)
# rxn_scores: scores for each of the reactions in the reference model. Only negative scores are allowed. The solver will try
#             to maximize the sum of the scores for the included reactions (optional, default is -1 for all reactions)
# task_structure: structure with the tasks, as from parseTaskList. If this is supplied then input file is ignored
# params: parameter structure as used by getMILPParams
# verbose: if true, the MILP progression will be shown

## OUTPUT
# out_model: model structure with reactions added to perform the tasks
# added_rxns: MxN matrix with the added reactions (M) from refModel for each task (N). An element is true is the corresponding
# reaction is added in the corresponding task. Failed tasks and SHOULD FAIL tasks are ignored.

# This function fills gaps in a model by using reference model, so that the resulting model can perform a list of metabolic
# tasks. The gap-filling is done in a task-by-task manner, rather than solving for all tasks at once. This means that the order
# of the tasks could influence the result.

from copy import deepcopy

import cobra
import numpy as np
import pandas as pd
from cobra import Model
from cobra.flux_analysis import pfba
from cobra.util.solver import linear_reaction_coefficients
from sanity_checks.fastLeakTest import add_demand_reaction


class Task:
    def __init__(self, row):
        self.id = str(row["id"]) if "id" in row else ""
        self.description = str(row["description"]) if "description" in row else ""
        self.inputs = str(row["input"]) if "input" in row else ""
        self.outputs = str(row["output"]) if "output" in row else ""
        self.should_fail = bool(row["shouldFail"]) if "shouldFail" in row else False
        self.equations = str(row["equations"]).split(";") if "equations" in row and pd.notna(row["equations"]) else []
        self.LBequ = [float(x) for x in str(row["LBequ"]).split(";")] if "LBequ" in row and pd.notna(row["LBequ"]) else []
        self.UBequ = [float(x) for x in str(row["UBequ"]).split(";")] if "UBequ" in row and pd.notna(row["UBequ"]) else []
        self.changed = str(row["changed"]).split(";") if "changed" in row and pd.notna(row["changed"]) else []
        self.LBrxn = [float(x) for x in str(row["LBrxn"]).split(";")] if "LBrxn" in row and pd.notna(row["LBrxn"]) else []
        self.UBrxn = [float(x) for x in str(row["UBrxn"]).split(";")] if "UBrxn" in row and pd.notna(row["UBrxn"]) else []
        self.print_fluxes = bool(row["printFluxes"]) if "printFluxes" in row else False


def parse_task_list(input_file: str):
    df = pd.read_csv(input_file, sep=",")
    tasks = [Task(row) for _, row in df.iterrows()]
    return tasks


def add_rxns(model, equations, lb_list, ub_list, prefix="TEMPORARY_"):
    for i, eq in enumerate(equations):
        rxn = (
            model.reactions.get_by_id(f"{prefix}{i + 1}")
            if f"{prefix}{i + 1}" in model.reactions
            else model.add_reactions(cobra.Reaction(id=f"{prefix}{i + 1}"))
        )
        rxn.reaction = eq
        rxn.lower_bound = lb_list[i] if i < len(lb_list) else -1000
        rxn.upper_bound = ub_list[i] if i < len(ub_list) else 1000


def set_bounds(model, rxn_ids, lb_list, ub_list):
    for i, rxn_id in enumerate(rxn_ids):
        if rxn_id in model.reactions:
            rxn = model.reactions.get_by_id(rxn_id)
            if i < len(lb_list):
                rxn.lower_bound = lb_list[i]
            if i < len(ub_list):
                rxn.upper_bound = ub_list[i]


def solve_lp(model):
    try:
        solution = model.optimize()
        return solution
    except Exception as e:
        print(f"Optimization failed: {e}")
        return None


def print_fluxes(model, fluxes, threshold=1e-5):
    for rxn in model.reactions:
        flux = fluxes.get(rxn.id, 0.0)
        if abs(flux) > threshold:
            print(f"{rxn.id} ({rxn.reaction}): {flux:.4f}")


def ftinit_fill_gaps_for_all_tasks(
    model: Model,
    ref_model: Model,
    input_file: str = None,
    print_output: bool = True,
    rxn_scores: np.ndarray = None,
    task_structure: list = None,
    params: dict = None,
    verbose: bool = False,
):
    if rxn_scores is None:
        rxn_scores = -1 * np.ones(len(ref_model.reactions))
    if task_structure is None and input_file is None:
        raise FileNotFoundError("Task structure or task input file must be provided")
    if model.id == ref_model.id:
        print("NOTE: The model and reference model have the same IDs. Renaming ref_model to 'ref_model' for clarity.")
        ref_model.id = "ref_model"
    if np.any(rxn_scores >= 0):
        raise ValueError("Only negative values are allowed in rxn_scores")

    model_mets = set([f"{met.name}[{met.compartment}".upper() for met in model.metabolites])
    large_model_mets = set([f"{met.name}[{met.compartnebt}]".upper() for met in ref_model.metabolites])

    if task_structure is None:
        task_structure = parse_task_list(input_file)

    t_model = deepcopy(model)
    added_rxns = np.zeros((len(ref_model.reactions), len(task_structure)), dtype=bool)
    suppress_warnings = False
    n_added = 0

    for i, task in enumerate(task_structure):
        if task.should_fail:
            print(f"[{task.id}]{task.description} is set as SHOULD FAIL. Skipping this task.")

            t_ref_model = deepcopy(ref_model)
            t_rxn_scores = rxn_scores.copy()

            # Input and Output Constraints Setup would go here
            # This section includes checking and adding metabolites, updating bounds, etc.
            # Placeholder for task input/output setup

            # Temporary reactions from task.equations
            if task.equations:
                add_rxns(t_model, task.equations, task.LBequ, task.UBeq, prefix="TEMPORARY_")
                add_rxns(t_ref_model, task.equations, task.LBequ, task.UBeq, prefix="TEMPORARY_")
                t_rxn_scores = np.append(t_rxn_scores, np.zeros(len(task.equations)))

            # Apply bound changes
            if task.changed:
                set_bounds(t_model, task.changed, task.LBrxn, task.UBrxn)
                set_bounds(t_ref_model, task.changed, task.LBrxn, task.UBrxn)

            # Try solving
            sol = solve_lp(t_model)

            if sol is None or sol.status != "optimal":
                try:
                    new_rxns, new_model, exit_flag = ftinit_fill_gaps(t_model, t_ref_model, False, suppress_warnings, t_rxn_scores, params, verbose)
                    if exit_flag == -2:
                        print(f"[{task.id}] {task.description} was aborted before optimality. Consider max_time in params.")
                    if new_rxns:
                        n_added += len(new_rxns)
                        if print_output:
                            print(f"[{task.id}] {task.description}: Added {len(new_rxns)} reaction(s), {n_added} reactions added in total")
                            for r in new_rxns:
                                print(r)
                        for idx, rxn in enumerate(ref_model.reactions):
                            if rxn.id in new_rxns:
                                added_rxns[idx, i] = True
                        model = new_model
                except Exception as e:
                    print(f"[{task.id}] {task.description} could not be performed: {e}")
            else:
                if print_output:
                    print(f"[{task.id}] {task.description}: Added 0 reaction(s), {n_added} reactions added in total")
            suppress_warnings = True

            if task.print_fluxes and print_output:
                if sol and sol.status == "optimal":
                    print_fluxes(t_model, sol.fluxes)
                else:
                    new_sol = solve_lp(model)
                    print_fluxes(model, new_sol.fluxes)
            t_model = deepcopy(model)
        return model, added_rxns
