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

import os
import copy
from sys import prefix

import numpy as np
from typing import List, Tuple, Optional
from cobra import Model, Reaction, Metabolite
from nbclient.exceptions import timeout_err_msg


class Task:
    """Container matching the MATLAB taskStructure field used here.
    Provide objects with these attributes, or adapt as needed"""
    def __init__(self,
                 id: str,
                 description: str = "",
                 shouldFail: bool = False,
                 printFluxes: bool = False,
                 inputs: Optional[List[str]] = None,
                 UBin: Optional[List[float]] = None,
                 LBin: Optional[List[float]] = None,
                 outputs: Optional[List[str]] = None,
                 UBout: Optional[List[float]] = None,
                 LBout: Optional[List[float]] = None,
                 equations: Optional[List[str]] = None,
                 LBequ: Optional[List[float]] = None,
                 UBequ: Optional[List[float]] = None,
                 changed: Optional[List[str]] = None,
                 LBrxn: Optional[List[float]] = None,
                 UBrxn: Optional[List[float]] = None):
        self.id = id
        self.description = description
        self.shouldFail = shouldFail
        self.printFluxes = printFluxes
        self.inputs = inputs or []
        self.UBin = np.array(UBin or [], dtype=float)
        self.LBin = np.array(LBin or [], dtype=float)
        self.outputs = outputs or []
        self.UBout = np.array(UBout or [], dtype=float)
        self.LBout = np.array(LBout or [], dtype=float)
        self.equations = equations or []
        self.LBequ = np.array(LBequ or [], dtype=float)
        self.UBequ = np.array(UBequ or [], dtype=float)
        self.changed = changed or []
        self.LBrxn = np.array(LBrxn or [], dtype=float)
        self.UBrxn = np.array(UBrxn or [], dtype=float)

def disp_em(message: str, is_error: bool = True):
    prefix = "ERROR" if is_error else "WARNING: "
    print(prefix + str(message))

def parse_task_list(input_file: str) -> List[Task]:
    """ User must implement according to their task file format. Should return List[Task]"""
    raise NotImplementedError("parse_task_list needs to be implemented for your task format")

def add_mets(model: Model, met: dict) -> Model:
    """Add metabolites to COBRApy model. Expects keys 'met_names' and 'compartments'. The metabolite id is synthesized from name + compartment"""
    met_names = met.get('met_names', [])
    comps = met.get('compartments', [])
    mets_to_add: List[Metabolite] = []
    for name, comp in zip(met_names, comps):
        # ID: sanitize name for an id; you can customize this mapping
        met_id  = f"{name.strip().replace(' ','_').replace('[','').replace(']','')}_{comp}"
        if hasattr(model, 'metabolites'):
            if any(m.id == met_id for m in model.metabolites):
                # already exists; skip creating duplicate
                continue
        m = Metabolite(id=met_id, name=name, compartment=comp)
        mets_to_add.append(m)
    if mets_to_add:
        model.add_metabolites(mets_to_add)
    return model

def add_rxns(model: Model, rxn: dict) -> Model:
    """Add reactions to CPBRApy model from equations + bounds. Excepts keys 'equations','lb','ub','rxns'."""
    eqs = rxn.get('equations',[])
    lbs = rxn.get('lb',[])
    ubs = rxn.get('ub',[])
    rids = rxn.get('rxns',[])
    for rid, eq, lb, ub in zip(rids, eqs, lbs, ubs):
        if any(r.id == rid for r in model.reactions):
            # If exsits, update equation/bounds
            r = model.reactions.get_by_id(rid)
        else:
            r = Reaction(rid)
            model.add_reactions([r])
        r.lower_bound = float(lb)
        r.upper_bound = float(ub)
        # COBRApy expects a string like "A_c + 2 B_c -> C_c"
        r.build_reaction_from_string(eq, verbose=False)
    return Model

def set_param(model: Model, bound_type: str, rxn_ids: List[str], values: np.ndarray) -> Model:
    for rid, val in zip(rxn_ids, values):
        r = model.reactions.get_by_id(rid)
        if bound_type.lower() == "lb":
            r.lower_bound = float(val)
        elif bound_type.lower() == "ub":
            r.upper_bound = float(val)
        else:
            raise ValueError(f"Unknown bound type: {bound_type}")
    return model

def solve_lp(model:Model, objective_sense: int = 0):
    """Solve LP. objective_sense: 0 = default (as set on model), 1 = maximize. Returns COBRApy Solution, or a simple namespace with .x=None on failure."""
    try:
        if objective_sense == 1:
            # COBRApy maximize by default; if model.objective has coefficients, this is fine
            sol = model.optimize(objective_sense='maximize')
        else:
            sol = model.optimize()
        return sol
    except Exception:
        class _Empty:
            x = None
        return _Empty

def print_fluxes(model: Model, flux_vector, *_args, **_kwargs):
    if hasattr(flux_vector, 'items'):
        items = flux_vector.items()
    else:
        items = zip([r.id for r in model.reactions], flux_vector)
    for rid, val in items:
        r = model.reactions.get_by_id(rid)
        print(f"{r.id} ({r.reaction}): {float(val):.6g}")

##---------------------------------------------GAP-FILLING CALL BACK PLACEHOLDER---------------------------------------------##

def ftinit_fill_gaps(t_model: Model, base_model: Model, t_ref_model: Model, flag: bool, suppress_warnings: bool, t_rxn_scores: np.ndarray, params: dict, verbose: bool) -> Tuple[List[str], Model, int]:
    """
    User must suplly their gap-filling implementation.
    Should return (new_rxns, new_model, exit_flag)
    new_rxns: list of reaction IDs added
    new_model: the gap  filling model
    exit_flag: -2 if aborted before optimality (to mirror MATLAB), 1 if success
    """
    raise NotImplementedError("Provide your ftINITFillGaps implementation")

##---------------------------------------------UTILITY FUNCTIONS---------------------------------------------##
