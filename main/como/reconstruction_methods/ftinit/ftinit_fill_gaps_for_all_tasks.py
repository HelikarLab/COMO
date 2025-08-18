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
from sys import prefix, excepthook

import numpy as np
from typing import List, Tuple, Optional
from cobra import Model, Reaction, Metabolite
from nbclient.exceptions import timeout_err_msg
from numba.np.new_arraymath import jit_np_isin
from umap.distances import named_distances_with_gradients


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
def _met_key(met: Metabolite) -> str:
    # Match MATLAB: NAME[comp]; using compartment key (e.g., 'c')
    return f"{met.name or met.id.upper()}[{met.compartment.upper()}"

def _upper_met_token(name_with_brackets: str) -> str:
    return name_with_brackets.upper()

def _extract_allmetsin_comp(token_upper: str) -> Optional[str]:
    # Extect from "ALLMETSIN[XXX]"
    if not token_upper.startwith('ALLMETSIN'):
        return None
    if '[' in token_upper and token_upper.endswith(']'):
        return token_upper[token_upper.index('[')+1:-1]
    return None

##---------------------------------------------ftinit_fill_gaps_for_all_tasks---------------------------------------------##
def ftinit_fill_gaps_for_all_tasks(model: Model, ref_model: Model, input_file: str, print_output: bool, rxn_scores: Optional[np.ndarray], task_structure: Optional[List[Task]], params: dict, verbose: bool) -> Tuple[Model, np.ndarray]:
    """Faithful translation of MATLAB ftINITFillGapsForAllTasks into Python(COBRApy)"""

    # Initialize rxn_scores
    if rxn_scores is None or len(rxn_scores) == 0:
        rxn_scores = np.ones(len(ref_model.reactions)) * -1
    if (not task_structure) and (not os.path.isfile(input_file)):
        raise FileNotFoundError(f"Task file {input_file} cannot be found")
    if getattr(model,'id',None) == getattr(ref_model, 'id',None):
        print('NOTE: The model and reference model have the same IDs. The ID for the reference modelwas set to "refModel" to track reaction origin.')
        setattr(ref_model, 'id','refModel')
    if np.any(rxn_scores >= 0):
        disp_em('Omly negative values are allowed in rxn_scores')

    # Prepare b matrix (cols: 0 = IN, 1 = OUT)
    # Attach as attribute; COBRApy won't use it internally but we keep it to parallel MATLAB
    model.b = np.zeros((len(model.metabolites), 2))

    model_mets = [_met_key(m) for m in model.metabolites]
    large_model_mets = [_met_key(m) for m in ref_model.metabolites]

    if not hasattr(model, 'unconstrained'):
        disp_em('Exchange metabolites should normally not be removed from the model when using check_tasks. Inputs and outputs are defined in the task '
                'file instead. Use a loader that preserves exchange reactions.', is_error=False)

    if not task_structure:
        task_structure = parse_task_list(input_file)

    t_model = copy.deepcopy(model)
    added_rxns = np.zeros((len(ref_model.reactions), len(task_structure)), dtype=bool)
    supress_warnings = False
    n_added = 0

    # Build a stable list of ref reaction IDs for indexing added_rxns
    ref_rxn_ids = [r.id for r in ref_model.reactions]

    for i, task in enumerate(task_structure):
        if not getattr(task, 'shouldFail', False):
            t_ref_model = copy.deepcopy(ref_model) # we also adjust constraints on this
            t_rxn_scores = copy.deepcopy(rxn_scores) # extended when reactions are added per task

            ##  Set the inputs
            if getattr(task, 'inputs', []):
                inputs_upper = [_upper_met_token(x) for x in task.inputs]
                I = [x in model_mets for x in inputs_upper]
                J = [model_mets.index(x) if x in model_mets else -1 for x in inputs_upper]
                I2 = [x in large_model_mets for x in inputs_upper]
                J2 = [large_model_mets.index(x) if x in large_model_mets else -1 for x in inputs_upper]
                K = [x == 'ALLMETS' for x in inputs_upper]
                L = [_extract_allmetsin_comp(x) is not None for x in inputs_upper]

                good_mets = [ a or b or c for a,b,c in zip(I, K, L)]
                if not all(good_mets):
                    missing_idx = [idx for idx, ok in enumerate(good_mets) if not ok]
                    missing_upper = [inputs_upper[idx] for idx in missing_idx]
                    found = [x in large_model_mets for x in missing_upper]
                    if not all(found):
                        disp_em(f"Could not find all inputs in \"[{task.id}] {task.description}\" in either model")
                    else:
                        met_match = [large_model_mets.index(x) for x in missing_upper]
                        met = {
                            'met_names': [ref_model.metabolites[k].name for k in met_match],
                            'compartments': [ref_model.metabolites[k].compartment for k in met_match]
                        }
                        model = add_mets(model,met)
                        t_model = add_mets(t_model, met)
                        model_mets.extend(missing_upper)
                        # Update J mapping for those newly added
                        for loc, midx in zip(missing_idx, range(len(model_mets) - len(missing_upper), len(model_mets))):
                            I[loc] = True
                            J[loc] = midx
                # check duplocates among indices J(I)
                selected = [J[idx] for idx, ok in enumerate(I) if ok and J[idx] >= 0]
                if len(set(selected)) != len(selected):
                    disp_em(f"The constraints on some input(s) in \"[{task.id}] {task.description}\" are defined more than one time")

                # ALLMETS: apply to all mets
                if any(K):
                    if not K[0]:
                        disp_em(f"ALLMETS is used as on inout in \"[{task.id}] {task.description}\" but it is not the first metabolite in the list. "
                                f"Constraints defined for the metabolites before it will be over-written", is_error=False)
                    ub_val = float(task.UBin[np.where(K)[0][0]]) * -1.0
                    t_model.b[:,0] = ub_val
                    t_ref_model.b = np.zeros((len(t_ref_model.metabolites),2)) if not hasattr(t_ref_model, 'b') else t_ref_model.b
                    t_ref_model.b[:,0] = ub_val

                # ALLMETSIN[comp]
                if any(L):
                    for lidx, flag in enumerate(L):
                        if not flag:
                            continue
                        compartment = _extract_allmetsin_comp(inputs_upper[lidx]) # already upper
                        # t_model
                        mask_t = [m.compartment.upper() == compartment for m in t_model.metabolites]
                        if any(mask_t):
                            t_model.b[np.where(mask_t)[0], 0] = float(task.UBin[lidx]) * -1.0
                        else:
                            disp_em(f"The compartment defined for ALLMETSIN in \"[{task.id}] {task.description}\" does not exits")

                        # t_ref_model
                        mask_r = [m.compartment.upper() == compartment for m in t_ref_model.metabolites]
                        t_ref_model.b = np.zeros((len(t_ref_model.metabolites), 2)) if not hasattr(t_ref_model, 'b') else t_ref_model.b
                        if any(mask_r):
                            t_ref_model.b[np.where(mask_r)[0], 0] = float(task.UBin[lidx]) * -1.0
                        else:
                            disp_em(f"The compartment defined for ALLMETSIN in \"[{task.id}] {task.description}\" does not exist")

                # Normal constraints
                for idx, ok in enumerate(I):
                    if ok and J[idx] >= 0:
                        t_model.b[J[idx], 0] = float(task.UBin[idx]) * -1.0
                        t_model.b[J[idx], 1] = float(task.LBin[idx]) * -1.0
                for idx, ok in enumerate(I2):
                    if ok and J2[idx] >= 0:
                        t_ref_model.b = np.zeros((len(t_ref_model.metabolites),2)) if not hasattr(t_ref_model, 'b') else t_ref_model.b
                        t_ref_model.b[J2[idx], 0] = float(task.UBin[idx]) * -1.0
                        t_ref_model.b[J2[idx], 1] = float(task.LBin[idx]) * -1.0

            ## set the outputs
            if getattr(task,'outputs', []):
                outputs_upper = [_upper_met_token(x) for x in task.outputs]
                I = [x in model_mets for x in outputs_upper]
                J = [model_mets.index(x) if x in model_mets else -1 for x in outputs_upper]
                I2 = [x in large_model_mets for x in outputs_upper]
                J2 = [large_model_mets.index(x) if x in large_model_mets else -1 for x in outputs_upper]
                K = [x == 'ALLMETS' for x in outputs_upper]
                L = [_extract_allmetsin_comp(x) is not None for x in outputs_upper]

                good_mets = [a or b or c for a,b,c in zip(I,K,L)]
                if not all(good_mets):
                    missing_idx = [idx for idx, ok in enumerate(good_mets) if not ok]
                    missing_upper = [outputs_upper[idx] for idx in missing_idx]
                    found = [x in large_model_mets for x in missing_upper]
                    if not all(found):
                        disp_em(f"Could not find all outputs in \"[{task.id}] {task.description}\" in either model")
                    else:
                        met_match = [large_model_mets.index(x) for x in missing_upper]
                        met = {
                            'met_names': [ref_model.metabolites[k].name for k in met_match],
                            'compartments': [ref_model.metabolites[k].compartment for k in met_match]
                        }
                        model = add_mets(model,met)
                        t_model = add_mets(t_model,met)
                        model_mets.extend(missing_upper)
                        for loc, midx in zip(missing_idx, range(len(model_mets) - len(missing_upper), len(model_mets))):
                            I[loc] = True
                            J[loc] = midx
                selected = [J[idx] for idx, ok in enumerate(I) if ok and J[idx] >= 0]
                if len(set(selected)) != len(selected):
                    disp_em(f"The constraints on some output(s) in \"[{task.id}] {task.description}\" are defined more than one time")

                # ALLMETS: apply to all mets (column 1 for OUT upper bound)
                if any(K):
                    if not K[0]:
                        disp_em(f"ALLMETS is used as an output in \"[{task.id}] {task.description}\" but it is not the first metabolite in the list. "
                                f"Constraints defined for the metabolites before it will be over-written", is_error=False)
                    ub_val = float(task.UBout[np.where(K)[0][0]])
                    t_model.b[:, 1] = ub_val
                    t_ref_model.b = np.zeros((len(t_ref_model.metabolites), 2)) if not hasattr(t_ref_model, 'b') else t_ref_model.b
                    t_ref_model.b[:,1] = ub_val

                # ALLMETSIN[comp]
                if any(L):
                    for lidx, flag in enumerate(L):
                        if not flag:
                            continue
                        compartment = _extract_allmetsin_comp(outputs_upper[lidx])
                        # t_model
                        mask_t = [m.compartment.upper() == compartment for m in t_model.metabolites]
                        if any(mask_t):
                            t_model.b[np.where(mask_t)[0], 1] = float(task.UBout[lidx])
                        else:
                            disp_em(f"The compartment defined for ALLMETSIN in \"[{task.id}] {task.description}\" does not exist")
                        # t_ref_model
                        mask_r = [m.compartment.upper() == compartment for m in t_ref_model.metabolites]
                        t_ref_model.b = np.zeros((len(t_ref_model.metabolites),2)) if not hasattr(t_ref_model, 'b') else t_ref_model.b
                        if any(mask_r):
                            t_ref_model.b[np.where(mask_r)[0], 1] = float(task.UBout[lidx])
                        else:
                            disp_em(f"The compartment defined for ALLMETSIN in \"[{task.id}] {task.description}\" does not exist")
                # Normal constraints with IN/OUT consistency check
                if any(I) and any(J[idx] >= 0 for idx in range(len(J))):
                    J_arr = np.array([J[idx] if I[idx] else -1 for idx in range(len(J))])
                    I_idx = np.where(I)[0]
                    # nonzero_LBin: t_model.b(:,2) < 0 in MATLAB (co; 1 in 0-based) -> here b[:,1] is OUT; they checked b(:,2) < 0 (that's OUT LB?)
                    # In MATLAB they: nonzero_LBin = tModel.b(J,2) < 0; but earlier they set tModel.b(J,2) = LBin*-1 for inputs => negative means required IN uptake.
                    # Here we mimic exactly:
                    nonzero_LBin = t_model.b[J_arr[I_idx], 1 ] < 0
                    nonzero_LBout = np.array([task.LBout[idx] for idx in I_idx]) > 0
                    if np.any(nonzero_LBin & nonzero_LBout):
                        disp_em(f"The IN LB and OUT LB in \"[{task.id}] {task.description}\" cannot be nonzero for the same metabolite")
                    # Apply
                    for loc_k, j_idx in zip(I_idx, J_arr[I_idx]):
                        if j_idx<0:
                            continue
                        if task.LBout[loc_k] > 0:
                            t_model.b[j_idx, 0] = float(task.LBout[loc_k]) # set IN lower bound (positive means required output? mirrors MATLAB
                        t_model.b[j_idx, 1] = float(task.UBout[loc_k])

                # and for t_ref_model
                if any(I2) and any(J2[idx] >= 0 for idx in range(len(J2))):
                    J2_arr = np.array([J2[idx] if I2[idx] else -1 for idx in range(len(J2))])
                    I2_idx = np.where(I2)[0]
                    t_ref_model.b = np.zeros((len(t_ref_model.metabolites),2)) if not hasattr(t_ref_model, 'b') else t_ref_model.b
                    nonzero_LBin = t_ref_model.b[J2_arr[I2_idx], 1] < 0
                    nonzero_LBout = np.array([task.LBout[idx] for idx in I2_idx]) > 0
                    if np.any(nonzero_LBin & nonzero_LBout):
                        disp_em(f"The IN LB and OUT LB in \"[{task.id}] {task.description}\" cannot be nonzero for the same metabolite")
                    for loc_k, j_idx in zip(I2_idx, J2_arr[I2_idx]):
                        if j_idx < 0:
                            continue
                        if task.LBout[loc_k] > 0:
                            t_ref_model.b[j_idx, 0] = float(task.LBout[loc_k])
                        t_ref_model.b[j_idx, 1] = float(task.UBout[loc_k])

            ## Add new reactions defined directly in the task
            if getattr(task, 'equations', []):
                rxn = {
                    'equations': task.equations,
                    'lb': task.LBequ.tolist() if hasattr(task.LBequ, 'tolist') else list(task.LBequ),
                    'ub': task.UBequ.tolist() if hasattr(task.UBequ, 'tolist') else list(task.UBequ),
                    'rxns': [f"TEMPORARY_{k+1}" for k in range(len(task.equations))]
                }
                t_model = added_rxns(t_model,rxn)
                t_ref_model = added_rxns(t_ref_model, rxn)
                # extend scores with zeros for thses temp reactions
                t_rxn_scores = np.concatenate([t_rxn_scores, np.zeros(len(rxn['lb']))])

            ## Apply changed bounds to existing reactions
            if getattr(task, 'changed', []):
                t_model = set_param(t_model, 'lb', task.changed, task.LBrxn)
                t_model = set_param(t_model, 'ub', task.changed, task.UBrxn)
                t_ref_model = set_param(t_ref_model, 'lb', task.changed, task.LBrxn)
                t_ref_model = set_param(t_ref_model, 'ub', task.changed, task.UBrxn)

            ## Solve and (if needed) gap-fill
            sol = solve_lp(t_model)
            sol_failed = not hasattr(sol, 'fluxes') or sol.status != 'optimal'

            if sol_failed:
                failed = False
                try:
                    new_rxns, new_model, exit_flag = ftinit_fill_gaps(t_model, model, t_ref_model, False, supress_warnings, t_rxn_scores, params, verbose)
                    if exit_flag == -2:
                        disp_em(f"\"[{task.id}] {task.description}\" was aborted before reaching optimality. Consider increasing params.max_time", is_error=False)
                except Exception:
                    disp_em(f"\"[{task.id}] {task.description}\" could not be performed for any set of reactions", is_error=False)
                    failed=True

                if not failed:
                    if new_rxns:
                        n_added += len(new_rxns)
                        # Print list of added reactions (per MATLAB behaivior)
                        print(f"task: {i+1}")
                        for rid in new_rxns:
                            print(rid)
                        # In MATLAB they previously merged; here we follow the final choice: assign model = new_model
                        model = new_model
                        # Mark added reactions against ref_model list
                        for rid in new_rxns:
                            if rid in ref_rxn_ids:
                                added_rxns[ref_rxn_ids.index(rid), i] = True
                    if print_output:
                        print(f"[{task.id}] {task.description}: Added {len(new_rxns)} reaction(s), {n_added} reactions added in total")
            else:
                if print_output:
                    print(f"[{task.id}] {task.description}: Added 0 reaction(s), {n_added} reactions added in total")

            supress_warnings =True

            ## Print fluxes if requested
            if getattr(task, 'print_fluxes', False) and print_output:
                if not sol_failed:
                    sol = solve_lp(t_model, 1) # mazimize
                    print_fluxes(t_model, sol.fluxes)
                    print()
                else:
                    # Use gap-filling model if available
                    try:
                        sol2 = solve_lp(model, 1)
                        print_fluxes(model, sol2.fluxes)
                        print()
                    except Exception:
                        pass
            # Reset temporary model to current base model (which may have changed via gap-filling)
            t_model = copy.deepcopy(model)
            # Update model_mets (new mets may have been introduced by gap-filled reactions)
            model_mets = [_met_key(m) for m in model.metabolites]
        else:
            disp_em(f"\"[{task.id}] {task.description}\ is set as SHOULD FAIL. Such tasks cannot be modelled using this approach and the task is therefore ignored", is_error=False)
    out_model = Model
    return out_model, added_rxns

