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

import numpy as np
from mpmath.libmp import mpi_delta


def ftinit_fill_gaps_for_all_tasks(model, ref_model, input_file, print_output, rxn_scores, task_structure, params, verbose):
    if rxn_scores is None or len(rxn_scores) == 0:
        rxn_scores = np.full(len(ref_model.reactions), -1.0)
    if not task_structure & input_file.isfile:
        ValueError(f"Task file cannot be found.{input_file}")
    if model.id.lower() == ref_model.id.lower():
        print(
            'NOTE: The model and reference model have the same IDs. The ID for the referece model was set to "refModel" '
            "in order to keep track of the origin of reactions.\n"
        )
    rxn_scores = np.array(rxn_scores)
    if np.any(rxn_scores >= 0):
        ValueError("Only negative values are allowed in rxnScores")
    # Prepare the input models a little
    # REWRITE
    model.b = np.zeros((len(model.metabolites),2))
    model_mets = [f"{met.name.upper()}[{met.compartment.upper()}" for met in model.metabolites]

    # This is the mets in the reference model. used if the tasks involve metabolites that doesn't exist in the model.
    large_model_mets = [f"{met.name.upper()}[{met.compartment.upper()}]" for met in ref_model.metabolites]

    if not model["unconstrained"]:
        print(
            "Exchange metabolites should normally not be removed from the model when using checkTasks. Inputs and outputs are defined in the task file "
            "instead. Use importModel(file, false) to import a model with exchange metabolites remaining"
        )

    if task_structure is None or len(task_structure) == 0:
        task_structure = parseTaskList(input_file)  ## parseTaskList function needed (not found as a MATLAB file)

    t_model = model
    added_rxns = np.zeros((len(ref_model.reactions), len(task_structure)), dtype=bool)
    supress_warnings = False
    n_added = 0

    for i in range(len(task_structure)):
        if task_structure[i].taskType == "SHOULD FAIL":
            t_refmodel = ref_model  # we need to add stuff to this one as well...
            t_rxn_scores = rxn_scores  # these need to be extended (with zeros) when rxns are added in tasks

            # Set the inputs
            if task_structure[i].inputs is None or len(task_structure[i].inputs) == 0:
                pass
                if not np.all(good_mets):
                    # Not all of the inputs could be found in the small model.
                    # Check if they exist in the large model.
                    pass
                    if not np.all(found):
                        pass
                    else:
                        # Otherwise add them to the model
                        """
                        code comes here
                        """
                        # Add the metabolite both to the base model and the model used in the current task
                        """
                        code comes here
                        """
                    # By now the indexes might be getting a bit confusing, but this is to update the indexes of the "real"
                    # metabolites to point to the newly added ones
                    """
                    code comes here
                    """
                # if numel(J(I))...
                # If all metabolites should be added
                if np.any(K):
                    # Check if ALLMETS is the first metabolite. Otherwise print a warning since it will write over any other constraints that are set
                    if K[0] == 0:
                        print(
                            f"ALLMETS is used as an input in [{task_structure[i].id}]{task_structure[i].description} but it is not the first "
                            f"metabolite in the list. Constraints defined for the metabolites before it will be over-written"
                        )
                    # Use the first match of ALLMETS. There should only be one, but still...
                    """
                    code comes here 
                    """
                if np.any(L):
                    L = find(L)
                    for j in range(len(L)):
                        # The compartment defined
                        #compartmen =
                        # Check if it exists in the model
                        # C = find()
                        if np.any(C):
                            # Match to metabolites
                            # t_model.b
                        else:
                            print(f"The compartment defined for ALLMETSIN in [{task_structure[i].id}] {task_structure[i].description} does not exist")
                # Then add the normal constraints
                if np.any(J(I)):
                    pass
                # for the t_refmodel as well
                if np.any(J2(I2)):
                    pass
            if task_structure[i].outputs is None:
                pass
    return out_model, added_rxns
