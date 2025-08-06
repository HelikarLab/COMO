# This is the Python version of ftINITFillGaps (used in ftinit_fill_gaps_for_all_tasks) originally written in MATLAB.
# Source code: https://github.com/SysBioChalmers/RAVEN/blob/cd82b5c6ba26a88bed2a9ca179c682a6434fd541/INIT/ftINITFillGaps.m

# ftINITFillGaps: Variant of fillGaps that specially adapted to speed up generation of ftINIT models

## INPUT
# t_model: model that contains that task-specific rxns
# orig_model: model without task-specific rxns
# t_ref_model: reference t_model - the full t_model, containing all rxns used for gap-filling of t_model + the task-specific rxns
# allow_net_production: true if net production of all metabolites is allowed. A reaction can be unable to carry flux because one of the reactants is
#                       unavailable or because one of tbe products can't ve further processed. If this parameter is true, only the first type of
#                       unconnectivity is considered (optional, default false)
# use_model_constraints: true if the constraints specified in the t_model structure should be used. If false, then reactions included from the
#                        template t_model(s) so that as many reactions as possible in t_model can carry flux (optional, default false)
# supress_warnings: false if warnings should be displayed (optional, defaule false)
# rxn_scores: scores for each of the reactions in the reference t_model. The solver will try to maximize the sum of the scores for the
#             included reactions.
# params: Parameter structure as used by get_milp_params
# verbose: if true, the MILP progression will be shown.

## OUTPUT
# added_rxns: the rxns added
# new_model: the t_model with reactions added to fill gaps
# exit_flag: 1: optimal solution found
#           -1: no feasible solution found
#           -2: optimization time out

# This method works by merging the t_model to the reference model and checking which reactions can carry flux. All reactions
# that can't carry flux are removed. If use_model_constraints is false it then solves the MILP problem of minimizing the number of active reactions
# from the reference model that are required to have flux in all the reactions in model. This requires that the input t_model has exchange reactions
# present for the nutrients that are needed for its metabolism. If use_model_constraints is true, then the problem is to include as few reactions as
# possible from the reference models in order to satisfy the t_model constraints. The intended use is that the user can attempt a general gap-filling
# using use_model_constraints=false or a more targeted gap-filling by setting constraints in the model structure and then use
# use_model_constraints=true. Say that the user want to include reactions so that all biomass components can be synthesized. He/She could then define
# a biomass equation and set the lower bound to > 0. Running this function with use_model_constraints=true would then give the smallest set of
# reactions that have to be included in order for the t_model to produce biomass.

from typing import Dict, List, Tuple

import numpy as np
from cobra import Model
from cobra.flux_analysis import pfba
from cobra.util.solver import linear_reaction_coefficients
from sympy import false
from sympy.printing.dot import template


def ftinit_fill_gaps(
    t_model: Model,
    orig_model: Model,
    t_ref_model: Model,
    allow_net_production: bool = False,
    use_model_constraints: bool = False,
    supress_warnings: bool = false,
    rxn_scores: np.ndarray = None,
    params: Dict = None,
    verbose: bool = False,
) -> Tuple[List[str], Model, int]:
    if rxn_scores is None:
        rxn_scores = -1 * np.ones(len(t_ref_model.reactions))
    # Add scores to t_ref_model reactions
    for rxn, score in zip(t_ref_model.reactions, rxn_scores):
        rxn.notes["score"] = score
    t_ref_model = simplify_model(t_ref_model, skip_constrained_rxns=True)

    # Set scores of t_model reactions to 0
    for rxn in t_model.reactions:
        rxn.notes["score"] = 0.0

    # Merge models: t_model + t_ref_model (assumed to be merged here manually if needed)
    full_model = t_ref_model.copy()

    if allow_net_production:
        #         Add larger upper bounds to b vector (RHS of stoichiometric system) COBRApy does not support modifying the b vector directly,
        #         so we handle this logically through bounds
        pass  # Placeholder, assuming flux variability analysis would model this indirectly

    # Zero objective function
    for rxn in full_model.reactions:
        rxn.objective_coefficient = 0.0

    # Identify template (non-t_model) reactions
    template_rxns = [rxn.id for rxn in full_model.reactions if rxn.id not in t_model.reactions]

    # Extract scores for template_rxns
    template_scores = [rxn.notes.get("score", 0) for rxn in full_model.reactions if rxn.id in template_rxns]

    # Solve MILP to identify which template reactions are needed
    selected_rxns, exit_flage = ftinit_fill_gaps_milp(full_model, template_rxns, template_scores, params, verbose)
