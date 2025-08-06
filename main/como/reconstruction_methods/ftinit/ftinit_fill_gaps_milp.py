# This is the Python version of ftINITFillGapsMILP (used in ftinit_fill_gaps) originally written in MATLAB.
# Source code: https://github.com/SysBioChalmers/RAVEN/blob/cd82b5c6ba26a88bed2a9ca179c682a6434fd541/INIT/ftINITFillGapsMILP.m

# ftinit_fill_gaps_milp: Returns the minimal set of fluxes that satisfy the model using mixed integer linear programming. This is an optimized
# variant of the old function "getMinNrFluxes" that is adapted to ftINIT. It deos not need to make an irrev model, which takes time. The problem
# also becomes smaller (fewer integers but larger matrix). Only tested with Gurobi.

## INPUT
# model: a model structure
# to_minimize: either a cell aray of reaction IDs, a logical vector with the same number of elements as reactions in the model, of a vector of
#              indexes for the reactions that should be minimized (optional, default model.rxns)
# params: psrameter structure as used by get_milp_params (optional)
# scores: vector of weights for the reactions. negative scores should not have flux. positive scores are not possible in this implementation, and
#         they are changed to max(scores(scores<0)). must have the same dimension as to_minimize (find(tominimize) if it is a logical vector)
#         (optional, default -1 for all reactions)
# verbose: if true, the MILP progression will be shown.

## OUTPUT
# x: the corresponding fluxes for the full model
# ind: the indexes of the reactions in to_minimize that were used in the solution
# exit_flag: 1: optimal solution found
#           -1: no feasible solution found
#           -2: optimization time out

# NOTE: Use 1000/mmol/gDW/h as a arbitary large flux. Could possibly cause problems if the fluxes in  the model are larger than that.

import numpy as np
from optlang import Model, Variable, Constraint, Objective
from cobra import Model as cobra_model

def ftinit_fill_gaps_milp(model: Model, to_minimize=None, params=None, scores=None, verbose: bool = False):
    if to_minimize is None:
        to_minimize = model.reactions
    else:
        if isinstance(to_minimize[0], str):
            to_minimize = [model.reactions.get_by_id(rid) for rid in to_minimize]
        elif isinstance(to_minimize[0], int):
            to_minimize = [model.reactions[i] for i in to_minimize]

    if params is None:
        params = {}

    indexes = [model.reactions.index(r) for r in to_minimize]
    rev_flags = np.array([r.reversibility for r in to_minimize])
    irrev_flags = ~rev_flags
    rev_indexes = np.where(rev_flags)[0]
    irrev_indexes = np.where(irrev_flags)[0]

    if scores is None:
        scores = -np.ones(len(to_minimize))
    else:
        scores = np.array(scores)
        if np.any(scores >= 0):
            scores[scores >= 0] = np.max(scores[scores < 0])
        scores = -scores

# The trick to make this possible using a reversible model is that for reversible reactions, we do the following. Set the flux Vi == Vipos - Vineg.
# What happens then isn that if the flux is positive, Vipos will have a nonzero value, and when the flux is negative, Vineg will have a positive
# value. In addition, we can add an arbitrary constant to both Vibos and Vineg. For example, if the flux Vi is -1, Vineg can be 4 and Vipos 3. This
# is however not a problem, because we can be sure that Vineg + Vipos >= abs(Vi). Since we are interested in forcing the ints to be pn when there is
# a flux, it doesn't matter if we overestimate the flux! So, we can simply constrain the boolean Yi to Yi*Maxflux >= Vineg + Vipos.

# The matrix then becomes as this:
#         S       p n int b     var
#         SSSSSSSS        -
#         SSSSSSSS         -
#         SSSSSSSS          -
#         SSSSSSSS           -
#         SSSSSSSS            -
#         SSSSSSSS             -
#         -       1 -
#               -  1 -
#           -          M         -
#             -         M         -
#                 - - M         -
#                  - -   M         -

# An example with 8 rxns and 6 metabolites. - means -1, M max flux, and S is the S matrix. 4 rxns are to be minimized (1,3,5,7) and 1,7 are
# reversible. The p and n are the Vipos and Vineg variables (2 rxns of each). The ints are the Yi for the variables that are to be minimized
# (the rest of the rxns doesn't have any). The mets here are the constraints, so right under the S matrix, you have Vi == Vipos - Vineg for the
# reactions 1 and 7 while the two next rows represent the non-reversible rxns 3 and 5, where we simply say that yi*M >= Vi. The last 2 rows are
# the reactions yi*M >= Vipos + Vineg. To the right, we first have a -I matrix for setting the b constraints, and under that we have rxns that are
# just varibles (rxns) between 0 and Inf to complete the constraints mentioned above.
# Ex: yi*M >= Vipos + Vineg is impl. as yi*M - Vipos - Vineg - var == 0,0 <= var <= Inf.

# All rows should be equal to zero, so we don't set the b vector in the problem
# The reactions should be constrained as follows:
# S - as given in model.lb and model.ub
# pos and neg - between 0 and inf
# ints - between 0 and 1
# b - as stated in the model.b vector - if this has one column, that is lb and ub (fixed value), if two columns, that is lb and ub
# var - between zero and inf

    # Add binary constraints in the following manner: - Add one unique "metabolite" for each integer reaction as a substrate.
    # These metabolites can have net production. Add reactions for the production of each of those metabolites. The amount produced in one reaction unit
    # must be larger than the largest possible flux in the model (but not too large to avoid bad scaling)

    opt_model = Model(name="ftinit_fill_gaps_milp")
    max_flux = 1000

    # Flux variable for all reactions
    flux_vars = {}
    for rxn in model.reactions:
        flux_vars[rxn.id] = Variable(rxn.id, lb=rxn.lower_bound, ub=rxn.upper_bound)
        opt_model.add(flux_vars[rxn.id])

    #Positive and negative flux variables for reversible reactions
    pos_vars = {}
    neg_vars = {}
    for i in rev_indexes:
        rxn = to_minimize[i]
        pos_vars[rxn.id] = Variable(f"{rxn.id}_pos", lb=0)
        neg_vars[rxn.id] = Variable(f"{rxn.id}_neg", lb=0)
        opt_model.add(pos_vars[rxn.id])
        opt_model.add(neg_vars[rxn.id])

    # Integer variables
    int_vars = {}
    for i in range(len(to_minimize)):
        rxn = to_minimize[i]
        int_vars[rxn.id] = Variable(f"y_{rxn.id}", type='binary')
        opt_model.add(int_vars[rxn.id])

    # Auxiliary variables (for MILP formulation)
    aux_vars = {}
    for i in range(len(to_minimize)):
        rxn = to_minimize[i]
        aux_vars[rxn.id] = Variable(f"aux_{rxn.id}", lb=0)
        opt_model.add(aux_vars[rxn.id])

    # Mass balance constraints
    for met in model.metabolites:
        expr = sum(coeff * flux_vars[rxn.id] for rxn, coeff in met._reaction.items())
        opt_model.add(Constraint(expr, lb=0,ub=0,name=f"mass_{met.id}"))

    # Reversible constraints: V = Vpos - Vneg
    for i in rev_indexes:
        rxn = to_minimize[i]
        expr = flux_vars[rxn.id] - pos_vars[rxn.id] + neg_vars[rxn.id]
        opt_model.add(Constraint(expr, lb=0,ub=0,name=f"rev_split{rxn.id}"))

    # Integer activation constraints
    for i in irrev_indexes:
        rxn = to_minimize[i]
        expr = flux_vars[rxn.id] - int_vars[rxn.id] * max_flux
        opt_model.add(Constraint(expr, ub=0, name=f"int_irrev_{rxn.id}"))

    for i in rev_indexes:
        rxn = to_minimize[i]
        expr = pos_vars[rxn.id] + neg_vars[rxn.id] - int_vars[rxn.id] * max_flux
        opt_model.add(Constraint(expr, ub=0, name=f"int_rev_{rxn.id})"))

    # Auxiliary constraint (optional warning mechanism)
    for i in range(len(to_minimize)):
        rxn = to_minimize[i]
        expr = int_vars[rxn.id] * max_flux - flux_vars[rxn.id] - aux_vars[rxn.id]
        opt_model.add(Constraint(expr,lb=0,ub=0,name=f"aux_con_{rxn.id}"))

    # Objective: minimize sum of scores * int_vars
    objective_expr = sum(scores[i] * int_vars[rxn.id] for i, rxn in enumerate(to_minimize))
    opt_model.objective = Objective(objective_expr, direction='min')

    # Solve
    status = opt_model.optimize()

    if status != 'optimal':
        return None, None, -1 if status == 'infeasible' else -2

    # Extract results
    x = np.array([flux_vars[rxn.id].primal for rxn in model.reactions])
    ind = np.array([int_vars[rxn.id].primal > 1e-3 for rxn in to_minimize])

    tmp = np.array([int_vars[rxn.id].primal for rxn in to_minimize])
    intermediate = (tmp > 1e-12) & (tmp < 0.5)

    if intermediate.any():
        print(f"Warning: {np.sum(intermediate)} MILP vars have intermediate values")

    return x,ind,1

