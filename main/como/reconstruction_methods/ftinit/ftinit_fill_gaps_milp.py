# This is the Python version of ftINITFillGapsMILP (used in ftinit_fill_gaps) originally written in MATLAB.
# Source code: https://github.com/SysBioChalmers/RAVEN/blob/cd82b5c6ba26a88bed2a9ca179c682a6434fd541/INIT/ftINITFillGapsMILP.m

# ftinit_fill_gaps_milp: Returns the minimal set of fluxes that satisfy the model using mixed integer linear programming. This is an optimized
# variant of the old function "getMinNrFluxes" that is adapted to ftINIT. It deos not need to make an irrev model, which takes time. The problem
# also becomes smaller (fewer integers but larger matrix). Only tested with Gurobi.

## INPUT
# model
# to_minimize:
# params:
# scores:

## OUTPUT
# x
# I
# exit_flag
