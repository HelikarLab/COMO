# This is the Python version of getinitsteps (used in ftinit) originally written in MATLAB.
# Source code: https://github.com/SysBioChalmers/RAVEN/blob/cd82b5c6ba26a88bed2a9ca179c682a6434fd541/INIT/getINITSteps.m

# getINITSteps
#   Converts a reaction score to the gene expression (CPM or TPM) required
#   to get that reaction score, if the GPR is only a single gene.
#   Useful function primarily in test cases, where you want to be able to
#   define the reaction scores of rxns, but need to send in gene expression.
#   Note that all combinations of steps will not work. In general, avoid 'exclude'
#   if you want to define new ways to run the algorithm.
#
#   metsToIgnore  Structure describing mets that can be removed from the model
#                 before running ftINIT, such as water etc.
#                 (optional, default [])
#       simpleMets
#           mets  Names of metabolites to remove
#           compsToKeep Compartments for which metabolites should be kept.
#   series        Describes the way to run ftINIT:
#                 '1+1'          Standard behavior. Step 1 and 2 described in
#                                the paper are merged into 1.
#                 '2+1'          The 3-step procedure described in the paper.
#                                Faster and slightly less accurate than '1+1 steps'
#                 '1+0'          Same as '1+1 steps', but skips step 3 in described
#                                in the paper. This will result in a model including
#                                a lot of reactions without GPRs. It is particularly
#                                useful for structural comparison, since the reactions
#                                removed in step 3 may be a bit random and doesn't
#                                really add any extra information. Faster than
#                                '1+1 steps'
#                 '2+0'          Same as '2+1 steps', but skips step 3 in described
#                                in the paper. Faster and slightly less accurate
#                                than '1+0 steps', but will yield similar results.
#                 'full'         1-step run - similar to the old tINIT version but
#                                without simplifications. Accurate, but very slow.
#                                This is mainly used for testing purposes.
#                 (optional, default '1+1')
#
#   steps         Cell array of steps, used as input to ftINIT
#
# Usage: steps = getINITSteps(metsToIgnore, series)
###----------------------------------------------------CODE----------------------------------------------------###
import inspect

from reconstruction_methods.ftinit.initstepdesc import INITStepDesc


def get_initstep(metsToIgnore, series):
    sig = inspect.signature(get_initstep)
    num_args = len(sig.parameters)

    if num_args < 1:
        metsToIgnore = []
    if num_args < 2:
        series = "1+1"

    if str(series) == series:
        params1 = {"MIPGap": 0.0004, "TimeLimit": 120}
        params2 = {"MIPGap": 0.0030, "TimeLimit": 5000}
        params = [params1, params2]
        params3 = {"MIPGap": 0.0004, "TimeLimit": 5}
        paramsStep3 = [params3, params1, params2]

        # The paramsStep3 involves a quick first run . The objective value is often small in the third step (~800),
        # and 0.0004 of that is a very small number.
        # With this first step, the rough value of the objective function will be estimated, which will generate
        # an absolute MIPGap limit that is much larger for the second iteration.
        steps = [
            INITStepDesc(False, False, "ignore", [1, 1, 1, 1, 1, 1, 1, 0], metsToIgnore, params, [[10], [20]]),
            INITStepDesc(False, False, "essential", [1, 0, 0, 0, 1, 0, 0, 0], metsToIgnore, paramsStep3, [[10], [10], [20]]),
        ]
        return steps

    elif str(series) == "2+1":
        params1 = {"MIPGap": 0.0004, "TimeLimit": 120}
        params2 = {"MIPGap": 0.0030, "TimeLimit": 5000}
        params = [params1, params2]
        params3 = {"MIPGap": 0.0004, "TimeLimit": 5}
        paramsStep3 = [params3, params1, params2]

        # The paramsStep3 involves a quick first run. The objective value is often small in the third step (~800), and 0.0004
        # of the test is a very small number. With this first step, the rough value of the objective function will be estimated,
        # which will generate an absolute MIPGap limit that is much larger for the second interation.

        steps = [
            INITStepDesc(True, True, "ignore", [1, 1, 1, 1, 1, 1, 1, 0], metsToIgnore, params, [[10], [20]]),
            INITStepDesc(False, False, "essential", [1, 1, 1, 1, 1, 1, 1, 0], metsToIgnore, params, [[10], [20]]),
            INITStepDesc(False, False, "essential", [1, 0, 0, 0, 1, 0, 0, 0], metsToIgnore, paramsStep3, [[10], [10], [20]]),
        ]
        return steps

    elif str(series) == "1+0":  # joins step 1 and 2, skips step3
        params1 = {"MIPGap": 0.0004, "TimeLimit": 120}
        params2 = {"MIPGap": 0.0030, "TimeLimit": 5000}
        params = [params1, params2]
        steps = [INITStepDesc(False, False, "ignore", [1, 1, 1, 1, 1, 1, 1, 0], metsToIgnore, params, [[10], [20]])]
        return steps

    elif str(series) == "2+0":  # Skips step 3
        params1 = {"MIPGap": 0.0004, "TimeLimit": 120}
        params2 = {"MIPGap": 0.0030, "TimeLimit": 5000}
        params = [params1, params2]
        steps = [
            INITStepDesc(True, True, "ignore", [1, 1, 1, 1, 1, 1, 1, 0], metsToIgnore, params, [[10], [20]]),
            INITStepDesc(False, False, "essential", [1, 1, 1, 1, 1, 1, 1, 0], metsToIgnore, params, [[10], [20]]),
        ]
        return steps
    elif str(series) == "full":  # Just one run, slow on large models, but this is the 'perfect' setup
        params1 = {"MIPGap": 0.0004, "TimeLimit": 10000}
        params = [params1]
        steps = [INITStepDesc(False, False, "ignore", [0, 0, 0, 0, 0, 0, 0, 0], [], params)]

    else:
        print(f"Invalid series in getINITSteps: {series}")  # not sure what dispEM does exactly (not function found so far in the original code
