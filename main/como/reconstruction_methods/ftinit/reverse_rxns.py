# This is the Python version of reverseRxns (used in ftinit) originally written in MATLAB.
# Source code: https://github.com/SysBioChalmers/RAVEN/blob/cd82b5c6ba26a88bed2a9ca179c682a6434fd541/INIT/reverseRxns.m

# reverse_rxns
# model: the model to change (input/output)
# rxns: the rxns to reverse

import numpy as np


def reverse_rxns(model, rxns):
    # Find indices of reactions to reverse
    rxn_ind = np.where(np.isin(model.reactions, rxns))[0]

    if len(rxn_ind) == 0:
        print("Warning: No matching reactions to reverse.")
        return model

    # Reverse the stoichiometric coefficients
    model['S'][:, rxn_ind] = model['S'][:, rxn_ind] * -1

    # Swap the bounds
    ub_temp = model['ub'][rxn_ind].copy()
    model['ub'][rxn_ind] = -model['lb'][rxn_ind]
    model['lb'][rxn_ind] = -ub_temp

    # Flip the objective coefficients
    model['c'][rxn_ind] = -model['c'][rxn_ind]

    return model

