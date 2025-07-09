# This function is required in the function in getinitsteps.py (used in ftINIT) originally written in MATLAB.
# Source code: https://github.com/SysBioChalmers/RAVEN/blob/cd82b5c6ba26a88bed2a9ca179c682a6434fd541/INIT/INITStepDesc.m

import numpy as np

"""
Need to be reworked. Learn about defining class in MATLAB and translating that into Python.
"""


class INITStepsDesc:
    def __init__(self, *args, posRevOff_, AllowMetSecr_, howToUsePrevResults_, rxnsToIgnoreMask_, metsToIgnore_, MILPParams_, absMIPGaps_):
        if len(args) > 0:
            self.PosRevOff = posRevOff_
        else:
            self.PosRevOff = False

        if len(args) > 1:
            self.AllowMetSecr = AllowMetSecr_
        else:
            self.AllowMetSecr = False

        if len(args) > 2:
            self.HowToUsePrevResults = howToUsePrevResults_
        else:
            self.HowToUsePrevResults = "essential"

        if len(args) > 3:
            self.RxnsToIgnoreMask = rxnsToIgnoreMask_
        else:
            self.RxnsToIgnoreMask = np.array([[1], [0], [0], [0], [0], [0], [0], [0]])

        if len(args) > 4:
            self.MetsToIgnore = metsToIgnore_
        else:
            self.MetsToIgnore = np.array([[0], [0], [0], [0], [0], [0], [0]])

        if len(args) > 5:
            self.MILPParams = MILPParams_
        else:
            params = {"MIPGap": 0.0004, "TimeLimit": 5000}
            self.MILPParams = params

        if len(args) > 6:
            self.AbsMIPGaps = absMIPGaps_
        else:
            self.AbsMIPGaps = 10
