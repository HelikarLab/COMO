# This function is required in the function in getinitsteps.py (used in ftINIT) originally written in MATLAB.
# Source code: https://github.com/SysBioChalmers/RAVEN/blob/cd82b5c6ba26a88bed2a9ca179c682a6434fd541/INIT/INITStepDesc.m

import numpy as np

"""
Need to be reworked. Learn about defining class in MATLAB and translating that into Python.
"""


class INITStepsDesc:
    def __init__(
        self,
        posRevOff_=False,
        allowMetSecr_=False,
        howToUsePrevResults_="essential",
        rxnsToIgnoreMasks_=None,
        metsToIgnore_=None,
        MILPParams_=None,
        absMIPGaps_=10,
    ):
        self.PosRevOff_ = posRevOff_
        self.AllowMetSecr = allowMetSecr_
        self.HowToUsePrevResults_ = howToUsePrevResults_

        if rxnsToIgnoreMasks_ is not None:
            self.RxnsToIgnoreMasks_ = rxnsToIgnoreMasks_
        else:
            self.RxnsToIgnoreMasks_ = np.array([[1], [0], [0], [0], [0], [0], [0], [0]])

        if metsToIgnore_ is not None:
            self.MetsToIgnore = metsToIgnore_
        else:
            self.MetsToIgnore = np.array([[1], [0], [0], [0], [0], [0], [0], [0]])

        if MILPParams_ is not None:
            self.MILPParams = MILPParams_
        else:
            params = {"MIPGap": 0.0004, "TimeLimit": 5000}
            self.MILPParams = [params]  # mimic cell array of struct

        self.AbsMIPGaps = absMIPGaps_
