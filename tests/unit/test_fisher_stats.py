import cobra
import numpy as np
from como.stats.fisher_exact_test import FisherExactTest


def test_fisher_stats():
    reference_model = cobra.io.load_matlab_model("main/data/reference_models/GeneralModelUpdatedV3.mat")
    scenario_model = cobra.io.read_sbml_model("tests/inputs/naiveB_model.xml")
    real = FisherExactTest.run(reference=reference_model, scenario=scenario_model, pathway="Glycolysis/gluconeogenesis")

    assert real.statistic == np.float64(4.321708185053381)
    assert real.pvalue == np.float64(1.2883495211648955e-05)
    assert real.a == 32
    assert real.b == 10
    assert real.c == 4496
    assert real.d == 6072
