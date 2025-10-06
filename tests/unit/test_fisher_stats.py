import contextlib
from pathlib import Path

import cobra
import numpy as np
from como.stats.fisher_exact_test import FisherExactTest


def test_fisher_stats():
    with (
        contextlib.redirect_stderr(Path("/dev/null").open("w", encoding="utf-8")),
        contextlib.redirect_stdout(Path("/dev/null").open("w", encoding="utf-8")),
    ):
        reference_model = cobra.io.load_matlab_model("main/data/reference_models/GeneralModelUpdatedV3.mat")
        scenario_model = cobra.io.read_sbml_model("tests/inputs/naiveB_model.xml")
    real = FisherExactTest.run(reference=reference_model, scenario=scenario_model, pathway="Glycolysis/gluconeogenesis")
    expected = FisherExactTest(
        pathway="Glycolysis/gluconeogenesis",
        statistic=np.float64(4.321708185053381),
        pvalue=np.float64(1.2883495211648955e-05),
        a=32,
        b=10,
        c=4496,
        d=6072,
    )

    assert real == expected


test_fisher_stats()
