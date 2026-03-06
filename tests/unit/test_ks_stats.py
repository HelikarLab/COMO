import numpy as np
import pandas as pd
import pytest

from como.stats.ks_test import KSTest


@pytest.fixture
def expected_result() -> pd.DataFrame:
    df = pd.read_csv("tests/inputs/expected_ks_results.csv", index_col=0)
    df["statistic_sign"] = df["statistic_sign"].astype(int)
    return df


@pytest.mark.parametrize("cores", [1, 2, 4])
def test_ks_stats(expected_result: pd.DataFrame, cores: int):
    seed = 123456789
    size = 1_000
    cols = 10
    gen = np.random.Generator(np.random.PCG64DXSM(seed=seed))

    df1 = pd.DataFrame(
        index=list(range(size)),
        columns=list(range(cols)),
        data=gen.normal(loc=0, size=size * cols).reshape(size, cols).astype(float),
    )
    df2 = pd.DataFrame(
        index=list(range(size)),
        columns=list(range(cols)),
        data=gen.normal(loc=1, size=size * cols).reshape(size, cols).astype(float),
    )

    real_result = KSTest.run(df1, df2, cores=cores)
    real_df = real_result.df
    # real_df.to_csv("/Users/joshl/Projects/COMO/tests/inputs/expected_ks_results.csv")
    # assert False
    assert len(expected_result.columns) == len(real_df.columns)
    assert len(expected_result) == len(real_df)
    assert all(col in real_df.columns for col in expected_result.columns)
    pd.testing.assert_frame_equal(expected_result, real_df)
