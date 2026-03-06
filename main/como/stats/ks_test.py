from dataclasses import dataclass
from typing import ClassVar, Literal

import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy.stats import ks_2samp

from como.stats._two_sample import KS_RESULT, T_ALTERNATIVE, BaseTwoSample

__all__ = ["KSTest"]


@dataclass(frozen=True, kw_only=True, slots=True)
class KSTest(BaseTwoSample[KS_RESULT]):
    _fields: ClassVar[dict[str, type]] = {
        "statistic": float,
        "pvalue": float,
        "statistic_location": float,
        "statistic_sign": int,
    }

    reaction_ids: list[str]
    statistic: npt.NDArray[np.floating]
    pvalue: npt.NDArray[np.floating]
    statistic_location: npt.NDArray[np.floating]
    statistic_sign: npt.NDArray[np.integer]

    @staticmethod
    def _worker(a: npt.NDArray[np.floating], b: npt.NDArray[np.floating], **kwargs) -> KS_RESULT:
        """Calculate the KS statistic.

        Args:
            a: First array
            b: Second array
            kwargs: Additional keyword arguments to pass to `ks_2samp`

        Returns:
            A tuple of (statistic, pvalue, statistic_location, statistic_sign)
        """
        res = ks_2samp(a, b, **kwargs)
        return res.statistic, res.pvalue, res.statistic_location, res.statistic_sign

    @classmethod
    def run(
        cls,
        df1: pd.DataFrame,
        df2: pd.DataFrame,
        alternative: T_ALTERNATIVE = "two-sided",
        method: Literal["auto", "exact", "asymp"] = "auto",
        axis: int = 0,
        nan_policy: Literal["raise", "propagate", "omit"] = "propagate",
        keepdims: bool = False,
        cores: int = 1,
    ) -> "KSTest":
        """Run the KS test on two dataframes.

        Args:
            df1: The first dataframe to process; obtained from running `cobra.sampling.sample`.
                Columns should be reaction IDs and rows should be samples.
            df2: The second dataframe to process; obtained from running `cobra.sampling.sample`.
                Columns should be reaction IDs and rows should be samples.
            alternative: The alternative hypothesis to test.
            method: The method to use for calculating the p-value.
            axis: The axis to perform the test along.
            nan_policy: The policy to use for handling NaNs.
            keepdims: Whether to keep the dimensions of the input arrays.
            cores: The number of CPU cores to use for multiprocessing.

        Returns:
            An instance of `KSTest` containing the results of the test.
        """
        all_reactions, results = cls._run(
            df1=df1,
            df2=df2,
            cores=cores,
            worker_kwargs={"alternative": alternative, "method": method, "axis": axis, "nan_policy": nan_policy, "keepdims": keepdims},
        )
        return cls(
            reaction_ids=all_reactions,
            statistic=results["statistic"].astype(float),
            pvalue=results["pvalue"].astype(float),
            statistic_location=results["statistic_location"].astype(float),
            statistic_sign=results["statistic_sign"].astype(int),
        )

    @property
    def df(self) -> pd.DataFrame:
        """DataFrame representation of the results.

        Returns:
            A DataFrame with columns "statistic", "pvalue", "statistic_location", and "statistic_sign".
        """
        return pd.DataFrame(
            {
                "statistic": self.statistic,
                "pvalue": self.pvalue,
                "statistic_location": self.statistic_location,
                "statistic_sign": self.statistic_sign,
            },
            index=pd.Index(name="reaction_id", data=self.reaction_ids),
        )
