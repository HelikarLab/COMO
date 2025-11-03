from dataclasses import dataclass
from typing import ClassVar, Literal

import numpy as np
import numpy.typing as npt
import pandas as pd
from scipy.stats import PermutationMethod, mannwhitneyu

from como.stats._two_sample import MW_RESULT, T_ALTERNATIVE, BaseTwoSample

__all__ = ["MannWhitneyUTest"]


@dataclass(frozen=True, kw_only=True, slots=True)
class MannWhitneyUTest(BaseTwoSample[MW_RESULT]):
    _fields: ClassVar[dict[str, type]] = {"statistic": float, "pvalue": float}

    reaction_ids: list[str]
    statistic: npt.NDArray[float]
    pvalue: npt.NDArray[float]

    @staticmethod
    def _worker(a: npt.NDArray[float], b: npt.NDArray[float], **kwargs) -> MW_RESULT:
        """Calculate the MWU statistic.

        Args:
            a: First array
            b: Second array
            kwargs: Additional keyword arguments to pass to `mannwhitneyu`

        Returns:
            A tuple of (statistic, pvalue)
        """
        res = mannwhitneyu(x=a, y=b, **kwargs)
        return float(res.statistic), float(res.pvalue)

    @classmethod
    def run(
        cls,
        df1: pd.DataFrame,
        df2: pd.DataFrame,
        alternative: T_ALTERNATIVE = "two-sided",
        use_continuity: bool = True,
        axis: int = 0,
        method: Literal["auto", "asymptotic", "exact"] | PermutationMethod = "auto",
        cores: int = 1,
    ) -> "MannWhitneyUTest":
        """Run the MWU test on two dataframes.

        Args:
            df1: The first dataframe to process; obtained from running `cobra.sampling.sample`.
                Columns should be reaction IDs and rows should be samples.
            df2: The second dataframe to process; obtained from running `cobra.sampling.sample`.
                Columns should be reaction IDs and rows should be samples.
            alternative: The alternative hypothesis to test.
            use_continuity: Whether to apply a continuity correction when using the asymptotic method.
            axis: The axis to perform the test along.
            method: The method to use for calculating the p-value.
            cores: The number of CPU cores to use for multiprocessing.

        Returns:
            An instance of `MannWhitneyUTest` containing the results of the test.
        """
        all_reactions, results = cls._run(
            df1=df1,
            df2=df2,
            cores=cores,
            worker_kwargs={"alternative": alternative, "use_continuity": use_continuity, "axis": axis, "method": method},
        )
        return cls(reaction_ids=all_reactions, statistic=results["statistic"].astype(float), pvalue=results["pvalue"].astype(float))

    @property
    def df(self) -> pd.DataFrame:
        """DataFrame representation of the results.

        Returns:
            A DataFrame with columns "statistic" and "pvalue".
        """
        return pd.DataFrame(
            {"statistic": self.statistic, "pvalue": self.pvalue},
            index=pd.Index(name="reaction_id", data=self.reaction_ids),
        )
