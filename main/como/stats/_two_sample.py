from abc import ABC, abstractmethod
from collections.abc import Mapping
from concurrent.futures import ProcessPoolExecutor, as_completed
from typing import ClassVar, Generic, Literal, TypeVar

import numpy as np
import numpy.typing as npt
import pandas as pd

T_BASE_SAMPLE = TypeVar("T_BASE_SAMPLE", bound="BaseTwoSample")
T_ALTERNATIVE = Literal["greater", "less", "two-sided"]
KS_RESULT = tuple[np.floating, np.floating, np.floating, np.int8]
MW_RESULT = tuple[np.floating, np.floating]
TEST_RESULT = TypeVar("TEST_RESULT", KS_RESULT, MW_RESULT)

__all__ = ["BaseTwoSample"]


class BaseTwoSample(ABC, Generic[TEST_RESULT]):
    _fields: ClassVar[dict[str, type]]

    @staticmethod
    @abstractmethod
    def _worker(a: npt.NDArray[np.floating], b: npt.NDArray[np.floating], **kwargs) -> TEST_RESULT: ...

    @property
    @abstractmethod
    def df(self) -> pd.DataFrame:
        """DataFrame representation of the results.

        Returns:
            A DataFrame with columns corresponding to the fields in `_fields`.
        """
        ...

    @classmethod
    def _run(
        cls: type[T_BASE_SAMPLE],
        df1: pd.DataFrame,
        df2: pd.DataFrame,
        cores: int = 1,
        worker_kwargs: dict | None = None,
    ) -> tuple[list[str], Mapping[str, npt.NDArray[np.float64 | np.uint8]]]:
        all_reactions = list(set(df1.columns) & set(df2.columns))
        array_a = df1[all_reactions].to_numpy(dtype=np.float64, copy=False)
        array_b = df2[all_reactions].to_numpy(dtype=np.float64, copy=False)
        n = len(all_reactions)

        results = {field: np.empty(n, dtype=np.dtype(dtype)) for field, dtype in cls._fields.items()}

        with ProcessPoolExecutor(max_workers=cores) as pool:
            futures = {pool.submit(cls._worker, array_a[:, i], array_b[:, i], **(worker_kwargs or {})): i for i in range(n)}
            for future in as_completed(futures):
                col_idx: int = futures[future]
                res: KS_RESULT | MW_RESULT = future.result()

                for (field, _), value in zip(cls._fields.items(), res, strict=True):
                    results[field][col_idx] = value

        return all_reactions, results
