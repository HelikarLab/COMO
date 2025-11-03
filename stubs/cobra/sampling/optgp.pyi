import pandas as pd
from .hr_sampler import HRSampler
from _typeshed import Incomplete
from cobra import Model

__all__ = ['OptGPSampler']

class OptGPSampler(HRSampler):
    processes: Incomplete
    center: Incomplete
    def __init__(self, model: Model, thinning: int = 100, processes: int | None = None, nproj: int | None = None, seed: int | None = None, **kwargs) -> None: ...
    def sample(self, n: int, fluxes: bool = True) -> pd.DataFrame: ...
