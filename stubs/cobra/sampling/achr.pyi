import pandas as pd
from .core import step as step
from .hr_sampler import HRSampler as HRSampler
from _typeshed import Incomplete
from cobra import Model as Model

class ACHRSampler(HRSampler):
    prev: Incomplete
    def __init__(self, model: Model, thinning: int = 100, nproj: int | None = None, seed: int | None = None, **kwargs) -> None: ...
    def sample(self, n: int, fluxes: bool = True) -> pd.DataFrame: ...
