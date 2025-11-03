import pandas as pd
from .achr import ACHRSampler as ACHRSampler
from .optgp import OptGPSampler as OptGPSampler
from cobra import Model as Model

def sample(model: Model, n: int, method: str = 'optgp', thinning: int = 100, processes: int = 1, seed: int | None = None) -> pd.DataFrame: ...
