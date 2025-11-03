import numpy as np
from .hr_sampler import HRSampler as HRSampler
from _typeshed import Incomplete

logger: Incomplete
MAX_TRIES: int

def step(sampler: HRSampler, x: np.ndarray, delta: np.ndarray, fraction: float | None = None, tries: int = 0) -> np.ndarray: ...
