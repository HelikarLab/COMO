import abc
import numpy as np
import pandas as pd
from _typeshed import Incomplete
from abc import ABC, abstractmethod
from cobra import Model as Model
from cobra.util import constraint_matrices as constraint_matrices, create_stoichiometric_matrix as create_stoichiometric_matrix, nullspace as nullspace
from typing import NamedTuple

logger: Incomplete

class Problem(NamedTuple):
    equalities: np.ndarray
    b: np.ndarray
    inequalities: np.ndarray
    bounds: np.ndarray
    variable_fixed: np.ndarray
    variable_bounds: np.ndarray
    nullspace: np.matrix
    homogeneous: bool

def shared_np_array(shape: tuple[int, int], data: np.ndarray | None = None, integer: bool = False) -> np.ndarray: ...

class HRSampler(ABC, metaclass=abc.ABCMeta):
    model: Incomplete
    feasibility_tol: Incomplete
    bounds_tol: Incomplete
    thinning: Incomplete
    nproj: Incomplete
    n_samples: int
    retries: int
    problem: Incomplete
    fwd_idx: Incomplete
    rev_idx: Incomplete
    warmup: Incomplete
    def __init__(self, model: Model, thinning: int, nproj: int | None = None, seed: int | None = None, **kwargs) -> None: ...
    n_warmup: int
    def generate_fva_warmup(self) -> None: ...
    @abstractmethod
    def sample(self, n: int, fluxes: bool = True) -> pd.DataFrame: ...
    def batch(self, batch_size: int, batch_num: int, fluxes: bool = True) -> pd.DataFrame: ...
    def validate(self, samples: np.matrix) -> np.ndarray: ...
