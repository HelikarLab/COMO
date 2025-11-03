import pandas as pd
from _typeshed import Incomplete
from cobra import Metabolite, Model, Reaction
from typing import Iterable

__all__ = ['Solution', 'get_solution']

class Solution:
    objective_value: Incomplete
    status: Incomplete
    fluxes: Incomplete
    reduced_costs: Incomplete
    shadow_prices: Incomplete
    def __init__(self, objective_value: float, status: str, fluxes: pd.Series, reduced_costs: pd.Series | None = None, shadow_prices: pd.Series | None = None, **kwargs) -> None: ...
    def __getitem__(self, reaction_id: str) -> float: ...
    get_primal_by_id = __getitem__
    def to_frame(self) -> pd.DataFrame: ...

def get_solution(model: Model, reactions: Iterable['Reaction'] | None = None, metabolites: Iterable['Metabolite'] | None = None, raise_error: bool = False) -> Solution: ...
