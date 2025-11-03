import pandas as pd
from _typeshed import Incomplete
from cobra import Model as Model, Reaction as Reaction, Solution as Solution
from cobra.flux_analysis import flux_variability_analysis as flux_variability_analysis, pfba as pfba
from cobra.summary import Summary as Summary

logger: Incomplete

class ReactionSummary(Summary):
    def __init__(self, *, reaction: Reaction, model: Model, solution: Solution | None = None, fva: float | pd.DataFrame | None = None, **kwargs) -> None: ...
    def to_string(self, names: bool = False, threshold: float | None = None, float_format: str = '.4G', column_width: int = 79) -> str: ...
    def to_html(self, names: bool = False, threshold: float | None = None, float_format: str = '.4G') -> str: ...
