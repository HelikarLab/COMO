import pandas as pd
from _typeshed import Incomplete
from cobra.core import Metabolite as Metabolite, Model as Model, Reaction as Reaction, Solution as Solution
from cobra.flux_analysis import flux_variability_analysis as flux_variability_analysis, pfba as pfba
from cobra.summary import Summary as Summary
from cobra.util.solver import linear_reaction_coefficients as linear_reaction_coefficients

logger: Incomplete

class ModelSummary(Summary):
    uptake_flux: pd.DataFrame | None
    secretion_flux: pd.DataFrame | None
    def __init__(self, *, model: Model, solution: Solution | None = None, fva: float | pd.DataFrame | None = None, **kwargs) -> None: ...
    def to_string(self, names: bool = False, element: str = 'C', threshold: float | None = None, float_format: str = '.4G', column_width: int = 79) -> str: ...
    def to_html(self, names: bool = False, element: str = 'C', threshold: float | None = None, float_format: str = '.4G') -> str: ...
