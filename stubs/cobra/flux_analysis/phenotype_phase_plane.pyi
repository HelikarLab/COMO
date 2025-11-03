import pandas as pd
from ..exceptions import OptimizationError as OptimizationError
from .helpers import normalize_cutoff as normalize_cutoff
from cobra import Model as Model, Reaction as Reaction
from optlang.interface import Objective as Objective

def production_envelope(model: Model, reactions: list['Reaction'], objective: dict | Objective | None = None, carbon_sources: list['Reaction'] | None = None, points: int = 20, threshold: float | None = None) -> pd.DataFrame: ...
