from .parsimonious import add_pfba as add_pfba
from .variability import flux_variability_analysis as flux_variability_analysis
from _typeshed import Incomplete
from cobra import Model as Model, Solution as Solution

logger: Incomplete

def geometric_fba(model: Model, epsilon: float = 1e-06, max_tries: int = 200, processes: int | None = None) -> Solution: ...
