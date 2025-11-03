from .helpers import normalize_cutoff as normalize_cutoff
from _typeshed import Incomplete
from cobra.core import Model as Model, Reaction as Reaction

logger: Incomplete
LARGE_VALUE: float

def fastcc(model: Model, flux_threshold: float = 1.0, zero_cutoff: float | None = None) -> Model: ...
