import pandas as pd
from .boundary_types import find_boundary_types as find_boundary_types
from _typeshed import Incomplete
from cobra import Model as Model, Reaction as Reaction

logger: Incomplete

def add_linear_obj(model: Model) -> None: ...
def add_mip_obj(model: Model) -> None: ...
def minimal_medium(model: Model, min_objective_value: float = 0.1, exports: bool = False, minimize_components: bool | int = False, open_exchanges: bool = False) -> pd.Series | pd.DataFrame | None: ...
