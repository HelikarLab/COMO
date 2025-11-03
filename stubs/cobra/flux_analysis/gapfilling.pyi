from ..core import Model as Model
from ..util import fix_objective_as_constraint as fix_objective_as_constraint, interface_to_str as interface_to_str
from _typeshed import Incomplete
from cobra import Reaction as Reaction

logger: Incomplete

class GapFiller:
    original_model: Incomplete
    lower_bound: Incomplete
    model: Incomplete
    integer_threshold: Incomplete
    universal: Incomplete
    penalties: Incomplete
    indicators: Incomplete
    costs: Incomplete
    def __init__(self, model: Model, universal: Model | None = None, lower_bound: float = 0.05, penalties: dict[str, int] | dict['Reaction', int] | None = None, exchange_reactions: bool = False, demand_reactions: bool = True, integer_threshold: float = 1e-06, **kwargs) -> None: ...
    def extend_model(self, exchange_reactions: bool = False, demand_reactions: bool = True) -> None: ...
    def update_costs(self) -> None: ...
    def add_switches_and_objective(self) -> None: ...
    def fill(self, iterations: int = 1) -> list[list['Reaction']]: ...
    def validate(self, reactions: list['Reaction']) -> bool: ...

def gapfill(model: Model, universal: Model | None = None, lower_bound: float = 0.05, penalties: dict[str, 'Reaction'] | None = None, demand_reactions: bool = True, exchange_reactions: bool = False, iterations: int = 1): ...
