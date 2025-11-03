from ..exceptions import OptimizationError as OptimizationError
from ..util.solver import check_solver_status as check_solver_status
from ..util.util import format_long_string as format_long_string
from .formula import elements_and_molecular_weights as elements_and_molecular_weights
from .species import Species as Species
from _typeshed import Incomplete
from cobra.core import Solution as Solution
from cobra.summary.metabolite_summary import MetaboliteSummary as MetaboliteSummary
from optlang.interface import Container as Container
from pandas import DataFrame as DataFrame

element_re: Incomplete

class Metabolite(Species):
    formula: Incomplete
    compartment: Incomplete
    charge: Incomplete
    def __init__(self, id: str | None = None, formula: str | None = None, name: str | None = '', charge: float | None = None, compartment: str | None = None) -> None: ...
    @property
    def constraint(self) -> Container: ...
    @property
    def elements(self) -> dict[str, int | float] | None: ...
    @elements.setter
    def elements(self, elements_dict: dict[str, int | float]) -> None: ...
    @property
    def formula_weight(self) -> int | float: ...
    @property
    def y(self) -> float: ...
    @property
    def shadow_price(self) -> float: ...
    def remove_from_model(self, destructive: bool = False) -> None: ...
    def summary(self, solution: Solution | None = None, fva: float | DataFrame | None = None) -> MetaboliteSummary: ...
