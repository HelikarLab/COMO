from .object import Object as Object
from _typeshed import Incomplete

element_re: Incomplete

class Formula(Object):
    formula: Incomplete
    elements: Incomplete
    def __init__(self, formula: str | None = None, **kwargs) -> None: ...
    def __add__(self, other_formula: Formula | str) -> Formula: ...
    def parse_composition(self) -> None: ...
    @property
    def weight(self) -> float: ...

elements_and_molecular_weights: Incomplete
