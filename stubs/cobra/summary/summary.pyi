import abc
import pandas as pd
from _typeshed import Incomplete
from abc import ABC, abstractmethod
from cobra import Model as Model, Solution as Solution

logger: Incomplete

class Summary(ABC, metaclass=abc.ABCMeta):
    def __init__(self, **kwargs) -> None: ...
    @property
    def tolerance(self) -> float: ...
    @abstractmethod
    def to_string(self, names: bool = False, threshold: float | None = None, float_format: str = '.4G', column_width: int = 79) -> str: ...
    @abstractmethod
    def to_html(self, names: bool = False, threshold: float | None = None, float_format: str = '.4G') -> str: ...
    def to_frame(self) -> pd.DataFrame: ...
