import pathlib
import types
from numbers import Number

from _typeshed import Incomplete

from .singleton import Singleton

__all__ = ["Configuration"]


class Configuration(metaclass=Singleton):
    tolerance: float
    lower_bound: Incomplete
    upper_bound: Incomplete
    processes: Incomplete
    max_cache_size: Incomplete
    cache_expiration: Incomplete

    def __init__(self, **kwargs) -> None: ...
    @property
    def solver(self) -> types.ModuleType: ...
    @solver.setter
    def solver(self, value) -> None: ...
    @property
    def bounds(self) -> tuple[Number | None, Number | None]: ...
    @bounds.setter
    def bounds(self, bounds: tuple[Number | None, Number | None]) -> None: ...
    @property
    def cache_directory(self) -> pathlib.Path: ...
    @cache_directory.setter
    def cache_directory(self, path: pathlib.Path | str) -> None: ...
