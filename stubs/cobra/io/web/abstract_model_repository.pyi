import abc
import httpx
from abc import ABC, abstractmethod

class AbstractModelRepository(ABC, metaclass=abc.ABCMeta):
    name: str
    def __init__(self, *, url: httpx.URL | str, **kwargs) -> None: ...
    @property
    def url(self) -> httpx.URL: ...
    @abstractmethod
    def get_sbml(self, model_id: str) -> bytes: ...
