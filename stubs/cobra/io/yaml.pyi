from .dict import model_from_dict as model_from_dict, model_to_dict as model_to_dict
from _typeshed import Incomplete
from cobra import Model as Model
from io import TextIOBase
from pathlib import Path
from ruamel.yaml.main import YAML
from typing import Any

YAML_SPEC: str

class CobraYAML(YAML):
    def dump(self, data: dict, stream: TextIOBase | None = None, **kwargs: Any) -> str: ...

yaml: Incomplete

def to_yaml(model: Model, sort: bool = False, **kwargs: Any) -> str: ...
def from_yaml(document: str) -> Model: ...
def save_yaml_model(model: Model, filename: str, sort: bool = False, **kwargs: Any) -> None: ...
def load_yaml_model(filename: str | Path) -> Model: ...
