from ...core import Configuration as Configuration
from .abstract_model_repository import AbstractModelRepository as AbstractModelRepository
from .bigg_models_repository import BiGGModels as BiGGModels
from .biomodels_repository import BioModels as BioModels
from .cobrapy_repository import Cobrapy as Cobrapy
from _typeshed import Incomplete
from cobra.core import Model as Model
from typing import Iterable

logger: Incomplete
configuration: Incomplete
DEFAULT_REPOSITORIES: Incomplete

def load_model(model_id: str, repositories: Iterable[AbstractModelRepository] = ..., cache: bool = True) -> Model: ...
def get_model_from_gzip_sbml(stream: bytes) -> Model: ...
