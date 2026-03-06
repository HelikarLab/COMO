import pydantic
from .abstract_model_repository import AbstractModelRepository as AbstractModelRepository

class BioModelsFile(pydantic.BaseModel):
    name: str
    size: int

class BioModelsFilesResponse(pydantic.BaseModel):
    main: list[BioModelsFile]

class BioModels(AbstractModelRepository):
    name: str
    def __init__(self, **kwargs) -> None: ...
    def get_sbml(self, model_id: str) -> bytes: ...
