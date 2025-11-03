import numpy as np
from ..core import Gene as Gene, Group as Group, Metabolite as Metabolite, Model as Model, Object as Object, Reaction as Reaction
from ..util import create_stoichiometric_matrix as create_stoichiometric_matrix
from ..util.solver import set_objective as set_objective
from _typeshed import Incomplete
from collections import OrderedDict
from pathlib import Path
from typing import IO

logger: Incomplete
MET_MATLAB_TO_PROVIDERS: Incomplete
MET_PROVIDERS_TO_MATLAB: Incomplete
MET_NOTES_TO_MATLAB: Incomplete
MET_MATLAB_TO_NOTES: Incomplete
RXN_MATLAB_TO_PROVIDERS: Incomplete
RXN_PROVIDERS_TO_MATLAB: Incomplete
CONFIDENCE_STR: str
RXN_MATLAB_TO_NOTES: Incomplete
RXN_NOTES_TO_MATLAB: Incomplete
GENE_MATLAB_TO_PROVIDERS: Incomplete
GENE_PROVIDERS_TO_MATLAB: Incomplete
DICT_GENE: str
DICT_GENE_REV: str
DICT_MET: str
DICT_MET_REV: str
DICT_MET_NOTES: str
DICT_MET_NOTES_REV: str
DICT_REACTION: str
DICT_REACTION_REV: str
DICT_REACTION_NOTES: str
DICT_REACTION_NOTES_REV: str
DICT_REPLACE: dict

def load_matlab_model(infile_path: str | Path | IO, variable_name: str | None = None, inf: float = ...) -> Model: ...
def save_matlab_model(model: Model, file_name: str | Path | IO, varname: str | None = None) -> None: ...
def mat_parse_annotations(target_list: list[Object], mat_struct: np.ndarray, d_replace: str = ...) -> None: ...
def mat_parse_notes(target_list: list[Object], mat_struct: np.ndarray, d_replace: str = ...) -> None: ...
def annotations_to_mat(mat_dict: OrderedDict, annotation_list: list[dict], d_replace: str = ...) -> None: ...
def notes_to_mat(mat_dict: OrderedDict, note_list: list[dict], d_replace: str = ...) -> None: ...
def create_mat_dict(model: Model) -> OrderedDict: ...
def from_mat_struct(mat_struct: np.ndarray, model_id: str | None = None, inf: float = ...) -> Model: ...
