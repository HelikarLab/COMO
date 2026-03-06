from ..core import GPR as GPR, Gene as Gene, Group as Group, Metabolite as Metabolite, Model as Model, Reaction as Reaction
from ..manipulation.validate import check_metabolite_compartment_formula as check_metabolite_compartment_formula
from ..util.solver import linear_reaction_coefficients as linear_reaction_coefficients, set_objective as set_objective
from _typeshed import Incomplete
from pathlib import Path
from typing import IO, NamedTuple, Pattern

class CobraSBMLError(Exception): ...

LOGGER: Incomplete
config: Incomplete
LOWER_BOUND_ID: str
UPPER_BOUND_ID: str
ZERO_BOUND_ID: str
BOUND_MINUS_INF: str
BOUND_PLUS_INF: str
SBO_FBA_FRAMEWORK: str
SBO_DEFAULT_FLUX_BOUND: str
SBO_FLUX_BOUND: str
SBO_EXCHANGE_REACTION: str
LONG_SHORT_DIRECTION: Incomplete
SHORT_LONG_DIRECTION: Incomplete

class Unit(NamedTuple):
    kind: Incomplete
    scale: Incomplete
    multiplier: Incomplete
    exponent: Incomplete

UNITS_FLUX: Incomplete
SBML_DOT: str
pattern_notes: Pattern
pattern_to_sbml: Pattern
pattern_from_sbml: Pattern
F_GENE: str
F_GENE_REV: str
F_SPECIE: str
F_SPECIE_REV: str
F_REACTION: str
F_REACTION_REV: str
F_GROUP: str
F_GROUP_REV: str
F_REPLACE: dict

def read_sbml_model(filename: str | IO | Path, number: type = ..., f_replace: dict = ..., **kwargs) -> Model: ...
def write_sbml_model(cobra_model: Model, filename: str | IO | Path, f_replace: dict = ..., **kwargs) -> None: ...

URL_IDENTIFIERS_PATTERN: Incomplete
URL_IDENTIFIERS_PREFIX: str
QUALIFIER_TYPES: Incomplete

def validate_sbml_model(filename: str | IO | Path, check_model: bool = True, internal_consistency: bool = True, check_units_consistency: bool = False, check_modeling_practice: bool = False, **kwargs) -> tuple[Model | None, dict]: ...
