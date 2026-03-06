from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass, field, fields
from enum import Enum
from pathlib import Path
from typing import ClassVar, NamedTuple, NotRequired, TypedDict

import cobra
import numpy as np
import numpy.typing as npt
import pandas as pd
from loguru import logger

PATH_TYPE = str | Path
LOG_FORMAT = "<green>{time:YYYY-MM-DD HH:mm:ss}</> | <level>{level:<8}</> | <cyan>{name}</>:<cyan>{line}</> - <level>{message}</>"


class AdjustmentMethod(str, Enum):
    """Adjustment method for expression requirement based on differences in number of provided data source types."""

    PROGRESSIVE = "PROGRESSIVE"
    REGRESSIVE = "REGRESSIVE"
    FLAT = "FLAT"
    CUSTOM = "CUSTOM"


class Algorithm(str, Enum):
    GIMME = "GIMME"
    FASTCORE = "FASTCORE"
    IMAT = "IMAT"
    TINIT = "TINIT"


class FilteringTechnique(str, Enum):
    """RNA sequencing filtering capabilities."""

    CPM = "CPM"
    ZFPKM = "ZFPKM"
    TPM = "TPM"
    UMI = "UMI"


class GeneIdentifier(str, Enum):
    ensembl_gene_id = "ensembl_gene_id"
    entrez_gene_id = "entrez_gene_id"
    gene_symbol = "gene_symbol"


class LogLevel(int, Enum):
    TRACE = 5
    DEBUG = 10
    INFO = 20
    SUCCESS = 25
    WARNING = 30
    ERROR = 40
    CRITICAL = 50
    NONE = 100


class RNAType(str, Enum):
    TRNA = "TOTAL"
    MRNA = "MRNA"
    SCRNA = "SCRNA"


class Solver(str, Enum):
    """Solver used to seed context specific model."""

    GLPK = "GLPK"
    GUROBI = "GUROBI"
    SCIPY = "SCIPY"
    GLPK_EXACT = "GLPK_EXACT"


class SourceTypes(str, Enum):
    trna = "trna"
    mrna = "mrna"
    scrna = "scrna"
    proteomics = "proteomics"


class CobraCompartments:
    """Convert from compartment "long-hand" to "short-hand".

    Shorthand from: https://cobrapy.readthedocs.io/en/latest/_modules/cobra/medium/annotations.html

    "Extracellular" -> "e"
    "golgi" -> "g"
    """

    SHORTHAND: ClassVar[dict[str, list[str]]] = {
        "ce": ["cell envelope"],
        "c": [
            "cytoplasm",
            "cytosol",
            "default",
            "in",
            "intra cellular",
            "intracellular",
            "intracellular region",
            "intracellular space",
        ],
        "er": ["endoplasmic reticulum"],
        "erm": ["endoplasmic reticulum membrane"],
        "e": [
            "extracellular",
            "extraorganism",
            "out",
            "extracellular space",
            "extra organism",
            "extra cellular",
            "extra-organism",
            "external",
            "external medium",
        ],
        "f": ["flagellum", "bacterial-type flagellum"],
        "g": ["golgi", "golgi apparatus"],
        "gm": ["golgi membrane"],
        "h": ["chloroplast"],
        "l": ["lysosome"],
        "im": ["mitochondrial intermembrane space"],
        "mm": ["mitochondrial membrane"],
        "m": ["mitochondrion", "mitochondria"],
        "n": ["nucleus"],
        "p": ["periplasm", "periplasmic space"],
        "x": ["peroxisome", "glyoxysome"],
        "u": ["thylakoid"],
        "vm": ["vacuolar membrane"],
        "v": ["vacuole"],
        "w": ["cell wall"],
        "s": ["eyespot", "eyespot apparatus", "stigma"],
    }

    _REVERSE_LOOKUP: ClassVar[dict[str, list[str]]] = {
        value.lower(): key for key, values in SHORTHAND.items() for value in values
    }

    @classmethod
    def get_shorthand(cls, longhand: str) -> str | None:
        """Get the short-hand compartment name from the long-hand name.

        Args:
            longhand: The long-hand compartment name (e.g., 'cytoplasm', 'extracellular').

        Returns:
            The short-hand compartment name if found, None otherwise.
        """
        return cls._REVERSE_LOOKUP.get(longhand.lower(), None)

    @classmethod
    def get_longhand(cls, shorthand: str) -> str | None:
        """Get the long-hand compartment name from the short-hand name.

        Args:
            shorthand: The short-hand compartment name (e.g., 'c', 'e', 'm').

        Returns:
            The long-hand compartment name if found, None otherwise.
        """
        longhand = cls.SHORTHAND.get(shorthand.lower(), None)
        return longhand[0] if longhand else None


class _BoundaryReactions(NamedTuple):
    """Boundary reactions to be used in the context specific model."""

    reactions: list[str]
    lower_bounds: list[float]
    upper_bounds: list[float]


@dataclass
class _BatchEntry:
    batch_num: int
    sample_names: list[str]
    _num_samples: int = field(init=False)

    def __post_init__(self):
        self._num_samples = len(self.sample_names)

    @property
    def num_samples(self):
        return self._num_samples


@dataclass
class _CombineOmicsInput:
    z_score_matrix: pd.DataFrame
    type: SourceTypes
    weight: int


class _BaseDataType:
    """Base class for common data types."""

    def __getitem__(self, value: str):
        """Access matrices using square bracket notation (e.g., `input_matrices['total_rna']`).

        Args:
            value: The name of the matrix to get ('trna', 'mrna', 'scrna', 'proteomics')

        Returns:
            The DataFrame if it exists, None otherwise.
        """
        self._validate_attribute(value)
        return getattr(self, value)

    def __setitem__(self, key, value):
        """Set matrices using square bracket notation (e.g., `input_matrices['total_rna'] = new_df`).

        :param key: The key to set
        :param value: The new value
        """
        self._validate_attribute(key)
        setattr(self, key, value)

    def _validate_attribute(self, key):
        if key not in {i.value for i in SourceTypes._member_map_.values()}:
            raise ValueError(f"{key} is not a valid attribute of {SourceTypes.__name__}; got '{key}'")

    def __iter__(self) -> Iterator[tuple[SourceTypes, pd.DataFrame | None]]:
        """Iterate over matrix fields and their names.

        Yields:
            A tuple containing (matrix_name, matrix_dataframe).

        """
        for field_ in fields(self):
            yield SourceTypes(field_.name), getattr(self, field_.name)


@dataclass
class _BatchNames(_BaseDataType):
    trna: list[_BatchEntry] = field(default_factory=list)
    mrna: list[_BatchEntry] = field(default_factory=list)
    scrna: list[_BatchEntry] = field(default_factory=list)
    proteomics: list[_BatchEntry] = field(default_factory=list)


@dataclass
class _InputMatrices(_BaseDataType):
    trna: pd.DataFrame | None = None
    mrna: pd.DataFrame | None = None
    scrna: pd.DataFrame | None = None
    proteomics: pd.DataFrame | None = None


@dataclass
class _OutputCombinedSourceFilepath(_BaseDataType):
    trna: Path | None
    mrna: Path | None
    scrna: Path | None
    proteomics: Path | None


@dataclass
class _SourceWeights(_BaseDataType):
    trna: int
    mrna: int
    scrna: int
    proteomics: int


@dataclass
class ModelBuildSettings:
    troppo_epsilon: float = 1e-4
    min_reaction_flux: float = 1e-7
    solver_timeout: int = 1800  # time in seconds, defaults to 30 minutes

    """
    Verbosity
        Type: int
        Default value: 0
        Range: [0, 3]
        From: https://docs.gurobi.com/projects/optimizer/en/current/reference/parameters.html#csclientlog
    """
    solver_verbosity: int = 0

    """
    Feasibility
        Type: double
        Default value: 1e-6
        Range: [1e-9, 1e-2]
        From: https://docs.gurobi.com/projects/optimizer/en/current/reference/parameters.html#feasibilitytol
    """
    solver_feasibility: float = 1e-6

    """
    OptimalityTol
        Type: double
        Default value: 1e-6
        Range: [1e-9, 1e-2]
        From: https://docs.gurobi.com/projects/optimizer/en/current/reference/parameters.html#optimalitytol
    """
    solver_optimality: float = 1e-6

    """
    IntFeasTol
        Type: double
        Default value: 1e-5
        Range: [1e-9, 1e-1]
        From: https://docs.gurobi.com/projects/optimizer/en/current/reference/parameters.html#intfeastol
    """
    solver_integrality: float = 1e-6

    """
    MIPGap
        Relative MIP optimality gap
        Default value: 1e-4
        Range: [0, inf)
        From: https://docs.gurobi.com/projects/optimizer/en/current/reference/parameters.html#mipgap
    """
    gurobi_mipgap: float = 1e-4

    """
    IntegralityFocus
        Type: int
        Default value: 0
        Range: [0, 1]
        From: https://docs.gurobi.com/projects/optimizer/en/current/reference/parameters.html#integralityfocus
    """
    gurobi_integrality_focus: int = 0

    """
    NumericFocus
        Type: int
        Default value: 0
        Range: [0, 3]
        From: https://docs.gurobi.com/projects/optimizer/en/current/reference/parameters.html#numericfocus
    """
    gurobi_numeric_focus: int = 2

    def __post_init__(self):  # noqa: C901
        """Validate provided arguments.

        :raises: ValueError if any check fails.
        """
        if self.troppo_epsilon < 0:
            raise ValueError("ModelBuildSettings: `troppo_epsilon` must be a non-negative float")
        if self.min_reaction_flux < 0:
            raise ValueError("ModelBuildSettings: `min_reaction_flux` must be a non-negative float")
        if self.solver_verbosity not in {0, 1, 2, 3} or not isinstance(self.solver_verbosity, int):
            raise ValueError("ModelBuildSettings: `solver_verbosity` must be an integer in the range [0, 3]")
        if self.solver_timeout < 0 or not isinstance(self.solver_timeout, int):
            raise ValueError("ModelBuildSettings: `solver_timeout` must be a non-negative integer")
        if not (1e-9 < self.solver_feasibility < 1e-2):
            raise ValueError("ModelBuildSettings: `solver_feasibility` must be a float in the range [1e-9, 1e-2]")
        if not (1e-9 < self.solver_optimality < 1e-2):
            raise ValueError("ModelBuildSettings: `solver_optimality` must be a float in the range [1e-9, 1e-2]")
        if not (1e-9 < self.solver_integrality < 1e-1):
            raise ValueError("ModelBuildSettings: `solver_integrality` must be a float in the range [1e-9, 1e-1]")
        if self.gurobi_mipgap < 0:
            raise ValueError("ModelBuildSettings: `gurobi_mipgap` must be a float in the range [1e-4, inf.)")
        if self.gurobi_integrality_focus not in {0, 1} or not isinstance(self.gurobi_integrality_focus, int):
            raise ValueError("ModelBuildSettings: `gurobi_integrality_focus` must be an integer in the range [0, 1]")
        if self.gurobi_numeric_focus not in {0, 1, 2, 3} or not isinstance(self.gurobi_numeric_focus, int):
            raise ValueError("ModelBuildSettings: `gurobi_numeric_focus` must be an integer in the range [0, 3]")
        if self.troppo_epsilon < self.min_reaction_flux * 1000:
            logger.warning(
                "ModelBuildSettings: `troppo_epsilon` and `min_reaction_flux` have similar values. "
                "This can cause inconsistent model builds. Consider ~1000x fold change between the two. "
            )
