from __future__ import annotations

from collections.abc import Iterator
from dataclasses import dataclass, field, fields
from enum import Enum
from pathlib import Path
from typing import ClassVar, NamedTuple

import cobra
import pandas as pd
from loguru import logger

PATH_TYPE = str | Path
LOG_FORMAT = (
    "<green>{time:YYYY-MM-DD HH:mm:ss}</> | "
    "<level>{level:<8}</> | "
    "<cyan>{name}</>:<cyan>{line}</> - <level>{message}</>"
)


class AdjustmentMethod(Enum):
    """Adjustment method for expression requirement based on differences in number of provided data source types."""

    PROGRESSIVE = "progressive"
    REGRESSIVE = "regressive"
    FLAT = "flat"
    CUSTOM = "custom"


class Algorithm(Enum):
    GIMME = "GIMME"
    FASTCORE = "FASTCORE"
    IMAT = "IMAT"
    TINIT = "TINIT"


class FilteringTechnique(Enum):
    """RNA sequencing filtering capabilities."""

    CPM = "cpm"
    ZFPKM = "zfpkm"
    TPM = "tpm"
    UMI = "umi"


class GeneIdentifier(Enum):
    ENSEMBL_GENE_ID = "ensembl_gene_id"
    ENTREZ_GENE_ID = "entrez_gene_id"
    GENE_SYMBOL = "gene_symbol"


class LogLevel(Enum):
    TRACE = 5
    DEBUG = 10
    INFO = 20
    SUCCESS = 25
    WARNING = 30
    ERROR = 40
    CRITICAL = 50


class PeakIdentificationParameters(NamedTuple):
    height: float
    distance: float


class RNAType(Enum):
    TRNA = "total"
    MRNA = "mrna"
    SCRNA = "scrna"


class Solver(Enum):
    """Solver used to seed context specific model."""

    GLPK = "GLPK"
    GUROBI = "GUROBI"
    SCIPY = "SCIPY"
    GLPK_EXACT = "GLPK_EXACT"


class SourceTypes(Enum):
    TRNA = "trna"
    MRNA = "mrna"
    SCRNA = "scrna"
    PROTEOMICS = "proteomics"


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
        """Get the short-hand compartment name from the long-hand name."""
        return cls._REVERSE_LOOKUP.get(longhand.lower(), None)

    @classmethod
    def get_longhand(cls, shorthand: str) -> str | None:
        """Get the long-hand compartment name from the short-hand name."""
        longhand = cls.SHORTHAND.get(shorthand.lower(), None)
        return longhand[0] if longhand else None


class _BuildResults(NamedTuple):
    """Results of building a context specific model."""

    model: cobra.Model
    expression_index_list: list[int]
    infeasible_reactions: pd.DataFrame


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

        :param value: The name of the matrix to get ('trna', 'mrna', 'scrna', 'proteomics')
        :returns: The DataFrame if it exists, None otherwise.
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
            # Unable to use como.utils._log_and_raise_error becuase it results in a circular import
            message = f"{key} is not a valid attribute of {SourceTypes.__name__}; got '{key}'"
            logger.warning(message)
            raise ValueError(message)

    def __iter__(self) -> Iterator[tuple[SourceTypes, pd.DataFrame | None]]:
        """Iterate over matrix fields and their names.

        Yields:
            A tuple containing (matrix_name, matrix_dataframe).

        """
        for field_ in fields(self):
            yield SourceTypes(field_.name), getattr(self, field_.name)


@dataclass
class _BatchNames(_BaseDataType):
    trna: list[_BatchEntry]
    mrna: list[_BatchEntry]
    scrna: list[_BatchEntry]
    proteomics: list[_BatchEntry]


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
