from __future__ import annotations

from enum import Enum
from pathlib import Path
from typing import Literal

PATH_TYPE = str | Path
LOG_FORMAT = (
    "<green>{time:YYYY-MM-DD HH:mm:ss}</> | "
    "<level>{level:<8}</> | "
    "<cyan>{name}</>:<cyan>{line}</> - <level>{message}</>"
)


class RNAPrepMethod(Enum):
    TOTAL = "total"
    MRNA = "mrna"
    SCRNA = "scrna"

    @staticmethod
    def from_string(value: str) -> RNAPrepMethod:
        """Build a preparation method object from a string."""
        match_value = "".join(c for c in value if c.isascii()).lower()

        match match_value:
            case "total" | "trna":
                return RNAPrepMethod.TOTAL
            case "mrna":
                return RNAPrepMethod.MRNA
            case "scrna":
                return RNAPrepMethod.SCRNA
            case _:
                possible_values = [t.value for t in RNAPrepMethod]
                raise ValueError(f"Filtering technique must be one of {possible_values}; got: {value}")
class Algorithm(Enum):
    GIMME = "GIMME"
    FASTCORE = "FASTCORE"
    IMAT = "IMAT"
    TINIT = "TINIT"


class FilteringTechnique(Enum):
    """RNA sequencing filtering capabilities."""

    cpm = "cpm"
    zfpkm = "zfpkm"
    tpm = "quantile"
    umi = "umi"

    @staticmethod
    def from_string(value: str) -> FilteringTechnique:
        """Create a filtering technique object from a string."""
        match value.lower():
            case "cpm":
                return FilteringTechnique.cpm
            case "zfpkm":
                return FilteringTechnique.zfpkm
            case "quantile":
                return FilteringTechnique.tpm
            case "umi":
                return FilteringTechnique.umi
            case _:
                possible_values = [t.value for t in FilteringTechnique]
                raise ValueError(f"Got a filtering technique of '{value}'; should be one of: {possible_values}")


PATH_TYPE = str | Path
RNA_TYPE = Literal["total", "mrna"]
LOG_LEVEL = Literal["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"]

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

