from __future__ import annotations

from enum import Enum
from pathlib import Path
from typing import Literal


class RNASeqPreparationMethod(Enum):
    # class RNAPrepMethod:
    TOTAL = "total"
    MRNA = "mrna"
    SCRNA = "scrna"

    @staticmethod
    def from_string(value: str) -> RNASeqPreparationMethod:
        """Build a preparation method object from a string."""
        match_value = "".join(c for c in value if c.isascii()).lower()

        match match_value:
            case "total" | "trna":
                return RNASeqPreparationMethod.TOTAL
            case "mrna":
                return RNASeqPreparationMethod.MRNA
            case "scrna":
                return RNASeqPreparationMethod.SCRNA
            case _:
                possible_values = [t.value for t in RNASeqPreparationMethod]
                raise ValueError(f"Filtering technique must be one of {possible_values}; got: {value}")


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


type_path = str | Path
type_rna = Literal["total", "mrna"]
