from __future__ import annotations

from abc import ABC
from enum import Enum
from pathlib import Path
from typing import Literal

from pydantic import BaseModel, ConfigDict
from pydantic_settings import BaseSettings, SettingsConfigDict


class _BaseModel(BaseModel, ABC):
    model_config = ConfigDict()


class _BaseArguments(BaseSettings, ABC):
    model_config = SettingsConfigDict(
        cli_parse_args=True,
        cli_kebab_case=True,
        nested_model_default_partial_update=True,
        case_sensitive=True,
        cli_avoid_json=True,
        cli_enforce_required=True,
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


type_path = str | Path
type_rna = Literal["total", "polya"]
