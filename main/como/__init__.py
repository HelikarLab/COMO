from enum import Enum

import pandas as pd

from como.como_utilities import stringlist_to_list
from como.project import Config

__all__ = ["stringlist_to_list", "Config", "RNASeqPreparationMethod"]
__version__ = "1.10.0"


class RNASeqPreparationMethod(Enum):
    TOTAL = "total"
    MRNA = "mrna"
    SCRNA = "scrna"

    @staticmethod
    def from_string(value: str) -> "RNASeqPreparationMethod":
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


def return_placeholder_data() -> pd.DataFrame:
    return pd.DataFrame(data=0, index=pd.Index(data=[0], name="entrez_gene_id"), columns=["expressed", "top"])
