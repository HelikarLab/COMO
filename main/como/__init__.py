from enum import Enum

from como.como_utilities import stringlist_to_list
from como.project import Config

__all__ = ["stringlist_to_list", "Config", "RNASeqPreparationMethod"]
__version__ = "1.10.0"


class RNASeqPreparationMethod(Enum):
    TOTAL = "total"
    MRNA = "mrna"
    SCRNA = "scrna"
