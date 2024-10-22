from enum import Enum

from como.como_utilities import stringlist_to_list
from como.project import Config

__all__ = ["stringlist_to_list", "Config", "RNASeqPreparationMethod"]


class RNASeqPreparationMethod(Enum):
    TOTAL = "total"
    MRNA = "mrna"
    SCRNA = "scrna"
