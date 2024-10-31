from enum import Enum


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
