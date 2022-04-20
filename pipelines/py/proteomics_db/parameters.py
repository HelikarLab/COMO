"""
|---------------------------------------------------------------------------------------------------------------------|
|  READ ME                                                                                                            |
|                                                                                                                     |
|  No "naked" code should exist in this file. It is meant to hold the API implementation ONLY.                        |
|  If you feel you must include code OTHER THAN that required for API implementations, YOU are wrong. Not this file.  |
|---------------------------------------------------------------------------------------------------------------------|


The purpose of this file is to hold the required parameters for ProteomicsDB API
--------------------
Class Implementation
_Parameters: A class devoted to holding required input parameters for each of the APIs available at ProteomicsDB
DownloadFormat: The format to download data in (JSON or XML)
LabelTypes: The label type for ProteinProteotypicity. I haven't been able to figure out how to place these values inside the class, so it exists as a standalone class right now
--------------------
"""

import sys
from enum import Enum


class DownloadFormat(Enum):
    """
    This is an Enum of the possible return values from ProteomicsDB API
    """

    JSON = "json"
    XML = "xml"


class FileSizeFormats(Enum):
    """
    This is a class to convert bytes to a different file size

    Dividing a byte value by any one of these results in that file size being returned
    """

    B = 1  # Use native byte-value
    KB = 1024  # 1024^1, kilobytes
    MB = 1024**2  # 1024^2, megabytes
    GB = 1024**3  # 1024^3, gigabytes
    TB = 1024**4  # 1024^4, terabytes


class LabelTypes(Enum):
    SILAC = "SILAC"
    DIMETHYL = "Dimethyl"
    ITRAQ4 = "iTRAQ4"
    ITRAQ8 = "iTRAQ8"


class Parameters:
    """
    This is the base class for additional API classes to utilize

    These are the required input parameters for each of the API endpoints of ProteomicsDB
    """

    def __init__(self):
        self._parameters: str = ""

    @property
    def INPUT_PARAMETERS(self) -> str:
        return self._parameters


class ProteinPeptideResult(Parameters):
    def __init__(self, protein_filter: str):
        """
        These are the required parameters for obtaining protein expression from ProteomeDB

        PROTEINFILTER: Uniprot Identifier of the protein of interest (String - 64), i.e. P00533
        """
        super().__init__()
        self._protein_filter: str = protein_filter
        self._parameters: str = f"PROTEINFILTER={self._protein_filter}"
        self._validate_arguments()

    def _validate_arguments(self) -> None:
        valid_arguments: bool = True
        if len(self._protein_filter) > 64:
            print("The length of the protein ID must be less than 64 characters")
            print(f"You have input a total length of {len(self._protein_filter)}")
            valid_arguments = False

        if not valid_arguments:
            sys.exit(2)


class ProteinsPerExperiment(Parameters):
    """
    These are the required parameters for obtaining protein information on an experiment-level from ProteomeDB

    EXPERIMENTFILTER: Experiment ID (Integer), the last digits of an experiment URL, i.e. /#projects/4476/132232242 -> 132232242
    """

    def __init__(self, experiment_filter: int):
        super().__init__()
        self._experiment_filter: int = experiment_filter
        self._parameters: str = f"EXPERIMENTFILTER={self._experiment_filter}"


class ProteinExpression(Parameters):
    """
    These are the required parameters for obtaining protein expression data

    PROTEINFILTER: Uniprot Identifier of the protein of interest (String - 64), i.e. P00533
    MSLevel: MS level used for quantification (Integer), i.e. 1 or 2
    TISSUE_CATEGORY_SELECTION: Category of tissues: tissue, fluid, cell line (String - 2048), i.e. "tissue" or "cell line;fluid"
    TISSUE_ID_SELECTION: List of specific biological sources see LINK TO TISSUE LIST (String - 20)
    SCOPE_SELECTION: Scope of the experiment: 1 - Proteome, 2 - Affinity (Integer)
    CALCULATION_METHOD: Calculation method used to estimate proteine abundance: 0 - iBAQ, 1 - TOP3 (Integer)
    GROUP_BY_TISSUE: Flag indicating whether multiple samples should grouped by biological origin: 1 - yes, 0 - no (Integer)
    EXP_ID: Experiment ID (Integer)
    """

    class CalculationMethod(Enum):
        """
        This is an enum of the possible calculation methods for mass spectrometry
        """

        IBAQ = 0
        TOP3 = 1

    class ScopeSelection(Enum):
        PROTEOME = 1
        AFFINITY = 2
        ANY = 0

    class GroupByTissue(Enum):
        TRUE = 1
        FALSE = 0

    class MSLevel(Enum):
        ONE = 1
        TWO = 2
        ANY = 0

    def __init__(
        self,
        protein_filter: str,
        ms_level: MSLevel,
        tissue_category_selection: str,
        tissue_id_selection: str,
        scope_selection: ScopeSelection,
        calculation_method: CalculationMethod,
        group_by_tissue: GroupByTissue,
        experiment_id: int,
    ):
        super().__init__()
        self._group_by_tissue: ProteinExpression.GroupByTissue = group_by_tissue.value
        self._scope_selection: ProteinExpression.ScopeSelection = scope_selection.value
        self._calculation_method: ProteinExpression.CalculationMethod = (
            calculation_method.value
        )

        self._ms_level: ProteinExpression.MSLevel = ms_level.value
        self._experiment_id: int = experiment_id
        self._protein_filter: str = f"'{protein_filter}'"
        self._tissue_id_selection: str = f"'{tissue_id_selection}'"
        self._tissue_category_selection: str = f"'{tissue_category_selection}'"

        self._validate_arguments()

        self._parameters: str = (
            f"PROTEINFILTER={self._protein_filter},"
            f"MS_LEVEL={self._ms_level},"
            f"TISSUE_CATEGORY_SELECTION={self._tissue_category_selection},"
            f"TISSUE_ID_SELECTION={self._tissue_id_selection},"
            f"SCOPE_SELECTION={self._scope_selection},"
            f"CALCULATION_METHOD={self._calculation_method},"
            f"GROUP_BY_TISSUE={self._group_by_tissue},"
            f"EXP_ID={self._experiment_id}"
        )

    def _validate_arguments(self) -> None:
        valid_arguments: bool = True
        if len(self._protein_filter) > 64:
            print("protein_filter must not be greater than 64 characters")
            print(f"{len(self._protein_filter)} characters were entered")
            valid_arguments = False
        elif len(self._tissue_category_selection) > 2048:
            print("tissue_category_selection must not be greater than 2048 characters")
            print(f"{len(self._tissue_category_selection)} characters were entered")
            valid_arguments = False
        elif len(self._tissue_id_selection) > 20:
            print("tissue_id_selection must not be greater than 20 characters")
            print(f"{len(self._tissue_id_selection)} characters were entered")
            valid_arguments = False

        if not valid_arguments:
            sys.exit(1)


class ProteinProteotypicityAvailableQuant(Parameters):
    """
    These are the required parameters for obtaining quantification types for a protein

    PROTEINFILTER: Uniprot Identifier of the protein of interest (String - 64), i.e. P00533
    """

    # PROTEINFILTER = "PROTEINFILTER"
    def __init__(self, protein_filter: str):
        super().__init__()
        self._protein_filter: str = protein_filter
        self._parameters: str = f"PROTEINFILTER={self._protein_filter}"

        self._validate_arguments()

    def _validate_arguments(self) -> None:
        valid_arguments = True

        if len(self._protein_filter) > 64:
            print("protein_filter must not be longer than 64 characters")
            print(f"{len(self._protein_filter)} characters were entered")
            valid_arguments = False

        if not valid_arguments:
            sys.exit(2)


class ProteinProteotypicity(Parameters):
    """
    These are the required parameters for obtaining experimental proteotypicity for peptides of a specific protein

    PROTEINFILTER: Uniprot Identifier of the protein of interest (String - 64), i.e. P00533
    LABEL_TYPES: Semicolon separated label types which should be used to calculte peptide proteotypicities (String - 128), i.e. SILAC;TMT


        The following is a list of quantification IDs (LABEL_TYPES) returned by ProteomicsDB. I have not been able to retrieve results with all of them, so some are not validated.
        SILAC quantitation analysis (VALID)
        Dimethyl                    (VALID)
        iTRAQ4                      (VALID)
        iTRAQ8                      (VALID)
        TMT quantitation analysis   (UNKNOWN)
        Label-free: Intensity       (UNKNOWN)
        TMT10-plex                  (UNKNOWN)
    """

    def __init__(self, protein_filter: str, label_type: str):
        super().__init__()
        self._protein_filter: str = protein_filter
        self._label_type: str = label_type
        self._parameters: str = (
            f"PROTEINFILTER={self._protein_filter}," f"LABEL_TYPES={self._label_type}"
        )

        self._valid_label_types: set[str] = {item.value for item in LabelTypes}

    def _validate_arguments(self) -> None:
        valid_arguments: bool = True
        if len(self._protein_filter) > 64:
            print("The length of the protein ID must be less than 64 characters")
            print(f"A total length of {len(self._protein_filter)} was entered")
            valid_arguments = False
        elif self._label_type not in self._valid_label_types:
            print("An invalid label_type was entered")
            print(f"Input label_type: {self._label_type}")
            print(f"Valid label_type: {self._valid_label_types}")
            valid_arguments = False

        if not valid_arguments:
            sys.exit(2)


class ProteinsPerTissue(Parameters):
    """
    This is a list of required inputs to retrieve proteins found on a per-tissue level

    TISSUE_ID: Brenda tissue ontology ID (String - 20), i.e. BTO:0000975
    CALCULATION_METHOD: 0: iBAQ; 1: Top3 (Integer)
    SWISSPROT_ONLY: 1: only returns Swissprot entries; 0: all other (Integer)
    NO_ISOFORM: 1: no isoforms are shown; 0: shows all (Integer)
    TAXCODE: Taxonomy code (Integer), i.e. 9606
    """

    def __init__(
        self,
        tissue_id: str,
        calculation_method: int,
        swissprot_only: bool,
        no_isoform: bool,
        taxcode: int,
    ):
        super().__init__()
        self._tissue_id: str = tissue_id
        self._calculation_method: int = calculation_method
        self._taxcode: int = taxcode
        self._swissprot_only: int = swissprot_only
        self._no_isoform: int = no_isoform

        self._parameters: str = (
            f"TISSUE_ID={self._tissue_id},"
            f"CALCULATION_METHOD={self._calculation_method},"
            f"SWISSPROT_ONLY={self._swissprot_only},"
            f"NO_ISOFORM={self._no_isoform},"
            f"TAXCODE={self._taxcode}"
        )

        self._validate_arguments()

    def _validate_arguments(self) -> None:
        valid_arguments: bool = True
        if self._calculation_method not in [0, 1]:
            print("calcuation_method must be 0 (iBAQ) or 1 (Top3)")
            print(f"A value of {self._calculation_method} was entered")
            valid_arguments = False
        elif self._swissprot_only not in [0, 1]:
            print("swissprot_only must be 0 or 1")
            print(f"0 = False, allow non-swissprot values")
            print("1 = True, only allow swissprot values")
            print(f"{self._swissprot_only} was entered")
            valid_arguments = False
        elif self._no_isoform not in [0, 1]:
            print("no_isoform must be 0 or 1")
            print(f"0 = False, show all values")
            print("1 = True, only show isoforms")
            print(f"{self._no_isoform} was entered")
            valid_arguments = False

        if not valid_arguments:
            sys.exit(2)


class TissueList(Parameters):
    """
    No input parameters are required for this API endpoint. It is defined to make sure it is not missed
    """

    def __init__(self):
        super().__init__()


class PeptideSearch(Parameters):
    """
    This is a list of required inputs to retrieve comprehensive information about a peptide

    PEPTIDE_SEQUENCE: Peptide sequence (String - 500), i.e. GSHQISLDNPDYQQDFFPK
    Q_VALUE_CUTOFF: q-value cutoff (Double), i.e. 0.01

    """

    def __init__(self, peptide_sequence: str, cutoff: float, taxcode: int):
        """
        :params sequence: This is the peptide sequence to search for. It has a maximum length of 500 amino acids
        :params cutoff: This is the q-value cutoff.
        """
        super().__init__()
        self._peptide_sequence: str = peptide_sequence
        self._cutoff: float = cutoff
        self._taxcode: int = taxcode

        self._parameters: str = (
            f"PEPTIDE_SEQUENCE={self._peptide_sequence},"
            f"Q_VALUE_CUTOFF={self._cutoff},"
            f"TAXCODE={self._taxcode}"
        )

        self._validate_arguments()

    def _validate_arguments(self) -> None:
        valid_arguments: bool = True

        if len(self._peptide_sequence) > 500:
            print("Maximum peptide_sequence is 500 characters")
            print(f"You have input a length of {len(self._peptide_sequence)}")
            valid_arguments = False
        elif len(str(self._taxcode)) != 4:
            print("The taxcode argument must be an integer with 4 digits")
            print("Example: 9606")
            print(f"A value of {self._taxcode} was entered")
            valid_arguments = False

        if not valid_arguments:
            sys.exit(2)


class PeptidesPerProtein(Parameters):
    """
    This is a list of input parameters required to get a comprehensive list of peptides and their modifications of a protein

    PROTEINFILTER: Uniprot identifier (String - 64), i.e. Q92769
    """

    def __init__(self, protein_filter: str):
        super().__init__()
        self._protein_filter = protein_filter

        self._parameters: str = f"PROTEINFILTER={self._protein_filter}"

        self._validate_arguments()

    def _validate_arguments(self) -> None:
        valid_arguments: bool = True
        if len(self._protein_filter) > 64:
            print("The length of the protein ID must be less than 64 characters")
            print(f"You have input a total length of {len(self._protein_filter)}")
            valid_arguments = False

        if not valid_arguments:
            sys.exit(2)
