"""
Setup for this file is located in Obsidian (local to my laptop)

Additional file documentation:
Enum information: https://stackoverflow.com/a/1695250
"""
import csv
import os
from enum import Enum
import json
import requests
import sys
import urllib.parse as urlparse
from pprint import pprint
import project


class _InputParameters:
    """
    These are the input parameters required for each of the API endpoints from ProteomeDB
    """

    class ProteinPeptideResult:
        """
        These are the required parameters for obtaining protein expression from ProteomeDB

        PROTEINFILTER: Uniprot Identifier of the protein of interest (String - 64), i.e. P00533
        """

        def __init__(self, protein_filter: str):
            self._protein_filter: str = protein_filter

            self._validate_arguments()

        def _validate_arguments(self) -> None:
            valid_arguments: bool = True
            if len(self._protein_filter) > 64:
                print("The length of the protein ID must be less than 64 characters")
                print(f"You have input a total length of {len(self._protein_filter)}")
                valid_arguments = False

            if not valid_arguments:
                sys.exit(2)

        @property
        def INPUT_PARAMETERS(self) -> str:
            return f"PROTEINFILTER={self._protein_filter}"

    class ProteinsPerExperiment:
        """
        These are the required parameters for obtaining protein information on an experiment-level from ProteomeDB

        EXPERIMENTFILTER: Experiment ID (Integer), the last digits of an experiment URL, i.e. /#projects/4476/132232242 -> 132232242
        """

        def __init__(self, experiment_filter: int):
            self._experiment_filter: int = experiment_filter

        @property
        def INPUT_PARAMETERS(self) -> str:
            return f"EXPERIMENTFILTER={self._experiment_filter}"

    class ProteinExpression:
        """
        These are the required parameters for obtaining protein expression data

        PROTEINFILTER: Uniprot Identifier of the protein of interest (String - 64), i.e. P00533
        MS_LEVEL: MS level used for quantification (Integer), i.e. 1 or 2
        TISSUE_CATEGORY_SELECTION: Category of tissues: tissue, fluid, cell line (String - 2048), i.e. "tissue" or "cell line;fluid"
        TISSUE_ID_SELECTION: List of specific biological sources see LINK TO TISSUE LIST (String - 20)
        SCOPE_SELECTION: Scope of the experiment: 1 - Proteome, 2 - Affinity (Integer)
        CALCULATION_METHOD: Calculation method used to estimate proteine abundance: 0 - iBAQ, 1 - TOP3 (Integer)
        GROUP_BY_TISSUE: Flag indicating whether multiple samples should grouped by biological origin: 1 - yes, 0 - no (Integer)
        EXP_ID: Experiment ID (Integer)
        """

        def __init__(
            self,
            protein_filter: str,
            ms_level: int,
            tissue_category_selection: str,
            tissue_id_selection: str,
            scope_selection: int,
            calculation_method: int,
            group_by_tissue: int,
            experiment_id: int,
        ):
            self._protein_filter: str = protein_filter
            self._ms_level: int = ms_level
            self._tissue_category_selection: str = tissue_category_selection
            self._tissue_id_selection: str = tissue_id_selection
            self._scope_selection: int = scope_selection
            self._calculation_method: int = calculation_method
            self._group_by_tissue: int = group_by_tissue
            self._experiment_id: int = experiment_id

            self._validate_arguments()

        def _validate_arguments(self) -> None:
            valid_arguments: bool = True
            if len(self._protein_filter) > 64:
                print("protein_filter must not be greater than 64 characters")
                print(f"{len(self._protein_filter)} characters were entered")
                valid_arguments = False
            elif self._ms_level not in [1, 2]:
                print("ms_level must be 1 or 2")
                print(f"{self._ms_level} was entered.")
                valid_arguments = False
            elif len(self._tissue_category_selection) > 2048:
                print(
                    "tissue_category_selection must not be greater than 2048 characters"
                )
                print(f"{len(self._tissue_category_selection)} characters were entered")
                valid_arguments = False
            elif len(self._tissue_id_selection) > 20:
                print("tissue_id_selection must not be greater than 20 characters")
                print(f"{len(self._tissue_id_selection)} characters were entered")
                valid_arguments = False
            elif self._scope_selection not in [1, 2]:
                print("scope_selection must be 1 or 2")
                print(f"{self._scope_selection} was entered")
                valid_arguments = False
            elif self._calculation_method not in [0, 1]:
                print("calculation_method must be 0 or 1")
                print(f"{self._calculation_method} was entered")
                valid_arguments = False
            elif self._group_by_tissue not in [0, 1]:
                print("group_by_tissue must be a 0 or 1")
                print(f"{self._group_by_tissue} was entered")
                valid_arguments = False

            if not valid_arguments:
                sys.exit(1)

        @property
        def INPUT_PARAMETERS(self) -> str:

            return (
                f"PROTEINFILTER={self._protein_filter},"
                f"MS_LEVEL={self._ms_level},"
                f"TISSUE_CATEGORY_SELECTION={self._tissue_category_selection},"
                f"TISSUE_ID_SELECTION={self._tissue_id_selection},"
                f"SCOPE_SELECTION={self._scope_selection},"
                f"CALCULATION_METHOD={self._calculation_method},"
                f"GROUP_BY_TISSUE={self._group_by_tissue},"
                f"EXP_ID={self._experiment_id}"
            )

    class ProteinProteotypicityAvailableQuant:
        """
        These are the required parameters for obtaining quantification types for a protein

        PROTEINFILTER: Uniprot Identifier of the protein of interest (String - 64), i.e. P00533
        """

        # PROTEINFILTER = "PROTEINFILTER"
        def __init__(self, protein_filter: str):
            self._protein_filter: str = protein_filter

        def _validate_arguments(self) -> None:
            valid_arguments = True

            if len(self._protein_filter) > 64:
                print("protein_filter must not be longer than 64 characters")
                print(f"{len(self._protein_filter)} characters were entered")
                valid_arguments = False

            if not valid_arguments:
                sys.exit(2)

        @property
        def INPUT_PARAMETERS(self) -> str:
            return f"PROTEINFILTER={self._protein_filter}"

    class ProteinProteotypicity:
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
            self._protein_filter: str = protein_filter
            self._label_type: str = label_type

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

        @property
        def INPUT_PARAMETERS(self) -> str:
            return (
                f"PROTEINFILTER={self._protein_filter},"
                f"LABEL_TYPES={self._label_type}"
            )

    class ProteinsPerTissue:
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
            self._tissue_id: str = tissue_id
            self._calculation_method: int = calculation_method
            self._taxcode: int = taxcode
            self._swissprot_only: int = swissprot_only
            self._no_isoform: int = no_isoform

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

        @property
        def INPUT_PARAMETERS(self) -> str:
            return (
                f"TISSUE_ID={self._tissue_id},"
                f"CALCULATION_METHOD={self._calculation_method},"
                f"SWISSPROT_ONLY={self._swissprot_only},"
                f"NO_ISOFORM={self._no_isoform},"
                f"TAXCODE={self._taxcode}"
            )

    class TissueList:
        """
        No input parameters are required for this API endpoint. It is defined to make sure it is not missed
        """

        pass

    class PeptideSearch:
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
            self._peptide_sequence: str = peptide_sequence
            self._cutoff: float = cutoff
            self._taxcode: int = taxcode

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

        @property
        def INPUT_PARAMETERS(self) -> str:
            return (
                f"PEPTIDE_SEQUENCE={self._peptide_sequence},"
                f"Q_VALUE_CUTOFF={self._cutoff},"
                f"TAXCODE={self._taxcode}"
            )

    class PeptidesPerProtein:
        """
        This is a list of input parameters required to get a comprehensive list of peptides and their modifications of a protein

        PROTEINFILTER: Uniprot identifier (String - 64), i.e. Q92769
        """

        def __init__(self, protein_filter: str):
            self._protein_filter = protein_filter

            self._validate_arguments()

        def _validate_arguments(self) -> None:
            valid_arguments: bool = True
            if len(self._protein_filter) > 64:
                print("The length of the protein ID must be less than 64 characters")
                print(f"You have input a total length of {len(self._protein_filter)}")
                valid_arguments = False

            if not valid_arguments:
                sys.exit(2)

        @property
        def INPUT_PARAMETERS(self) -> str:
            return f"PROTEINFILTER={self._protein_filter}"


class DownloadFormat(Enum):
    """
    This is an Enum of the possible return values from ProteomicsDB API
    """

    JSON = "json"
    XML = "xml"


class LabelTypes(Enum):
    SILAC = "SILAC"
    DIMETHYL = "Dimethyl"
    ITRAQ4 = "iTRAQ4"
    ITRAQ8 = "iTRAQ8"


class URLBuilder:
    def __init__(
        self,
        save_format: DownloadFormat,
        file_save_location: str,
        force_download: bool,
    ):
        """
        import urllib.parse
        url = 'https://example.com/somepage/?'
        params = {'var1': 'some data', 'var2': 1337}
        print(url + urllib.parse.urlencode(params))

        """
        self._base_url: str = "https://www.proteomicsdb.org/proteomicsdb/logic/api/"
        self._force_download: bool = force_download
        self.download_format: str = save_format.name.lower()
        self.file_save_location: str = file_save_location

    def _download_file(self, url: str) -> requests.Response | None:
        # Only download the file content if we have a place to save it, OR we are going to force download
        if self.file_save_location == "" and self._force_download == False:
            print("No file save path provided, not downloading!")
            request = None
        else:
            request = requests.get(url)

        return request

    def _get_file_size(self, url) -> str:
        """
        This function will get the file size of the file to be downloaded
        It will return a string containing the size of the file with an appropriate suffix (KB, MB, etc.)

        :return str: The size of the file
        """
        suffix_dict: dict[str, int] = {
            "B": 1,
            "KB": 1024,
            "MB": 1024**2,  # 1024^2
            "GB": 1024**3,  # 1024^3
            "TB": 1024**4,  # 1024^4
        }

        headers = requests.get(
            url, stream=True, headers={"Accept-Encoding": "identity"}
        ).headers

        # If Content-Length is available to read
        if headers["Content-Length"]:
            file_size = float(headers["Content-Length"])
            for suffix in suffix_dict.keys():
                readable_file_size = file_size / suffix_dict[suffix]

                if readable_file_size < 1024:
                    file_size = readable_file_size
                    file_size = round(file_size, 2)
                    file_size = f"{file_size} {suffix}"
                    break

        # if Content-Length is unavailable
        else:
            file_size = None

        return file_size

    def _write_to_file(self, request: requests.Response) -> None:
        """
        This function will write content from a requests.get() to a file

        It will only write the data if the request was OK and a save location is available
        :param
        """
        if self.file_save_location != "" and request.ok:
            if self._validate_file_extension(
                self.file_save_location, self.download_format
            ):
                with open(self.file_save_location, "wb") as o_stream:
                    o_stream.write(request.content)

    def _validate_file_extension(self, file_path: str, expected_extension: str) -> bool:
        """
        This function will validate the expected file extension matches the requsted return format
        :param file_path: The file path. This can be the full path, or the filename only
        :param expected_extension: The anticipated file extension of the file_path
        """
        correct_extension: bool
        file_extension: str = file_path.split(".")[-1]
        if file_extension == expected_extension:
            correct_extension = True
        else:
            correct_extension = False
        return correct_extension


class GetProteinPeptideResult(URLBuilder):
    """
    This class is responsible for obtaining an aggregated view of protein data
    For a given protein identifier, this API returns the best scored peptide per experiment

    For more information about this api, find the section titled "Get all peptide identifications per given protein"
        at https://www.proteomicsdb.org/#api

    This service provides a highly aggregated view on observed peptides.
    Only the best (in terms of score) peptide per experiment is reported, independent of it's modifications.
    For a more detailed view please check API#9.
    """

    def __init__(
        self,
        protein_identifier: str,
        save_format: DownloadFormat,
        file_save_location: str = "",
        force_download: bool = False,
    ):
        """
        :param protein_identifier: The unique protein identifier from proteomicsdb
        :param save_format: The format to download data in; possible options are JSON or XML
        :param file_save_location: The location to save the file. If no path is provided, the file will not be downladed
        """
        super().__init__(
            save_format=save_format,
            file_save_location=file_save_location,
            force_download=force_download,
        )

        self._parameter = _InputParameters.ProteinPeptideResult(
            protein_identifier
        ).INPUT_PARAMETERS

        self._selection: str = (
            "ENTRY_NAME,PROTEIN_NAME,UNIQUE_IDENTIFIER,TAXCODE,"
            "CHROMOSOME_NAME,GENE_NAME,STRAND,PEPTIDE_SEQUENCE,"
            "PEPTIDE_MASS,START_POSITION,END_POSITION,PEPTIDE_SCORE,"
            "PROTEIN_Q_VALUE,PEPTIDE_Q_VALUE,SEARCH_ENGINE,ISUNIQUE,"
            "ISUNIQUE_PROTEIN,PROJECT_NAME,PROJECT_DESCRIPTION,"
            "EXPERIMENT_NAME,EXPERIMENT_ID,EXPERIMENT_DESCRIPTION,PUBMEDID"
        )
        self._url: str = urlparse.urljoin(
            self._base_url,
            f"proteinpeptideresult.xsodata/InputParams({self._parameter})"
            f"/Results?$select={self._selection}&$format={self.download_format}",
        )
        self._file_size = self._get_file_size(self._url)
        self._request = self._download_file(self._url)

    @property
    def file_size(self) -> str:
        return self._file_size


class GetProteinsPerExperiment(URLBuilder):
    def __init__(
        self,
        experiment_id: int,
        save_format: DownloadFormat,
        file_save_location: str = "",
        force_download: bool = False,
    ):
        """
        This class is responsible for collecting data from the ProteinsPerExperiment ProteomicsDB endpoint
        """
        super().__init__(
            save_format=save_format,
            file_save_location=file_save_location,
            force_download=force_download,
        )
        self._request: requests.Response

        self._parameters = _InputParameters.ProteinsPerExperiment(
            experiment_id
        ).INPUT_PARAMETERS
        self._selection: str = (  # The data to retrieve from the API
            "CHROMOSOME_NAME,COUNT_SPECTRA,DECOY,"
            "END_POSITION,ENTRY_NAME,EVIDENCE,"
            "EXPERIMENT_DESCRIPTION,EXPERIMENT_ID,"
            "EXPERIMENT_NAME,GENE_END,GENE_NAME,"
            "GENE_START,ISUNIQUE,ISUNIQUE_PROTEIN,"
            "OBSERVED_SEQUENCE_COVERAGE_ALL,"
            "PEPTIDE_MASS,PEPTIDE_Q_VALUE,PEPTIDE_SCORE,"
            "PEPTIDE_SEQUENCE,PROJECT_DESCRIPTION,PROJECT_NAME,"
            "PROTEIN_ID,PROTEIN_MASS,PROTEIN_NAME,"
            "PROTEIN_Q_VALUE,PUBMEDID,SEARCH_ENGINE,"
            "SEQUENCE_LENGTH,SPECIES,START_POSITION,"
            "STRAND,TAXCODE,UNIQUE_IDENTIFIER"
        )

        self._url: str = urlparse.urljoin(
            self._base_url,
            f"proteinsperexperiment.xsodata/InputParams({self._parameters})/Results?$select={self._selection}&$format={self.download_format}",
        )

        # Get the file size. This is a fast operation, and may be useful for large downloads
        self._request = self._download_file(self._url)
        self._write_to_file(self._request)

    @property
    def request_content_as_json(self) -> list[dict]:
        print(type(self._request))
        # format = json.loads(content)
        # format = format["d"]["results"]

        return [{}]


class GetProteinExpression(URLBuilder):
    def __init__(
        self,
        protein_filter: str,
        ms_level: int,
        tissue_category_selection: str,
        tissue_id_selection: str,
        scope_selection: int,
        calculation_method: int,
        group_by_tissue: int,
        experiment_id: int,
        download_format: DownloadFormat,
        file_save_location: str = "",
        force_download: bool = False,
    ):
        """
        This function is responsible for colleting protein expression through the proteinexpression ProteomicsDB API
        """
        super().__init__(
            save_format=download_format,
            file_save_location=file_save_location,
            force_download=force_download,
        )
        self._request: requests.Response

        self._parameters = _InputParameters.ProteinExpression(
            protein_filter=protein_filter,
            ms_level=ms_level,
            tissue_category_selection=tissue_category_selection,
            tissue_id_selection=tissue_id_selection,
            scope_selection=scope_selection,
            calculation_method=calculation_method,
            group_by_tissue=group_by_tissue,
            experiment_id=experiment_id,
        ).INPUT_PARAMETERS

        self._selection: str = (  # The data to retrieve from the API
            "UNIQUE_IDENTIFIER,TISSUE_ID,TISSUE_NAME,"
            "TISSUE_SAP_SYNONYM,SAMPLE_ID,SAMPLE_NAME,"
            "AFFINITY_PURIFICATION,EXPERIMENT_ID,EXPERIMENT_NAME,"
            "EXPERIMENT_SCOPE,EXPERIMENT_SCOPE_NAME,PROJECT_ID,"
            "PROJECT_NAME,UNNORMALIZED_INTENSITY,NORMALIZED_INTENSITY,"
            "MIN_NORMALIZED_INTENSITY,MAX_NORMALIZED_INTENSITY,"
            "SAMPLES,PROJECT_STATUS"
        )
        # /proteomicsdb/logic/api/proteinexpression.xsodata/InputParams(
        self._url: str = urlparse.urljoin(
            self._base_url,
            f"proteinexpression.xsodata/InputParams({self._parameters})/Results?$select={self._selection}&$format={self.download_format}",
        )
        self._file_size = self._get_file_size(self._url)

        # Get the file size
        print(f"ProteinExpression file size of {self._file_size}")

    @property
    def file_size(self) -> str:
        return self._file_size


def read_experiment_ids(file_path: str) -> list[int]:
    """
    This function is responsible for reading experiment ID values from a file

    :param file_path: The full file path to the experiment IDs (/my/full/path/proteomicdb_experiment_ids.csv)
    """
    experiment_ids: list[int] = []

    with open(file_path, "r") as i_stream:
        reader = csv.DictReader(i_stream, delimiter=",")
        for row in reader:
            try:
                experiment_ids.append(int(row["experiment_id"]))
            except ValueError:
                print(
                    f"Unable to convert experiment ID value {row['experiment_id']} to an integer!"
                )

    return experiment_ids


def download_controller(experiment_id: int):
    """
    This function is responsible for controlling the download flow of ProteinsPerExperiment -> ProteinExpression
    """
    file_save_dir = "/Users/joshl/Library/Application Support/JetBrains/PyCharm2022.1/scratches/proteomics"
    experiment_file_save = (
        f"{file_save_dir}/{experiment_id}_proteins_per_experiment.json"
    )

    proteins_per_experiment = GetProteinsPerExperiment(
        experiment_id=experiment_id,
        save_format=DownloadFormat.JSON,
        force_download=True
        # file_save_location=experiment_file_save,
    )


def main():
    """
    Empty values can be input for items that are not required to filter by
    The following API command is able to get ALL MS_LEVEL information pertaining to Experiment ID 13207 and Protein Identifier P00533
    https://www.proteomicsdb.org/proteomicsdb/logic/api/proteinexpression.xsodata/InputParams(PROTEINFILTER=''P00533'',MS_LEVEL=0,TISSUE_ID_SELECTION='''',TISSUE_CATEGORY_SELECTION='''',SCOPE_SELECTION=1,GROUP_BY_TISSUE=0,CALCULATION_METHOD=0,EXP_ID=13207)/Results?$select=UNIQUE_IDENTIFIER,TISSUE_ID,TISSUE_NAME,TISSUE_SAP_SYNONYM,SAMPLE_ID,SAMPLE_NAME,AFFINITY_PURIFICATION,EXPERIMENT_ID,EXPERIMENT_NAME,EXPERIMENT_SCOPE,EXPERIMENT_SCOPE_NAME,PROJECT_ID,PROJECT_NAME,PROJECT_STATUS,UNNORMALIZED_INTENSITY,NORMALIZED_INTENSITY,MIN_NORMALIZED_INTENSITY,MAX_NORMALIZED_INTENSITY,SAMPLES&$format=json
    """
    config = project.configs
    experiment_file = os.path.join(config.configdir, "proteomicdb_experiment_ids.csv")
    experiment_ids = read_experiment_ids(experiment_file)

    for i, exp_id in enumerate(experiment_ids):
        download_controller(exp_id)


if __name__ == "__main__":
    main()
