from pathlib import Path
import requests
import os
import json
from parameters import *


class URLBuilder:
    def __init__(
        self,
        api_endpoint: str,
        parameters: type(Parameters),
        selection: str,
        download_format: DownloadFormat,
        file_save_location: Path,
        force_download: bool,
        session: requests.Session,
        print_file_size: bool = True,
    ):
        """
        This is the base class for the ProteomicsDB API
        It is responsible for the following high-level actions:
            - Downloading data
            - Writing data to files

        This class should not exit the program on its own. It should simply return success/error values
            Any subclasses created from this class should handle errors on their own

        """
        self._parameters: str = parameters.INPUT_PARAMETERS
        self._force_download: bool = force_download
        self._print_file_size = print_file_size
        self._download_format: str = download_format.name.lower()
        self._file_save_location: Path = file_save_location
        self._api_url: str = f"https://www.proteomicsdb.org/proteomicsdb/logic/api/{api_endpoint}.xsodata/InputParams({self._parameters})/Results?$select={selection}&$format={download_format.value}"
        if session is None:
            self._session = requests.Session()
        else:
            self._session = session

        self._byte_file_size: int = self._get_file_size()
        self._readable_file_size: str = self._format_file_size()

        self._request: requests.Response = self._download_file()

    def _download_file(self) -> requests.Response:
        # Only download the file content if we have a place to save it, OR we are going to force download
        if self._file_save_location == "" and self._force_download == False:
            print("No file save path provided, not downloading!")
            download_request = None
        else:
            file_size = self._get_file_size()
            formatted_size = self._format_file_size()
            if self._print_file_size:
                print(f"Starting download, file size is: {formatted_size}")
            download_request = self._session.get(self._api_url)

        return download_request

    def _get_file_size(self) -> int:
        """
        This function will get the file size of the file to be downloaded
        It will return an integer of the file size in bytes

        :return int: The binary size of the file
        """
        file_size = 0
        headers = self._session.get(
            self._api_url, stream=True, headers={"Accept-Encoding": "identity"}
        ).headers

        if headers["Content-Length"]:
            file_size = int(headers["Content-Length"])

        return file_size

    def _write_to_file(self) -> bool:
        """
        This function will write content from a requests.get() to a file

        It will only write the data if the response was OK and a save location is available
        :param
        """
        file_written: bool = False
        file_extension = str(self._file_save_location).split(".")[-1]
        valid_file_extension = self._validate_file_extension(
            file_extension, self._download_format
        )

        # Validate the response was OK and the correct file extension was requested/downloaded
        if self._request.ok and valid_file_extension:
            json_object = self._request.json()["d"]["results"]

            os.makedirs(os.path.dirname(self._file_save_location), exist_ok=True)
            with open(self._file_save_location, "w") as o_stream:
                json.dump(json_object, o_stream)
                file_written = True

        elif not self._request.ok:
            print(
                f"An error occurred trying to save the proteomics data. "
                f"The API response returned a status code of {self._request.status_code}"
            )
            file_written = False
        elif not valid_file_extension:
            print(
                f"The requested download format does not match the expected file save format."
            )
            print(f"The download was requested in a format of {self._download_format}")
            print(
                f"The file is going to be saved with an extension of '{file_extension}'"
            )
            file_written = False

        return file_written

    def _read_from_file(self) -> list[dict]:
        file_data = open(self._file_save_location, "r").read()
        json_data = json.loads(file_data)
        return json_data

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

    def _format_file_size(self) -> str:
        """
        This function will return a string of the new formatted size and its size description

        For example: A size of 2048 bytes is passed in (original_size = 2014)
        This function will loop through possible sizes until a size under 1024 (file base sizes) is reached
        Then, it will make a new string containing this new file size (2) with the suffix (KB)
        """
        formatted_size: float = 0.0
        suffix: str = "bytes"
        for size in FileSizeFormats:
            formatted_size = self._byte_file_size / size.value
            if formatted_size <= 1024:
                suffix = size.name
                break

        formatted_size = round(formatted_size, 2)
        return_size = f"{formatted_size} {suffix}"
        return return_size

    def __len__(self):
        """
        Return the length of the JSON response
        This function is called when using the "len()" function
        """
        return len(self.json_response)

    @property
    def file_byte_size(self) -> float:
        return self._byte_file_size

    @property
    def file_readable_size(self) -> str:
        return self._format_file_size()

    @property
    def json_response(self) -> list[dict]:
        return self._request.json()["d"]["results"]

    @property
    def response_object(self) -> requests.Response:
        return self._request

    @property
    def save_path(self):
        return self._file_save_location


class GetProteinsPerExperiment(URLBuilder):
    def __init__(
        self,
        experiment_id: int | list[int],
        download_format: DownloadFormat,
        file_save_location: Path = Path(),
        force_download: bool = False,
        session: requests.Session = None,
        **kwargs,
    ):
        """
        This class is responsible for collecting data from the ProteinsPerExperiment ProteomicsDB endpoint
        """
        self._experiment_id: int = experiment_id
        self._download_format: DownloadFormat = download_format
        self._file_save_location: Path = file_save_location
        self._parameters: ProteinsPerExperiment = ProteinsPerExperiment(
            self._experiment_id
        )
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

        super().__init__(
            api_endpoint="proteinsperexperiment",
            parameters=self._parameters,
            selection=self._selection,
            download_format=self._download_format,
            file_save_location=self._file_save_location,
            force_download=force_download,
            session=session,
            **kwargs,
        )

    def save_data(self):
        """
        This is simply a wrapper around the URLBuilder _write_to_file() function.
        It is done this way so that we can simply call save_data() from the class instance, without having to worry about extra arguments
        """
        file_written: bool = super()._write_to_file()
        if not file_written:
            sys.exit(2)

    @property
    def experiment_id(self):
        return self._experiment_id


class GetProteinExpression(URLBuilder):
    def __init__(
        self,
        protein_filter: str,
        experiment_id: int,
        download_format: DownloadFormat,
        calculation_method: ProteinExpression.CalculationMethod = ProteinExpression.CalculationMethod.TOP3,
        scope_selection: ProteinExpression.ScopeSelection = ProteinExpression.ScopeSelection.ANY,
        group_by_tissue: ProteinExpression.GroupByTissue = ProteinExpression.GroupByTissue.FALSE,
        tissue_category_selection: str = "",
        tissue_id_selection: str = "",
        ms_level: ProteinExpression.MSLevel = ProteinExpression.MSLevel.ANY,
        file_save_location: Path = Path(),
        force_download: bool = False,
        session: requests.Session = None,
        **kwargs,
    ):
        """
        This function is responsible for colleting protein expression through the proteinexpression ProteomicsDB API
        """
        self._protein_id = protein_filter
        self._experiment_id = experiment_id
        self._download_format = download_format
        self._file_save_location = file_save_location
        self._parameters = ProteinExpression(
            protein_filter=protein_filter,
            ms_level=ms_level,
            tissue_category_selection=tissue_category_selection,
            tissue_id_selection=tissue_id_selection,
            scope_selection=scope_selection,
            calculation_method=calculation_method,
            group_by_tissue=group_by_tissue,
            experiment_id=experiment_id,
        )

        self._selection: str = (  # The data to retrieve from the API
            "UNIQUE_IDENTIFIER,TISSUE_ID,TISSUE_NAME,"
            "TISSUE_SAP_SYNONYM,SAMPLE_ID,SAMPLE_NAME,"
            "AFFINITY_PURIFICATION,EXPERIMENT_ID,EXPERIMENT_NAME,"
            "EXPERIMENT_SCOPE,EXPERIMENT_SCOPE_NAME,PROJECT_ID,"
            "PROJECT_NAME,UNNORMALIZED_INTENSITY,NORMALIZED_INTENSITY,"
            "MIN_NORMALIZED_INTENSITY,MAX_NORMALIZED_INTENSITY,"
            "SAMPLES,PROJECT_STATUS"
        )
        super().__init__(
            api_endpoint="proteinexpression",
            parameters=self._parameters,
            selection=self._selection,
            download_format=self._download_format,
            file_save_location=self._file_save_location,
            force_download=force_download,
            session=session,
            **kwargs,
        )

    def save_data(self):
        """
        This is a wrapper around URLBuilder's _write_to_file() function
        This is done this way so we can call save_data() from any class instance, and easily change its implementation
        """
        file_written: bool = super()._write_to_file()
        if not file_written:
            sys.exit(2)

    @property
    def protein_id(self):
        return self._protein_id

    @property
    def experiment_id(self):
        return self._experiment_id
