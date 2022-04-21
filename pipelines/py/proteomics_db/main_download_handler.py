from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import argparse
import csv
import endpoints
import json
import os
import requests
import requests.adapters
import sys
from .. import project


class DownloadHandler:
    def __init__(self, args: argparse.Namespace):
        """
        This function will handle the following functions:
            - Read experiment IDs from the config_file (from 'args' function)
            - Download experiment data using experiment IDs
                - Save the resulting data in a semi-readable format (one CSV file per experiment)
            - Download protein expression using each of the protein IDs and experiment IDs
                - Save the resulting data in a semi-readable format (one CSV file )
        """
        self._config_file: str = args.config_file
        self._file_type: endpoints.DownloadFormat = args.file_type
        self._output_dir: str = args.output_dir
        self._experiment_ids: list[int] = self.read_experiment_ids()

        # Download proteins per experiment
        self._protein_per_experiment_data: list[
            Path
        ] = self.save_proteins_per_experiment()

        self._protein_expression_data: list[Path] = self.save_protein_expression()

    def read_experiment_ids(self) -> list[int]:
        """
        This function is responsible for reading experiment ID values from a file
        """
        experiment_ids: list[int] = []

        with open(self._config_file, "r") as i_stream:
            reader = csv.DictReader(i_stream, delimiter=",")
            for row in reader:
                try:
                    experiment_ids.append(int(row["experiment_id"]))
                except ValueError:
                    print(
                        f"Unable to convert experiment ID value {row['experiment_id']} to an integer! Validate the experiment IDs found in {self._config_file}."
                    )

        return experiment_ids

    def parse_experiment_json(self, experiment_path: Path) -> dict:
        """
        This function will accept a JSON file, and retrieve the following keys:
            - UNIQUE_IDENTIFIER (Corresponds to PROTEINFILTER in ProteomicsDB API)
            - EXPERIMENT_ID
        """
        json_data: dict = {
            "UNIQUE_IDENTIFIER": [],
            "EXPERIMENT_ID": [],
        }
        with open(experiment_path, "r") as i_stream:
            json_object = json.load(i_stream)

        for entry in json_object:
            json_data["UNIQUE_IDENTIFIER"].append(str(entry["UNIQUE_IDENTIFIER"]))
            json_data["EXPERIMENT_ID"].append(int(entry["EXPERIMENT_ID"]))

        return json_data

    def save_proteins_per_experiment(self) -> list[Path]:
        """
        This function is responsible for downloading the protein data associated with each experiment ID

        :return list[str]: A list of file save locations
        """
        session = requests.Session()
        file_save_locations: list[Path] = []
        for exp_id in self._experiment_ids:
            file_extension = str(self._file_type.value).lower()
            download_file_name = f"experiment_{exp_id}.{file_extension}"
            save_location: Path = Path(
                self._output_dir, "proteomics_db", "experiment_data", download_file_name
            )
            file_save_locations.append(save_location)

            experiment_proteins = endpoints.GetProteinsPerExperiment(
                experiment_id=exp_id,
                download_format=self._file_type,
                file_save_location=save_location,
                session=session,
                print_file_size=False,
            )

            experiment_proteins.save_data()

        return file_save_locations

    def get_protein_expression(self, arguments):
        return endpoints.GetProteinExpression(**arguments)

    def save_protein_expression(
        self,
        **kwargs,
    ) -> list[Path]:
        """
        TODO: Implement multiple-downloading
        https://stackoverflow.com/a/68583332/13885200

        This function will handle the actual download of protein expression data
        At the very least, this function requires the following information:
            - Protein ID (from UniProt)
            - Experiment ID: int,
            - Download Format (i.e. DownloadFormat.JSON)

        :param kwargs: Any additional arguments as defined by the ProteomicsDB ProteinExpression API
        :return: A Path object containing the saved file location
        """
        max_workers = os.cpu_count()
        session: requests.Session = requests.Session()
        session.mount(
            "https://",
            requests.adapters.HTTPAdapter(
                pool_maxsize=max_workers, max_retries=3, pool_block=False
            ),
        )

        file_save_paths: list[Path] = []
        file_extension: str = str(self._file_type.value).lower()

        multiprocess_data: list[dict] = []
        # Build the items required for getting protein data
        print("Building dictionary for multiprocess. . .")
        for experiment_path in self._protein_per_experiment_data:
            experiment_data = self.parse_experiment_json(experiment_path)

            for i, (protein_id, experiment_id) in enumerate(
                zip(*experiment_data.values())
            ):
                download_file_name = (
                    f"experiment_{experiment_id}-protein_{protein_id}.{file_extension}"
                )
                save_path = Path(
                    self._output_dir,
                    "proteomics_db",
                    "protein_expression",
                    download_file_name,
                )

                # If the save file does not exist, continue
                if not save_path.is_file():
                    multiprocess_data.append(
                        {
                            "protein_filter": protein_id,
                            "experiment_id": experiment_id,
                            "download_format": self._file_type,
                            "file_save_location": save_path,
                            "session": session,
                            "print_file_size": False,
                        }
                    )

        # Use multithreaded requests; from: https://stackoverflow.com/a/68583332/13885200
        use_data = multiprocess_data
        print(f"Data pieces: {len(use_data)}")
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            print("Mapping workers to data")
            for response in list(executor.map(self.get_protein_expression, use_data)):
                if len(response) != 0:
                    response.save_data()
                    file_save_paths.append(response.save_path)
                    print(
                        f"Saving {response.protein_id}, size of {response.file_readable_size}"
                    )

        # Return save location if needed for further use
        return file_save_paths


def parse_arguments(argv: list[str]) -> argparse.Namespace:
    """
    This function is used to parse arguments from the command line. It is placed in a separate function to limit the amount of lines in main()
    """
    parser = argparse.ArgumentParser(
        prog="proteomics_db.py",
        description="Download and analyze proteomics data from ProteomicsDB",
        epilog="For additional help, please post questions/issues in the MADRID GitHub repo at "
        "https://github.com/HelikarLab/MADRID or email babessell@gmail.com",
    )
    parser.add_argument(
        "-c",
        "--config-file",
        required=False,
        dest="config_file",
        metavar="config_file.csv",
        help="The proteomicsdb_experiment_ids.csv file location",
    )
    parser.add_argument(
        "-d",
        "--output-dir",
        required=False,
        dest="output_dir",
        metavar="/data/proteomics/",
        help="The output directory to save proteomics results",
    )
    file_type_group = parser.add_mutually_exclusive_group(required=True)
    file_type_group.add_argument(
        "-j",
        "--json",
        dest="json",
        action="store_true",
        help="Download in JSON format",
    )
    file_type_group.add_argument(
        "-x",
        "--xml",
        dest="xml",
        action="store_true",
        help="Download in XML format",
    )

    # Get arguments
    args: argparse.Namespace = parser.parse_args(argv)

    # ----- Set Default Arguments -----
    # Args set this way to use '.' operator when accessing arguments
    # Set file_type to json or xml
    if args.json:
        args.file_type = endpoints.DownloadFormat.JSON
    else:
        args.file_type = endpoints.DownloadFormat.XML

    # Set default config_file if none provided
    if args.config_file is None:
        args.config_file = os.path.join(
            project.configs.configdir, "proteomicdb_experiment_ids.csv"
        )

    # Set default output_dir if none provided
    if args.output_dir is None:
        args.output_dir = project.configs.outputdir
    # ----- Done Setting Args -----

    return args


def main(argv: list[str]):
    # Get command line arguments
    args: argparse.Namespace = parse_arguments(argv)

    # Retrieve the experiment IDs from config file
    downloads = DownloadHandler(args)

    # Download proteins per experiment
    # downloads.save_proteins_per_experiment()


if __name__ == "__main__":
    argv: list[str] = sys.argv[1:]
    main(argv)
