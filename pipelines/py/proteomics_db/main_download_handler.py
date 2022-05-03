import argparse
import csv
import endpoints
import json
import os
from pathlib import Path
import requests
import requests.adapters
import sys
from tqdm import tqdm
from pprint import pprint


# Add parent directory to path, allows us to import the "project.py" file
# From: https://stackoverflow.com/a/30536516/13885200
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import project


class ProteinData:
    def __init__(self, args: argparse.Namespace, **kwargs):
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
        self._experiment_ids: list[int] = self._read_experiment_ids()

        # fmt: off
        # Get proteins per experiment
        self._protein_per_experiment_data: list[endpoints.GetProteinsPerExperiment] = self._get_proteins_per_experiment(**kwargs)

        # Get protein expression
        self._protein_expression_data: list[endpoints.GetProteinExpression] = self._get_protein_expression(**kwargs)
        # fmt: on

    def _read_experiment_ids(self) -> list[int]:
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
                    print(f"Unable to convert experiment ID value {row['experiment_id']} to an integer! Validate the experiment IDs found in {self._config_file}.")  # fmt: skip

        return experiment_ids

    def _parse_experiment_json(self, experiment_path: Path) -> dict:
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
            # Create a set of unique items
            if entry["UNIQUE_IDENTIFIER"] not in json_data["UNIQUE_IDENTIFIER"]:
                json_data["UNIQUE_IDENTIFIER"].append(str(entry["UNIQUE_IDENTIFIER"]))
                json_data["EXPERIMENT_ID"].append(int(entry["EXPERIMENT_ID"]))

        return json_data

    def _get_proteins_per_experiment(
        self,
        session: requests.Session() = None,
        **kwargs,
    ) -> list[endpoints.GetProteinsPerExperiment]:
        """
        This function is responsible for downloading the protein data associated with each experiment ID

        :param kwargs: Any additional parameters as defined by the endpoints.GetProteinsPerExperiment class
        :return list[str]: A list of file save locations
        """
        # Define default parameters for this function
        if session is None:
            session = requests.Session()

        experiment_data: list[endpoints.GetProteinsPerExperiment] = []

        for exp_id in self._experiment_ids:

            file_extension = str(self._file_type.value).lower()
            download_file_name = f"experiment_{exp_id}.{file_extension}"
            save_path: Path = Path(self._output_dir, "proteomics_db", "experiment_data", download_file_name)  # fmt: skip

            experiment_proteins = endpoints.GetProteinsPerExperiment(
                experiment_id=exp_id,
                download_format=self._file_type,
                file_save_location=save_path,
                session=session,
                **kwargs,
            )
            experiment_proteins.save_data()

            experiment_data.append(experiment_proteins)

        return experiment_data

    def _get_protein_expression(
        self,
        session: requests.Session() = None,
        **kwargs,
    ) -> list[endpoints.GetProteinExpression]:
        """
        TODO: Implement multi-threaded downloading
        https://stackoverflow.com/a/68583332/13885200

        This function will handle the actual download of protein expression data
        At the very least, this function requires the following information:
            - Protein ID (from UniProt)
            - Experiment ID: int,
            - Download Format (i.e. DownloadFormat.JSON)

        :param kwargs: Any additional arguments as defined by the endpoints.GetProteinExpression class
        :return: A Path object containing the saved file location
        """
        if session is None:
            session = requests.Session()

        protein_expression_data: list[endpoints.GetProteinExpression] = []
        print(f"{len(self._protein_per_experiment_data)} experiment(s) collected")

        for experiment in self._protein_per_experiment_data:
            print(f"Working on experiment ID: {experiment.id}")

            experiment_data = self._parse_experiment_json(experiment.save_path)
            file_extension = str(self._file_type.value).lower()

            print(f"{len(experiment_data['UNIQUE_IDENTIFIER'])} proteins to get.")

            # Pair the items in each list found in experiment_data
            for i, (protein_id, experiment_id) in tqdm(
                enumerate(zip(*experiment_data.values()))
            ):
                download_file_name = f"experiment_{experiment_id}-protein_{protein_id}.{file_extension}"  # fmt: skip

                save_path = Path(
                    self._output_dir,
                    "proteomics_db",
                    "protein_expression",
                    download_file_name,
                )

                # Only save data if the file does not exist
                if not save_path.is_file():
                    protein_data = endpoints.GetProteinExpression(
                        protein_filter=protein_id,
                        experiment_id=experiment_id,
                        file_save_location=save_path,
                        download_format=self._file_type,
                        **kwargs,
                    )

                else:
                    protein_data = self._load_json_data(save_path)

                # Only save protein data if data exists
                if len(protein_data) != 0:
                    protein_data.save_data()
                    protein_expression_data.append(protein_data)

        return protein_expression_data

    def _load_json_data(self, file_path: Path) -> list[dict]:
        """
        This function is responsible for loading the JSON data of experiment or protein data which is already downloaded
        """
        lines = open(file_path, "r").read()
        json_object = json.loads(lines)
        return json_object

    @property
    def experiment_ids(self) -> list[int]:
        return self._experiment_ids

    @property
    def protein_per_experiment_data(self) -> list[endpoints.GetProteinsPerExperiment]:
        return self._protein_per_experiment_data

    @property
    def protein_expression_data(self) -> list[endpoints.GetProteinExpression]:
        return self._protein_expression_data


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
        args.config_file = os.path.join(project.configs.configdir, "proteomicdb_experiment_ids.csv")  # fmt: skip

    # Set default output_dir if none provided
    if args.output_dir is None:
        args.output_dir = project.configs.outputdir
    # ----- Done Setting Args -----

    return args


def main(argv: list[str]):
    # Get command line arguments
    args: argparse.Namespace = parse_arguments(argv)

    # Gather all protein information
    protein_data = ProteinData(args=args, print_file_size=False)


if __name__ == "__main__":
    argv: list[str] = sys.argv[1:]
    main(argv)
