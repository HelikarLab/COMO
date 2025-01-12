from __future__ import annotations

import argparse
import csv
import os
from pathlib import Path

from loguru import logger

from como.data_types import LogLevel
from como.proteomics import Crux, FileInformation, FTPManager
from como.utils import _log_and_raise_error


class ArgParseFormatter(argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter):
    """Use the RawTextHelpFormatter and the ArgumentDefaultsHelpFormatter in a single argparse parser()."""

    pass


class ParseCSVInput:
    def __init__(self, input_csv_file: Path):
        """Parse input CSV into two fields.

        1. proteomeXchange URLs
        2. Cell Type
        3. Replicate (optional)

        This class is meant to make it easier to access each of these things


        {
            "naiveB": {
                "url": [one, two, three],
                "replicate": [A, B, C],
            },
            "nucleophile": {
                "url": [four, five, six],
                "replicate": [D, E, F],
            }
        }
        """
        self._input_csv_file: Path = input_csv_file
        self._data: [str, dict[str, list[str]]] = {}

        # Get data from CSV
        with self._input_csv_file.open("w") as i_stream:
            reader = csv.reader(i_stream)
            next(reader)
            for line in reader:
                if line == "" or line[0][0] == "#":  # Skip 'comments' and empty lines
                    continue
                else:
                    url = line[0]
                    cell_type = line[1]
                    try:
                        study = line[2]
                    except IndexError:
                        study = ""

                    if cell_type not in self._data:
                        self._data[cell_type] = {"url": [], "study": []}

                    self._data[cell_type]["url"].append(url)
                    self._data[cell_type]["study"].append(study)

        # Convert from 'old' /pride/data/archive to 'new' /pride-archive
        for key in self._data:
            urls = self._data[key]["url"]
            for i, url in enumerate(urls):
                urls[i] = url.replace("/pride/data/archive", "/pride-archive")

    @property
    def ftp_urls(self) -> list[str]:
        """Return a list of FTP URLs contained in the input CSV.

        Example: ftp://ftp.my_server.com
        """
        master_urls: list[str] = []
        for cell_type in self._data:
            urls = self._data[cell_type]["url"]
            master_urls.extend(urls)

        return urls

    @property
    def input_cell_types(self) -> list[str]:
        """Return the cell types as defined in the input CSV file.

        TODO: Match folder paths to correlate S1R1, S1R2, etc.?
        """
        cell_types: list[str] = []
        for key in self._data:
            # Get the number of URLs per cell type to find the amount of cell types input
            num_urls: int = len(self._data[key]["url"])
            cell_types.extend([key] * num_urls)
        return cell_types

    @property
    def studies(self) -> list[str]:
        """Return the replicates as defined in the input CSV file."""
        master_studies: list[str] = []
        for cell_type in self._data:
            replicates = self._data[cell_type]["study"]
            master_studies.extend(replicates)
        return master_studies

    @property
    def csv_dict(self) -> dict[str, dict[str, list[str]]]:
        """Return the CSV information as a dictionary.

        It contains data in the following format
        {
            CELL_TYPE_1: {
                "url": ["url_one", "url_two", 'url_three', "url_four"],
                "replicate": ["S1R1", "S1R2", "S2R1", "S3R1"]
            },
            CELL_TYPE_2: {
                "url": ["url_five", "url_six", 'url_seven', "url_eight"],
                "replicate": ["S1R1", "S1R2", "S2R1", "S2R2"]
            }
        }
        """
        return self._data


class PopulateInformation:
    def __init__(
        self,
        file_information: list[FileInformation],
        csv_data: ParseCSVInput,
        skip_download: bool,
        preferred_extensions: list[str] | None = None,
    ):
        """Populate FileInformation list with data from the input CSV file."""
        self.file_information: list[FileInformation] = file_information
        self._csv_data: ParseCSVInput = csv_data
        self._csv_data: dict[str, dict[str, list[str]]] = csv_data.csv_dict
        self._skip_download: bool = skip_download

        # Set default value for extensions to search for
        self._preferred_extensions: list[str] = preferred_extensions
        if self._preferred_extensions is None:
            self._preferred_extensions = ["raw"]

        self._gather_data()
        self._set_replicate_numbers()

        if self._skip_download is False:
            self.print_download_size()

    def _gather_data(self):
        # Iterate through the cell type and corresponding list of URLS
        # cell_type: naiveB
        # ftp_urls: ["url_1", "url_2"]
        for cell_type in self._csv_data:
            ftp_urls: list[str] = self._csv_data[cell_type]["url"]
            studies: list[str] = self._csv_data[cell_type]["study"]
            url_count = 0

            # Iterate through the URLs available
            for url, study in zip(ftp_urls, studies):
                ftp_data: FTPManager.Reader = FTPManager.Reader(
                    root_link=url, file_extensions=self._preferred_extensions
                )

                urls = list(ftp_data.files)
                sizes = list(ftp_data.file_sizes)
                url_count += len(urls)

                # Iterate through all files and sizes found for url_##
                for file, size in zip(urls, sizes):
                    self.file_information.append(
                        FileInformation(cell_type=cell_type, download_url=file, file_size=size, study=study)
                    )

    def print_download_size(self):
        """Print the total size to download if we must download data."""
        total_size: int = 0
        for information in self.file_information:
            total_size += information.file_size

        # Convert to MB
        total_size = total_size // 1024**2
        logger.info(f"Total size to download: {total_size} MB")

    def _set_replicate_numbers(self):
        instances: dict[str, list[FileInformation]] = {}
        for information in self.file_information:
            if information.cell_type not in instances:
                instances[information.cell_type] = FileInformation.filter_instances(information.cell_type)

        for cell_type in instances:
            replicate_num: int = 1
            for i, file_information in enumerate(instances[cell_type]):
                current_info: FileInformation = file_information
                previous_info: FileInformation = instances[cell_type][i - 1] if i > 0 else None

                # Do not modify the replicate value if we are on the very first iteration of this cell type
                if i == 0:
                    pass
                # If the current_info cell type and study match the previous, increment the replicate_num by one
                elif current_info.cell_type == previous_info.cell_type and current_info.study == previous_info.study:
                    replicate_num += 1
                else:
                    replicate_num: int = 1

                replicate_value: str = f"R{replicate_num}"
                current_info.set_replicate(replicate_value)

    def _collect_cell_type_information(self, cell_type: str) -> list[FileInformation]:
        """Collect all FileInformation objects of a given cell type."""
        return [information for information in self.file_information if information.cell_type == cell_type]


def parse_args() -> argparse.Namespace:
    """Parse arguments from the command line."""
    parser = argparse.ArgumentParser(
        prog="proteomics_preprocess.py",
        description="Download and analyze proteomics data from proteomeXchange\n"
        "Comments can be added to the csv file by starting a line with a '#'\n"
        "The input file should be formatted as the following example:\n"
        "\n"
        "url,cell_type,study\n"
        "# This is a comment\n"
        "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2022/02/PXD026140,naiveB,S1\n"
        "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2022/02/PXD017564,m0Macro,S1\n"
        "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2022/02/PXD017987,naiveB,S2\n",
        epilog="For additional help, please post questions/issues in the MADRID GitHub repo at "
        "https://github.com/HelikarLab/MADRID or email babessell@gmail.com",
        formatter_class=ArgParseFormatter,
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        dest="input_csv",
        metavar="/home/USER/data/proteomics_urls.csv",
        help="The proteomeXchange CSV file location",
    )
    parser.add_argument(
        "-d",
        "--database",
        required=True,
        dest="database",
        metavar="/home/USER/data/database_file.fasta",
        help="The fasta database to search for protein identification",
    )
    parser.add_argument(
        "-e",
        "--extensions",
        dest="extensions",
        required=False,
        default="raw",
        help="A list of file extensions to download from the FTP server",
        metavar="raw,txt,another_extension",
    )
    parser.add_argument(
        "--skip-download",
        required=False,
        dest="skip_download",
        default=False,
        type=bool,
        help="If this action is passed in, FTP data will not be downloaded. "
        "This assumes you have raw data under the folder specified by the option '--ftp-out-dir'",
    )
    parser.add_argument(
        "--skip-mzml",
        required=False,
        dest="skip_mzml",
        type=bool,
        default=False,
        help="If this action is passed in, files will not be converted from RAW to mzML format. "
        "This assumes you have mzML files under the folder specified by the option '--mzml-out-dir'. "
        "This will continue the workflow from SQT file creation -> CSV ion intensity creation. "
        "If this option is passed in, FTP data will also not be downloaded.",
    )
    parser.add_argument(
        "--skip-sqt",
        required=False,
        dest="skip_sqt",
        type=bool,
        default=False,
        help="If this action is passed in, SQT files will not be created. "
        "This assumes you have SQT files under the folder specified by the option '--sqt-out-dir'. "
        "This will only read data from SQT files and create a CSV ion intensity file. "
        "If this option is passed in, FTP data will not be downloaded, "
        "RAW files will not be converted, and SQT files will not be created.",
    )
    parser.add_argument(
        "-c",
        "--cores",
        required=False,
        dest="core_count",
        metavar="cores",
        default=os.cpu_count() // 2,
        help="This is the number of threads to use for downloading files. "
        "It will default to the minimum of: half the available CPU cores available, "
        "or the number of input files found. "
        "It will not use more cores than necessary. "
        "Options are an integer or 'all' to use all available cores. "
        "Note: Downloading will use a MAX of 2 threads at once, "
        "as some FTP servers do not work well with multiple connections from the same IP address at once.",
    )
    # TODO: Add option to delete intermediate files (raw, mzml, sqt)

    args: argparse.Namespace = parser.parse_args()
    args.extensions = args.extensions.split(",")

    # Validte the input file exists
    if not Path(args.input_csv).is_file():
        _log_and_raise_error(f"Input file {args.input} does not exist!", error=FileNotFoundError, level=LogLevel.ERROR)

    if args.core_count == "all":
        args.core_count = os.cpu_count()
    elif not str(args.core_count).isdigit():
        _log_and_raise_error(
            f"Invalid option '{args.core_count}' for option '--cores'. Enter an integer or 'all' to use all cores",
            error=ValueError,
            level=LogLevel.ERROR,
        )

    else:
        args.core_count = int(args.core_count)
        if args.core_count > os.cpu_count():
            logger.info(
                f"{args.core_count} cores not available, system only has {os.cpu_count()} cores. "
                f"Setting '--cores' to {os.cpu_count()}"
            )
            args.core_count = os.cpu_count()

    return args


def _main():
    file_information: list[FileInformation] = []
    args: argparse.Namespace = parse_args()
    csv_data = ParseCSVInput(args.input_csv)

    """
        This comment is for the logic surrounding "skipping" a step in the workflow
        1. skip_download (download FTP data - FTPManager)
        2. skip_mzml_conversion (Convert raw to mzML - Crux)
        3. skip_sqt_creation (Convert mzML to SQT - Crux)

        If args.skip_sqt is True, do not perform steps 1, 2, or 3
        If args.skip_mzml is True, do not perform step 1 or 2
        If args.skip_download is True, do not perform step 1

        Ultimately, this results in if-statements that look like:
            if __ is False:
                do_tasks
        Because we are performing tasks if the "skip" is False
        """
    if args.skip_sqt:
        args.skip_mzml = True
        args.skip_download = True
    elif args.skip_mzml:
        args.skip_download = True

    # Populate the file_information list
    PopulateInformation(
        file_information=file_information,
        csv_data=csv_data,
        skip_download=args.skip_download,
        preferred_extensions=args.extensions,
    )

    # Download data if we should not skip anything
    if args.skip_download is False:
        # Start the download of FTP data
        FTPManager.Download(
            file_information=file_information,
            core_count=args.core_count,
        )

    if args.skip_mzml is False:
        # Convert raw to mzML and then create SQT files
        Crux.RAWtoMZML(
            file_information=file_information,
            core_count=args.core_count,
        )

    if args.skip_sqt is False:
        # Convert mzML to SQT
        Crux.MZMLtoSQT(
            file_information=file_information,
            fasta_database=args.database,
            core_count=args.core_count,
        )

    # Create CSV file from SQT files
    Crux.SQTtoCSV(
        file_information=file_information,
        core_count=args.core_count,
    )


if __name__ == "__main__":
    _main()
