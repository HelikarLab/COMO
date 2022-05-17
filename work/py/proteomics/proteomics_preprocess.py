"""
This is the main driver-file for downloading proteomics data
"""
import argparse
import csv
import os
from pathlib import Path
import sys

# Our classes
import Crux
import FTPManager


# Add parent directory to path, allows us to import the "project.py" file from the parent directory
# From: https://stackoverflow.com/a/30536516/13885200
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import project


class ArgParseFormatter(
    argparse.RawTextHelpFormatter, argparse.ArgumentDefaultsHelpFormatter
):
    """
    This allows us to use the RawTextHelpFormatter and the ArgumentDefaultsHelpFormatter in a single argparse parser()
    """

    pass


class Defaults:
    """
    This class is responsible for holding default values for argument parsing
    """
    # TODO: Move "/work/data/{outputdir}" to "/work/data/results/CELL_TYPE/proteomics"
    # TODO: Output "output-csv" to data-matrices/protein_abundance_matrix_CELL_TYPE.csv
    ftp_links: Path = Path(project.configs.configdir, "proteomeXchange_urls.csv")
    output_intensity_csv: Path = Path(project.configs.outputdir, "proteomics", "proteomics_intensity.csv")  # fmt: skip
    raw_file_directory: Path = Path(project.configs.outputdir, "proteomics", "raw_files")  # fmt: skip
    mzml_directory: Path = Path(project.configs.outputdir, "proteomics", "mzml")
    sqt_directory: Path = Path(project.configs.outputdir, "proteomics", "sqt")
    fasta_database_file: Path = Path(project.configs.datadir, "human_proteome_UP000005640_database.fasta")  # fmt: skip
    core_count: int = os.cpu_count() // 2
    
    @property
    def collect_raw_files(self) -> list[Path]:
        """
        This function is responsible for recursively collecting all RAW files from raw_file_directory
        """
        raw_files: list[Path] = []
        for file in self.raw_file_directory.rglob("*"):
            if file.suffix == ".raw":
                raw_files.append(file)
        return raw_files
    
    @property
    def collect_sqt_files(self) -> list[Path]:
        """
        This function is responsible for recursively collecting all SQT files from sqt_directory
        """
        sqt_files: list[Path] = []
        # Recursively collect all SQT files from sqt_directory
        for file in self.sqt_directory.rglob("*"):
            if file.suffix == ".sqt":
                sqt_files.append(file)
        return sqt_files
    
    @property
    def collect_mzml_files(self) -> list[Path]:
        """
        This function is responsible for recursively collecting all mzML files from the mzml_directory
        """
        mzml_files: list[Path] = []
        for file in self.mzml_directory.rglob("*"):
            if file.suffix == ".mzml":
                mzml_files.append(file)
        return mzml_files


class ParseInputCSV:
    """
    This class is responsible for parsing the input CSV into two fields
    1. proteomeXchange URLs
    2. Cell Type
    
    This class is meant to make it easier to access each of these things
    """
    def __init__(self, input_csv_file: Path):
        self._input_csv_file: Path = input_csv_file
        self._proteome_urls: list[str] = []
        self._cell_types: list[str] = []
        
        with open(self._input_csv_file, "r") as i_stream:
            reader = csv.reader(i_stream)
            header = next(reader)
            for line in reader:
                if line[0][0]== "#":  # Skip 'comments'
                    continue
                else:
                    self._proteome_urls.append(line[0])
                    self._cell_types.append(line[1])

    @property
    def ftp_urls(self) -> list[str]:
        """
        This will return a list of FTP URLs contained in the input CSV
        
        Example: ftp://ftp.my_server.com
        """
        return self._proteome_urls
    
    @property
    def cell_types(self) -> list[str]:
        """
        This will return the cell types as defined in the input CSV file
        """
        return self._cell_types


def parse_args(args: list[str]) -> argparse.Namespace:
    """
    This function is used to parse arguments from the command line

    :param args: The list of arguments collected from the command line
    """

    parser = argparse.ArgumentParser(
        prog="proteomics_preprocess.py",
        description="Download and analyze proteomics data from proteomeXchange\n"
        "Comments can be added to the csv file by starting a line with a '#'\n"
        "The input file should be formatted as the following example:\n"
        "\n"
        "url\n"
        "# This is a comment\n"
        "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2022/02/PXD026140\n"
        "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2022/02/PXD017564\n",
        epilog="For additional help, please post questions/issues in the MADRID GitHub repo at "
        "https://github.com/HelikarLab/MADRID or email babessell@gmail.com",
        formatter_class=ArgParseFormatter,
    )
    parser.add_argument(
        "-i",
        "--input",
        required=False,
        dest="input_csv",
        default=Defaults.ftp_links,
        metavar="proteomics_urls.csv",
        help="The proteomeXchange CSV file location",
    )
    parser.add_argument(
        "-o",
        "--output-csv",
        required=False,
        default=Defaults.output_intensity_csv,
        dest="output_csv",
        metavar="proteomics_intensity.csv",
        help="The location to write the output intensity CSV file",
    )
    parser.add_argument(
        "-f",
        "--ftp-out-dir",
        required=False,
        dest="ftp_output_dir",
        default=Defaults.raw_file_directory,
        metavar="ftp_directory",
        help="The output directory to raw mass spectrometry data",
    )
    parser.add_argument(
        "-m",
        "--mzml-out-dir",
        required=False,
        dest="mzml_output_dir",
        default=Defaults.mzml_directory,
        metavar="mzml_output_dir",
        help="The output directory to save mzML files. These are generated by converting raw files",
    )
    parser.add_argument(
        "-s",
        "--sqt-out-dir",
        required=False,
        dest="sqt_output_dir",
        default=Defaults.sqt_directory,
        metavar="sqt_output_directory",
        help="The output directory to save SQT files. These are generated by analyzing mzML files",
    )
    parser.add_argument(
        "-d",
        "--database-search",
        required=False,
        dest="database",
        default=Defaults.fasta_database_file,
        metavar="database_file.fasta",
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
        "-p",
        "--print-only",
        required=False,
        dest="print_only",
        action="store_true",
        help="Should the files found within the FTP links only be printed to the console? Enabling this option will not download any files.",
    )
    parser.add_argument(
        "--skip-download",
        required=False,
        dest="skip_download",
        action="store_true",
        default=False,
        help="If this action is passed in, FTP data will not be downloaded.\nThis assumes you have raw data under the folder specified by the option '--ftp-out-dir'",
    )
    parser.add_argument(
        "--skip-mzml",
        required=False,
        dest="skip_mzml",
        action="store_true",
        default=False,
        help="If this action is passed in, files will not be converted from RAW to mzML format.\n"
        "This assumes you have mzML files under the folder specified by the option '--mzml-out-dir'.\n"
        "This will continue the workflow from SQT file creation -> CSV ion intensity creation.\n"
        "If this option is passed in, FTP data will also not be downloaded.",
    )
    parser.add_argument(
        "--skip-sqt",
        required=False,
        dest="skip_sqt",
        action="store_true",
        default=False,
        help="If this action is passed in, SQT files will not be created. This assumes you have SQT files under the folder specified by the option '--sqt-out-dir'.\n"
        "This will only read data from SQT files and create a CSV ion intensity file.\n"
        "If this option is passed in, FTP data will not be downloaded, RAW files will not be converted, and SQT files will not be created.",
    )
    parser.add_argument(
        "-c",
        "--cores",
        required=False,
        dest="core_count",
        metavar="cores",
        default=Defaults.core_count,
        help="This is the number of threads to use for downloading files. It will default to the minimum of: half the available CPU cores available, or the number of input files found.\nIt will not use more cores than necessary\nOptions are an integer or 'all' to use all available cores",
    )
    # TODO: Add "keep" flag to optionally keep the downloaded intermediate files (raw_files, mzml, sqt)

    args: argparse.Namespace = parser.parse_args(args)
    args.extensions = args.extensions.split(",")

    # Validte the input file exists
    if not Path(args.input_csv).is_file():
        print(f"Input file {args.input} does not exist!")
        raise FileNotFoundError

    try:
        # Try to get an integer, fails if "all" input
        args.core_count = int(args.core_count)
        if args.core_count > os.cpu_count():
            print(f"{args.core_count} cores not available, system only has {os.cpu_count()} cores. Setting '--cores' to {os.cpu_count()}")  # fmt: skip
            args.core_count = os.cpu_count()
    except ValueError:
        if args.core_count == "all":
            args.core_count = os.cpu_count()
        else:
            print(f"Invalid option '{args.core_count}' for option '--cores'. Enter an integer or 'all' to use all cores")  # fmt: skip
            exit(2)

    # We have a default output directory set (default_output)
    # Because of this, even if the user inputs their own output directory, we can just make the output
    # It will be located at "pipelines/output/proteomics"
    os.makedirs(args.ftp_output_dir, exist_ok=True)
    os.makedirs(args.mzml_output_dir, exist_ok=True)
    os.makedirs(args.sqt_output_dir, exist_ok=True)

    return args


def main(args: list[str]):
    """
    This is the main driver function

    :param args: The list of arguments collected from the command line
    """
    args: argparse.Namespace = parse_args(args)

    # Create variables so they are always defined
    # This is required in case a "skip-..." option is passed in
    csv_parser = ParseInputCSV(args.input_csv)
    ftp_links: list[str] = csv_parser.ftp_urls
    csv_cell_types: list[str] = csv_parser.cell_types
    raw_file_paths: list[Path] = Defaults().collect_raw_files
    sqt_file_paths = Defaults().collect_sqt_files
    mzml_file_paths = Defaults().collect_mzml_files
    raw_file_cell_types: list[str] = []
    
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
        
    # If skipping download, must define cell types
    if args.skip_download is True:
        for file in raw_file_paths:
            prefix = file.stem.split("_")
            prefix = f"{prefix[0]}_{prefix[1]}"
            raw_file_cell_types.append(prefix)

    # Download data if we should not skip anything
    if args.skip_download is False:
        # Start the download of FTP data
        ftp_manager = FTPManager.Download(
            ftp_links=ftp_links,
            output_dir=args.ftp_output_dir,
            cell_types=csv_cell_types,
            preferred_extensions=args.extensions,
            skip_download=args.skip_download,
            core_count=args.core_count,
        )
        raw_file_paths: list[Path] = ftp_manager.raw_files
        raw_file_cell_types = ftp_manager.collected_cell_types

    if args.skip_mzml is False:
        # Convert raw to mzML and then create SQT files
        raw_to_mzml = Crux.RAWtoMZML(
            raw_file_input=raw_file_paths,
            mzml_output_dir=args.mzml_output_dir,
            core_count=args.core_count,
        )
        mzml_file_paths: list[Path] = raw_to_mzml.mzml_file_paths

    if args.skip_sqt is False:
        # Convert mzML to SQT
        mzml_to_sqt = Crux.MZMLtoSQT(
            mzml_file_paths=mzml_file_paths,
            fasta_database=args.database,
            sqt_output_dir=args.sqt_output_dir,
            cell_types=raw_file_cell_types,
            core_count=args.core_count,
        )
        sqt_file_paths: list[Path] = mzml_to_sqt.sqt_file_paths

    # Create CSV file from SQT files
    conversion_manager = Crux.SQTtoCSV(
        cell_types=raw_file_cell_types,
        sqt_input_files=sqt_file_paths,
        output_file=args.output_csv
    )

    print(f"Gene ID output saved at {conversion_manager.output_csv}")


if __name__ == "__main__":
    args = sys.argv[1:]
    main(args)
