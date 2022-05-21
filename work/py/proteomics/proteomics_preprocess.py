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
from Defaults import DefaultValues
from FileInformation import FileInformation

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


class ParseCSVInput:
    """
    This class is responsible for parsing the input CSV into two fields
    1. proteomeXchange URLs
    2. Cell Type
    
    This class is meant to make it easier to access each of these things
    """
    def __init__(self, input_csv_file: Path):
        self._input_csv_file: Path = input_csv_file
        self._urls: list[str] = []
        self._input_cell_types: list[str] = []
        
        # Get data from CSV
        with open(self._input_csv_file, "r") as i_stream:
            reader = csv.reader(i_stream)
            header = next(reader)
            for line in reader:
                if line[0][0] == "#":  # Skip 'comments'
                    continue
                else:
                    self._urls.append(line[0])
                    self._input_cell_types.append(line[1])

        # Convert from 'old' /pride/data/archive to 'new' /pride-archive
        for i, url in enumerate(self._urls):
            self._urls[i] = url.replace("/pride/data/archive", "/pride-archive")
        
    @property
    def ftp_urls(self) -> list[str]:
        """
        This will return a list of FTP URLs contained in the input CSV
        
        Example: ftp://ftp.my_server.com
        """
        return self._urls
    
    @property
    def input_cell_types(self) -> list[str]:
        """
        This will return the cell types as defined in the input CSV file
        """
        return self._input_cell_types


class PopulateInformation:
    def __init__(
            self,
            file_information: list[FileInformation],
            csv_data: ParseCSVInput
    ):
        self.file_information: list[FileInformation] = file_information
        self._csv_data: ParseCSVInput = csv_data
        
        for i, (root_url, cell_type) in enumerate(zip(self._csv_data.ftp_urls, self._csv_data.input_cell_types)):
            print(f"Collecting file information for cell type: {cell_type} ({i+1} / {len(self._csv_data.input_cell_types)})")
            ftp_files: FTPManager.Reader = FTPManager.Reader(root_link=root_url, file_extensions=["raw"])
        
            for j, (file, size) in enumerate(zip(ftp_files.files, ftp_files.file_sizes)):
                replicate_name: str = f"S1R{j + 1}"
                self.file_information.append(
                    FileInformation(
                        cell_type=cell_type,
                        download_url=file,
                        file_size=size,
                        replicate=replicate_name,
                    )
                )
                
        
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
        required=True,
        dest="input_csv",
        metavar="/home/USER/data/proteomics_urls.csv",
        help="The proteomeXchange CSV file location",
    )
    parser.add_argument(
        "-d",
        "--database-search",
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
        default=os.cpu_count() // 2,
        help="This is the number of threads to use for downloading files. It will default to the minimum of: half the available CPU cores available, or the number of input files found.\nIt will not use more cores than necessary\nOptions are an integer or 'all' to use all available cores",
    )
    # parser.add_argument(
    #     "--delete",
    #     required=False,
    #     dest="delete",
    #     action="store_true",
    #     default=False,
    #     help="If this action is passed in, all intermediate files will be deleted. This includes: raw files, mzML files, and SQT files.\n"
    #          "Only the final proteomics intensities CSV file will be kept.",
    # )

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

    return args
        

def main(args: list[str]):
    """
    This is the main driver function

    :param args: The list of arguments collected from the command line
    """
    file_information: list[FileInformation] = []
    args: argparse.Namespace = parse_args(args)
    csv_data = ParseCSVInput(args.input_csv)
    
    PopulateInformation(file_information, csv_data)
    
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

    # Download data if we should not skip anything
    if args.skip_download is False:
        # Start the download of FTP data
        ftp_manager = FTPManager.Download(
            file_information=file_information,
            core_count=args.core_count,
        )

    if args.skip_mzml is False:
        # Convert raw to mzML and then create SQT files
        raw_to_mzml = Crux.RAWtoMZML(
            file_information=file_information,
            core_count=args.core_count,
        )

    if args.skip_sqt is False:
        # Convert mzML to SQT
        mzml_to_sqt = Crux.MZMLtoSQT(
            file_information=file_information,
            fasta_database=args.database,
            core_count=args.core_count,
        )

    # Create CSV file from SQT files
    conversion_manager = Crux.SQTtoCSV(
        file_information=file_information,
    )

    print(f"Protein intensities saved")
    

if __name__ == "__main__":
    args = sys.argv[1:]
    main(args)
