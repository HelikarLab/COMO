"""
This is the main driver-file for downloading proteomics data
"""

import argparse
import csv
import os
from pathlib import Path
import sys

# Our classes
from Database import Database
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


def get_ftp_urls(input_csv_file: str) -> list[str]:
    """
    This function is responsible for collecting the FTP URLs from the input csv file

    :param input_csv_file: The input file to read from
    """
    proteome_urls: list[str] = []
    with open(input_csv_file, "r") as i_stream:
        reader = csv.DictReader(i_stream)
        for line in reader:
            url: str = line["url"]

            # If the line has a "#", do not add it. Consider "#" as comments
            if url[0] != "#":
                proteome_urls.append(url)

    return proteome_urls


def parse_args(args: list[str]) -> argparse.Namespace:
    """
    This function is used to parse arguments from the command line

    :param args: The list of arguments collected from the command line
    """
    # Define several default values
    default_csv: Path = Path(project.configs.configdir, "proteomeXchange_urls.csv")
    default_output: Path = Path(project.configs.outputdir, "proteomics")

    parser = argparse.ArgumentParser(
        prog="main_proteomics.py",
        description="Download and analyze proteomics data from proteomeXchange\n"
        "Comments can be added to the csv file by starting a line with a '#'\n"
        "The input file should be formatted as the following example:\n"
        "\n"
        "url\n"
        "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2022/02/PXD026140\n"
        "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2022/02/PXD017564\n"
        "# This is a comment",
        epilog="For additional help, please post questions/issues in the MADRID GitHub repo at "
        "https://github.com/HelikarLab/MADRID or email babessell@gmail.com",
        formatter_class=ArgParseFormatter,
    )
    parser.add_argument(
        "-i",
        "--input",
        required=False,
        dest="input_csv",
        default=default_csv,
        metavar=str(default_csv),
        help="The proteomeXchange CSV file location",
    )
    parser.add_argument(
        "-o",
        "--out-dir",
        required=False,
        dest="output_dir",
        default=default_output,
        metavar=str(default_output),
        help="The output directory to save data",
    )
    # Accept a list of arguments. From: https://stackoverflow.com/a/15753721/13885200
    parser.add_argument(
        "-e",
        "--extensions",
        dest="extensions",
        required=False,
        default="raw,fasta",
        help="A list of file extensions to download",
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

    args: argparse.Namespace = parser.parse_args(args)
    args.extensions = args.extensions.split(",")

    # Validte the input file exists
    if not Path(args.input_csv).is_file():
        print(f"Input file {args.input} does not exist!")
        exit(1)

    # We have a default output directory set (default_output)
    # Because of this, even if the user inputs their own output directory, we can just make the output
    # It will be located at "pipelines/output/proteomics"
    os.makedirs(args.output_dir, exist_ok=True)

    return args


def main(args: list[str]):
    """
    This is the main driver function

    :param args: The list of arguments collected from the command line
    """
    args: argparse.Namespace = parse_args(args)
    ftp_links: list[str] = get_ftp_urls(args.input_csv)
    FTPManager.Download(
        ftp_links=ftp_links,
        output_dir=args.output_dir,
        print_only=args.print_only,
        preferred_extensions=args.extensions,
    )


if __name__ == "__main__":
    args = sys.argv[1:]
    main(args)
