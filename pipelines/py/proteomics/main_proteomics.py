import argparse
import csv
import os
from pathlib import Path
import sys

# Our classes
from Database import Database
import FTPManager
import MaxQuant
import XMLEditor

# Add parent directory to path, allows us to import the "project.py" file from the parent directory
# From: https://stackoverflow.com/a/30536516/13885200
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import project


def get_urls_from_config(args: argparse.Namespace) -> list[str]:
    proteome_urls: list[str] = []
    with open(args.input_file, "r") as i_stream:
        reader = csv.DictReader(i_stream)
        for line in reader:
            proteome_urls.append(line["url"])

    return proteome_urls


def parse_args(args) -> argparse.Namespace:
    """
    This function is used to parse arguments from the command line
    """
    parser = argparse.ArgumentParser(
        prog="main_proteomics.py",
        description="Download and analyze proteomics data from proteomeXchange\n"
        "The input file should be formatted as the following example:\n"
        "\n"
        "url\n"
        "http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD026142\n"
        "http://proteomecentral.proteomexchange.org/cgi/GetDataset?ID=PXD02475",
        epilog="For additional help, please post questions/issues in the MADRID GitHub repo at "
        "https://github.com/HelikarLab/MADRID or email babessell@gmail.com",
        formatter_class=argparse.RawTextHelpFormatter,
    )
    parser.add_argument(
        "-i",
        "--input",
        required=True,
        dest="input_file",
        metavar="proteome_urls.csv",
        help="The proteomeXchange CSV file location",
    )
    parser.add_argument(
        "-o",
        "--output",
        required=False,
        dest="output_dir",
        metavar="/home/john/output_directory",
        help="The output directory to save data",
    )

    args: argparse.Namespace = parser.parse_args(args)

    # Validte the input file exists
    if not Path(args.input_file).is_file():
        print(f"Input file {args.input} does not exist!")

    # If no output directory given, create our own
    if args.output_dir is None:
        args.output_dir = os.path.join(project.configs.outputdir, "proteomics")
    os.makedirs(args.output_dir, exist_ok=True)

    return args


def build_ftp_objects(args: argparse.Namespace, urls: list[str]):
    database_list: list[Database] = []
    for url in urls:
        database_list.append(Database(proteomeXchange_url=url))


def test():
    databases: list[Database] = []
    for i in range(10):
        databases.append(
            Database(
                proteomeXchange_url=f"{i}_proteomeXchange",
                ftp_url=f"{i}_ftp",
                raw_file_count=i,
            )
        )

    for i in databases:
        print(i.ftp_url)


def main(args):
    args: argparse.Namespace = parse_args(args)
    proteome_urls: list[str] = get_urls_from_config(args)
    FTPManager.Download(proteome_urls, args)


if __name__ == "__main__":
    args = sys.argv[1:]
    main(args)
