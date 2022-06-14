"""
This is the main driver-file for downloading proteomics data
"""
import argparse
import Crux
import csv
from FileInformation import FileInformation
import FTPManager
import os
import sys
from pathlib import Path


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
    def __init__(self, input_csv_file: Path):
        self._input_csv_file: Path = input_csv_file
        self._data: [str, dict[str, list[str]]] = {}
        
        # Get data from CSV
        with open(self._input_csv_file, "r") as i_stream:
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
        """
        This will return a list of FTP URLs contained in the input CSV
        
        Example: ftp://ftp.my_server.com
        """
        master_urls: list[str] = []
        for cell_type in self._data.keys():
            urls = self._data[cell_type]["url"]
            master_urls.extend(urls)
        
        return urls
    
    @property
    def input_cell_types(self) -> list[str]:
        """
        This will return the cell types as defined in the input CSV file
        TODO: Match folder paths to correlate S1R1, S1R2, etc.?
        """
        cell_types: list[str] = []
        for key in self._data.keys():
            # Get the number of URLs per cell type to find the amount of cell types input
            num_urls: int = len(self._data[key]["url"])
            cell_types.extend([key] * num_urls)
        return cell_types
    
    @property
    def studies(self) -> list[str]:
        """
        This will return the replicates as defined in the input CSV file
        """
        master_studies: list[str] = []
        for cell_type in self._data.keys():
            replicates = self._data[cell_type]["replicate"]
            master_studies.extend(replicates)
        return master_studies
    
    @property
    def csv_dict(self) -> dict[str, dict[str, list[str]]]:
        """
        This function returns the data dictionary
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
            preferred_extensions: list[str] = None,

    ):
        self.file_information: list[FileInformation] = file_information
        self._csv_data: ParseCSVInput = csv_data
        self._csv_data: dict[str, dict[str, list[str]]] = csv_data.csv_dict
        self._skip_download: bool = skip_download

        # Set default value for extensions to search for
        self._preferred_extensions: list[str] = preferred_extensions
        if self._preferred_extensions is None:
            self._preferred_extensions = ["raw"]
        
        self._gather_data()
        self._new_set_replicate_numbers()
        
        if self._skip_download is False:
            self.print_download_size()
            
    def _gather_data(self):
        # Iterate through the cell type and corresponding list of URLS
        # cell_type: naiveB
        # ftp_urls: ["url_1", "url_2"]
        for i, cell_type in enumerate(self._csv_data.keys()):
            ftp_urls: list[str] = self._csv_data[cell_type]["url"]
            studies: list[str] = self._csv_data[cell_type]["study"]
            url_count = 0
        
            # Iterate through the URLs available
            for j, (url, study) in enumerate(zip(ftp_urls, studies)):
            
                # Print updates to he screen
                print(
                    f"\rParsing cell type {i + 1} of {len(self._csv_data.keys())} ({cell_type}) | {j + 1} of {len(ftp_urls)} folders navigated",
                    end="",
                    flush=True
                )
                
                ftp_data: FTPManager.Reader = FTPManager.Reader(
                    root_link=url,
                    file_extensions=self._preferred_extensions
                )
            
                urls = [url for url in ftp_data.files]
                sizes = [size for size in ftp_data.file_sizes]
                url_count += len(urls)
            
                # Iterate through all files and sizes found for url_##
                for k, (file, size) in enumerate(zip(urls, sizes)):
                    # Set the study if it was not set in the input CSV
                    if study == "":
                        replicate_name: str = f"S{j + 1}R{k + 1}"
                    else:
                        replicate_name: str = f"{study}R{k + 1}"
                
                    self.file_information.append(
                        FileInformation(
                            cell_type=cell_type,
                            download_url=file,
                            file_size=size,
                            study=study
                        )
                    )
        
            # Print number of URLs found for each cell type
            # This is done after appending because some cell types may have more than 1 root URL, and it messes up the formatting
            if url_count == 1:
                print(f" | 1 file found")
            else:
                print(f" | {url_count} files found")
            
    def print_download_size(self):
        # Print the total size to download if we must download data
        total_size: int = 0
        for information in self.file_information:
            total_size += information.file_size
        
        # Convert to MB
        total_size = total_size // 1024 ** 2
        print(f"Total size to download: {total_size} MB")
        
    def _set_replicate_numbers(self):
        """
        This function is responsible for setting the "R#" portion of each "S#R#" batch
        This is done after all the data is gathered so each batch can be collated together
        """
        replicate_num: int = 1
        for i, information in enumerate(self.file_information):
            current_info = information
            previous_info = self.file_information[i - 1] if i > 0 else None
            
            if i > 0:
                print(f"Testing: {current_info.cell_type}_{current_info.study} - {previous_info.cell_type}_{previous_info.study}")
            # Set the initial current_info replicate to 1
            if i == 0:
                replicate_value = f"R{replicate_num}"
            # If the current info cell type and study match the previous, increment the replicate by one. Otherwise reset it to 1
            elif current_info.cell_type == previous_info.cell_type and current_info.study == previous_info.study:
                replicate_num += 1
                replicate_value = f"R{replicate_num}"
            else:
                replicate_num = 1
                replicate_value = f"R{replicate_num}"
            
            current_info.set_replicate(replicate_value)
            
    def _new_set_replicate_numbers(self):
        instances: dict[str, list[FileInformation]] = {}
        for information in self.file_information:
            if information.cell_type not in instances.keys():
                instances[information.cell_type] = FileInformation.filter_instances(information.cell_type)
        
        for cell_type in instances.keys():

            replicate_num: int = 1
            for i, file_information in enumerate(instances[cell_type]):
                current_info: FileInformation = file_information
                previous_info: FileInformation = instances[cell_type][i - 1] if i > 0 else None
                
                if i > 0:
                    print(f"Testing: {current_info.cell_type}_{current_info.study} - {previous_info.cell_type}_{previous_info.study}")
                
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
        """
        This function is responsible for collecting all FileInformation objects of a given cell type
        """
        file_information_list: list[FileInformation] = []
        for information in self.file_information:
            if information.cell_type == cell_type:
                file_information_list.append(information)
                
        return file_information_list
            
    
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
        "url,cell_type\n"
        "# This is a comment\n"
        "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2022/02/PXD026140,naiveB\n"
        "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2022/02/PXD017564,m0Macro\n",
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
    # TODO: Add option to delete intermediate files (raw, mzml, sqt)

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
    
    PopulateInformation(
        file_information=file_information,
        csv_data=csv_data,
        skip_download=args.skip_download,
        preferred_extensions=args.extensions
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

    # Get the root folder of output CSV file
    root_folders: set[Path] = set([i.intensity_csv.parent for i in file_information])
    print(f"\nProtein intensities saved under:")
    for folder in root_folders:
        print(folder)
    

if __name__ == "__main__":
    args = sys.argv[1:]
    main(args)
