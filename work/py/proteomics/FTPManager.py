"""
This file is responsible downloading data found at FTP links

TODO: Find a way to mark a file as "downloaded"
    - Keep a list of file names in a ".completd" hidden folder?
    -
"""
from ftplib import FTP
import multiprocessing
from multiprocessing.sharedctypes import Synchronized
import numpy as np
from pathlib import Path
import time
from urllib.parse import urlparse
from urllib.parse import ParseResult


class Download:
    def __init__(
        self,
        ftp_links: list[str],
        output_dir: Path,
        cell_types: list[str],
        preferred_extensions: list[str],
        skip_download: bool,
        core_count: int = 1,
    ):
        """
        This function is responsible for downloading items from the FTP urls passed into it

        :param ftp_links: A list of FTP links to download from
        :param output_dir: The directory to save data into
        :param preferred_extensions: The file extension to look for. If not set, it will recursively find every file
        """
        self._input_cell_types: list[str] = cell_types
        self._ftp_links: list[str] = ftp_links
        self._output_dir: Path = output_dir
        self._extensions: tuple[str] = tuple(preferred_extensions)
        self._core_count: int = core_count
        self._raw_file_cell_types: list[str] = []

        # Create a manager so each process can append data to variables
        # From: https://stackoverflow.com/questions/67974054
        self._save_locations: list[Path] = multiprocessing.Manager().list()
        self._download_counter: Synchronized = multiprocessing.Value("i", 1)
        self._finished_counter: Synchronized = multiprocessing.Value("i", 1)

        # These values are modified by self.get_files_to_download()
        self._files_to_download: list[str] = []
        self._file_sizes: list[int] = []

        # Find files to download
        self.get_files_to_download()
        self.download_data_wrapper()

    def get_files_to_download(self):
        """
        This function is responsible for asynchronously downloading FTP files from the self._urls list
        """

        # Connect to the FTP server
        print("Collecting files to download...", end=" ", flush=True)
        for cell_type, url in enumerate(self._ftp_links):
            parsed_url: ParseResult = urlparse(url)
            scheme: str = parsed_url.scheme
            host: str = parsed_url.hostname
            folder: str = parsed_url.path

            # Attempt to connect to host
            connection_successful: bool = False
            max_attempts: int = 5
            attempts: int = 1
            while not connection_successful:
                try:
                    client: FTP = FTP(host=host, user="anonymous", passwd="guest", timeout=5)
                    client.set_pasv(False)
                    connection_successful = True
                except ConnectionResetError:
                    # If connection is reset, wait a bit and try again
                    print(
                        f"\nConnection reset. Retrying... [{max_attempts - attempts} attempt(s) remaining]",
                        end="",
                        flush=True,
                    )
                    time.sleep(1)
                    attempts += 1
                finally:
                    if attempts > 5:
                        raise ConnectionError(f"Could not connect to {host}. Please check your internet connection.")  # fmt: skip

            for file_path in client.nlst(folder):
                if file_path.endswith(self._extensions):
                    download_url: str = f"{scheme}://{host}{file_path}"
                    file_size: int = client.size(file_path)
                    self._raw_file_cell_types.append(self._input_cell_types[cell_type])
                    self._files_to_download.append(download_url)
                    self._file_sizes.append(file_size)
                    
            client.quit()

        print(f"{len(self._files_to_download)} files found\n")

    def download_data_wrapper(self):
        """
        This function is responsible for using multiprocessing to download as many files at once in parallel
        """
        # Calculate the number of cores ot use
        # We are going to use half the number of cores, OR the number of files to download, whichever is less
        # Otherwise if user specified the number of cores, use that

        # Split list into chunks to process separately
        file_chunks: list[str] = np.array_split(
            self._files_to_download, self._core_count
        )
        file_size_chunks: list[int] = np.array_split(self._file_sizes, self._core_count)

        # Create a list of jobs
        jobs: list[multiprocessing.Process] = []

        for i, (files, sizes) in enumerate(zip(file_chunks, file_size_chunks)):
            # Append a job to the list
            job = multiprocessing.Process(
                target=self.download_data,
                args=(files, sizes,),
            )
            jobs.append(job)

        [job.start() for job in jobs]  # Start jobs
        [job.join() for job in jobs]  # Wait for jobs to finish
        [job.terminate() for job in jobs]  # Terminate jobs

    def download_data(
        self,
        files_to_download: list[str],
        file_sizes: list[int],
    ):

        # Start processing the files
        
        for i, (file, size) in enumerate(zip(files_to_download, file_sizes)):
            # Define FTP-related variables
            parsed_url = urlparse(file)
            host = parsed_url.hostname
            path = parsed_url.path

            # Convert file_size from byte to MB
            size_mb: int = size // (1024**2)

            # Determine the file name and save path
            replicate_name: str = f"{self._raw_file_cell_types[i]}_S1R{i + 1}"
            file_name: str = f"{replicate_name}_{Path(path).name}"
            save_path: Path = Path(self._output_dir, file_name)
            self._save_locations.append(save_path)

            # Connect to the host, login, and download the file
            client: FTP = FTP(host=host, user="anonymous", passwd="guest")

            # Get the lock and print file info
            self._download_counter.acquire()
            print(f"Started {self._download_counter.value:02d} / {len(self._files_to_download):02d} ({size_mb} MB) - {file_name}")  # fmt: skip
            self._download_counter.value += 1
            self._download_counter.release()

            # Download the file
            client.retrbinary(f"RETR {path}", open(save_path, "wb").write)
            client.quit()

            # Show finished status
            self._finished_counter.acquire()
            print(f"\tFinished {self._finished_counter.value:02d} / {len(self._files_to_download):02d}")  # fmt: skip
            self._finished_counter.value += 1
            self._finished_counter.release()

    @property
    def ftp_links(self) -> list[str]:
        """
        This function returns the FTP links used to download data
        """
        return self._ftp_links

    @property
    def output_dir(self) -> Path:
        """
        This function returns the output directory used to save data
        """
        return self._output_dir

    @property
    def raw_files(self) -> list[Path]:
        """
        This function returns a list of downloaded file locations
        """
        return self._save_locations
    
    @property
    def collected_cell_types(self) -> list[str]:
        """
        This function returns a list of cell types
        Index 0 matches Index 0 of ftp_links
        """
        return self._raw_file_cell_types


if __name__ == "__main__":
    print("Use the proteomics_preprocess.py file")
