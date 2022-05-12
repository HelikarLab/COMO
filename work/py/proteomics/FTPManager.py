"""
This file is responsible for viewing and downloading information found at an FTP link
"""
from ftplib import FTP
import multiprocessing
import numpy as np
import os
from pathlib import Path
import sys
from urllib.parse import urlparse
from urllib.parse import ParseResult

# Add parent directory to path, allows us to import the "project.py" file from the parent directory
# From: https://stackoverflow.com/a/30536516/13885200
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import project


class Download:
    def __init__(
        self,
        ftp_links: list[str],
        output_dir: Path,
        print_only: bool = False,
        preferred_extensions: list[str] = None,
        core_count: int = 1,
    ):
        """
        This function is responsible for downloading items from the FTP urls passed into it

        :param ftp_links: A list of FTP links to download from
        :param output_dir: The directory to save data into
        :param print_only: Should the files only be printed to console? If True, nothing will be downloaded, even if an output directory is specified
        :param preferred_extensions: The file extension to look for. If not set, it will recursively find every file
        """
        self._ftp_links: list[str] = ftp_links
        self._output_dir: Path = output_dir
        self._print_only: bool = print_only
        self._extensions: list[str] = preferred_extensions
        self._core_count: int = core_count

        # Create a manager and lock to print values onto screen
        self._file_lock = multiprocessing.Lock()
        self._download_counter = multiprocessing.Value("i", 1)

        # If no preferred extensions are passed in, we must convert it to a string before using it
        # Without doing so, str().endswith() returns False. We want it to return True.
        if not self._extensions:
            self._extensions: str = ""
        else:
            # If extensions exist, convert it to a tuple of strings
            self._extensions: tuple[str] = tuple(self._extensions)

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
        for url in self._ftp_links:
            parsed_url: ParseResult = urlparse(url)
            scheme: str = parsed_url.scheme
            host: str = parsed_url.hostname
            folder: str = parsed_url.path

            client: FTP = FTP(host=host, user="anonymous", passwd="guest")
            for file_path in client.nlst(folder):
                if file_path.endswith(self._extensions):
                    download_url: str = f"{scheme}://{host}{file_path}"
                    file_size: int = client.size(file_path)

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
                args=(files, sizes, self._file_lock, self._download_counter),
            )
            jobs.append(job)

        [job.start() for job in jobs]  # Start jobs
        [job.join() for job in jobs]  # Wait for jobs to finish
        [job.terminate() for job in jobs]  # Terminate jobs

    def download_data(
        self,
        files_to_download: list[str],
        file_sizes: list[int],
        lock,
        download_counter,
    ):

        # Start processing the files
        for i, (file, size) in enumerate(zip(files_to_download, file_sizes)):
            # Define FTP-related variables
            parsed_url = urlparse(file)
            host = parsed_url.hostname
            path = parsed_url.path

            # Convert file_size from byte to MB
            size_mb: int = size // (1024**2)

            # Determine the file name using os.path.basename
            file_name: str = os.path.basename(path)
            save_path: Path = Path(self._output_dir, file_name)

            # Connect to the host, login, and download the file
            client: FTP = FTP(host=host, user="anonymous", passwd="guest")

            # Get the lock and print file info
            lock.acquire()
            print(f"Started {download_counter.value:02d} / {len(self._files_to_download):02d} ({size_mb} MB) - {file_name}")  # fmt: skip
            download_counter.value += 1
            lock.release()

            # Download the file
            client.retrbinary(f"RETR {path}", open(save_path, "wb").write)
            client.quit()


if __name__ == "__main__":
    print("Use the main_proteomics.py file!")
