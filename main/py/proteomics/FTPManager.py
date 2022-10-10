"""
This file is responsible downloading data found at FTP links

TODO: Find a way to mark a file as "downloaded"
    - Keep a list of file names in a ".completd" hidden folder?
"""
import ftplib
from ftplib import FTP
import multiprocessing
from multiprocessing.sharedctypes import Synchronized
import numpy as np
import time
import typing
from urllib.parse import urlparse

from .FileInformation import FileInformation
from .FileInformation import clear_print


def ftp_client(host: str, max_attempts: int = 3, port: int = 21, user: str = "anonymous", passwd: str = "guest") -> FTP:
    """
    This class is responsible for creating a "client" connection
    """
    connection_successful: bool = False
    attempt_num: int = 1

    # Attempt to connect, throw error if unable to do so
    while not connection_successful and attempt_num <= max_attempts:
        try:
            client: FTP = FTP()
            client.connect(host=host, port=port)
            client.login(user=user, passwd=passwd)
            connection_successful = True

        except (ConnectionResetError, ftplib.error_temp) as error:

            if attempt_num > max_attempts:
                print("")
                raise error

            # Make sure the error print statement is on a new line on the first error
            if attempt_num == 1:
                print("")
            
            # Line clean: https://stackoverflow.com/a/5419488/13885200
            clear_print(f"Attempt {attempt_num} of {max_attempts} failed to connect, waiting 5 seconds before trying again")
            attempt_num += 1
            time.sleep(5)
    
    return client


class Reader:
    def __init__(self, root_link: str, file_extensions: list[str], max_attempts: int = 3, port: int = 21, user: str = "anonymous", passwd: str = "guest"):
        """
        This class is responsible for reading data about root FTP links
        """
        self._root_link: str = root_link
        self._extensions: list[str] = file_extensions
        self._max_attempts: int = max_attempts
        self._port: int = port
        self._user: str = user
        self._passwd: str = passwd

        self._files: list[str] = []
        self._file_sizes: list[int] = []
        
        self._get_info()

    # Define a magic function to yeild a tuple containing files and file_sizes
    def __iter__(self):
        yield from zip(self._files, self._file_sizes)
        
    def _get_info(self):
        """
        This function is responsible for getting all files under the root_link
        """
        scheme: str = urlparse(self._root_link).scheme
        host = urlparse(self._root_link).hostname
        folder = urlparse(self._root_link).path

        with ftp_client(
                host=host,
                max_attempts=self._max_attempts,
                port=self._port,
                user=self._user,
                passwd=self._passwd
        ) as client:
            for file_path in client.nlst(folder):
                if file_path[0] != "/":
                    file_path = "/" + file_path

                if file_path.endswith(tuple(self._extensions)):
                    download_url: str = f"{scheme}://{host}{file_path}"
                    self._files.append(download_url)
                    try:
                        self._file_sizes.append(client.size(file_path))
                    except ftplib.error_perm:
                        self._file_sizes.append(0)

    @property
    def files(self) -> typing.Iterator[str]:
        yield self._files

    @property
    def file_names(self) -> typing.Iterator[str]:
        yield from [file.split("/")[-1] for file in self._files]

    @property
    def file_sizes(self) -> typing.Iterable[str]:
        yield self._file_sizes


class Download:
    def __init__(
        self,
        file_information: list[FileInformation],
        core_count: int = 1,
    ):
        """
        This function is responsible for downloading items from the FTP urls passed into it
        """
        self._file_information: list[FileInformation] = file_information
        self._core_count: int = min(core_count, 2)  # Limit to 2 downloads at a time
        self._download_counter: Synchronized = multiprocessing.Value("i", 1)

        # Find files to download
        self.download_data_wrapper()

    def download_data_wrapper(self):
        """
        This function is responsible for using multiprocessing to download as many files at once in parallel
        """
        # Calculate the number of cores ot use
        # We are going to use half the number of cores, OR the number of files to download, whichever is less
        # Otherwise if user specified the number of cores, use that
        print("Starting file download")

        # Split list into chunks to process separately
        file_chunks: list[FileInformation] = np.array_split(self._file_information, self._core_count)

        # Create a list of jobs
        jobs: list[multiprocessing.Process] = []
        
        for i, information in enumerate(file_chunks):

            # Append a job to the list
            job = multiprocessing.Process(
                target=self.download_data,
                args=(information,),
            )
            jobs.append(job)

        [job.start() for job in jobs]       # Start jobs
        [job.join() for job in jobs]        # Wait for jobs to finish
        [job.terminate() for job in jobs]   # Terminate jobs

    def download_data(self, file_information: list[FileInformation]):

        # Start processing the files
        for i, information in enumerate(file_information):
        
            # Define FTP-related variables
            parsed_url = urlparse(information.download_url)
            host = parsed_url.hostname
            path = parsed_url.path

            # Convert file_size from byte to MB
            size_mb: int = round(information.file_size / (1024**2))

            # Get the lock and print file info
            self._download_counter.acquire()
            print(f"Started download {self._download_counter.value:02d} / {len(self._file_information):02d} ({size_mb} MB) - {information.raw_file_name}")  # fmt: skip
            self._download_counter.value += 1
            self._download_counter.release()

            # Download the file
            information.raw_file_path.parent.mkdir(parents=True, exist_ok=True)
            with ftp_client(host=host) as client:
                client.retrbinary(f"RETR {path}", open(information.raw_file_path, "wb").write)


if __name__ == "__main__":
    print("Use the proteomics_preprocess.py file")
