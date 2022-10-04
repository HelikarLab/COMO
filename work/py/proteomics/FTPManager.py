"""
This file is responsible downloading data found at FTP links

TODO: Find a way to mark a file as "downloaded"
    - Keep a list of file names in a ".completd" hidden folder?
"""
from FileInformation import FileInformation
from FileInformation import clear_print

import aioftp
import asyncio
import multiprocessing
from multiprocessing.sharedctypes import Synchronized
import time
from urllib.parse import urlparse


async def aioftp_client(host: str, username: str = "anonymous", password: str = "guest", max_attempts: int = 3) -> aioftp.Client:
    """
    This class is responsible for creating a "client" connection
    """
    connection_successful: bool = False
    attempt_num: int = 1
    # Attempt to connect, throw error if unable to do so
    while not connection_successful and attempt_num <= max_attempts:
        try:
            client: aioftp.Client = aioftp.Client()
            await client.connect(host, 21)
            await client.login(user=username, password=password)
            connection_successful = True
        except ConnectionResetError:

            # Make sure this print statement is on a new line on the first error
            if attempt_num == 1:
                print("")

            # Line clean: https://stackoverflow.com/a/5419488/13885200
            clear_print(f"Attempt {attempt_num} of {max_attempts} failed to connect")
            attempt_num += 1
            time.sleep(5)
    if not connection_successful:
        print("")
        raise ConnectionResetError("Could not connect to FTP server")

    return client


class Reader:
    def __init__(self, root_link: str, file_extensions: list[str]):
        """
        This class is responsible for reading data about root FTP links
        """
        self._root_link: str = root_link
        self._extensions: list[str] = file_extensions
        self._files: list[str] = []
        self._file_sizes: list[int] = []
        
        self._get_info_wrapper()

    def _get_info_wrapper(self):
        event_loop = asyncio.new_event_loop()
        asyncio.set_event_loop(event_loop)
        async_tasks = [self._get_info()]

        event_loop.run_until_complete(asyncio.wait(async_tasks))
        event_loop.close()

    async def _get_info(self):
        """
        This function is responsible for getting all files under the root_link
        """
        scheme: str = urlparse(self._root_link).scheme
        host = urlparse(self._root_link).hostname
        folder = urlparse(self._root_link).path

        client = await aioftp_client(host)
        for path, info in (await client.list(folder, recursive=True)):
            if str(path).endswith(tuple(self._extensions)):
                download_url: str = f"{scheme}://{host}{path}"
                self._files.append(download_url)
                self._file_sizes.append(int(info["size"]))

    @property
    def files(self):
        return self._files
    
    @property
    def file_sizes(self):
        return self._file_sizes


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
        self._semaphore = asyncio.Semaphore(self._core_count)

        # Find files to download
        self._download_data_wrapper()

    def _download_data_wrapper(self):
        """
        This function is responsible for "kicking off" asynchronous data downloading
        """
        print("Starting file download")

        event_loop = asyncio.new_event_loop()
        asyncio.set_event_loop(event_loop)
        async_tasks = []
        for file_information in self._file_information:
            async_tasks.append(
                self._aioftp_download_data(
                    file_information=file_information,
                    semaphore=self._semaphore,
                )
            )

        # Await all the tasks
        event_loop.run_until_complete(asyncio.wait(async_tasks))
        event_loop.close()

    async def _aioftp_download_data(self, file_information: FileInformation, semaphore: asyncio.Semaphore) -> None:

        # Define FTP-related variables
        parsed_url = urlparse(file_information.download_url)
        host = parsed_url.hostname
        path = parsed_url.path

        # Convert file size from byte to MB
        size_mb: int = round(file_information.file_size / (1024 ** 2))

        # Use a semaphore so only N number of tasks can be started at once
        async with semaphore:
            client = await aioftp_client(host)
            self._download_counter.acquire()
            print(f"Started download {self._download_counter.value:02d} / {len(self._file_information):02d} ({size_mb} MB) - {file_information.raw_file_name}")  # fmt: skip
            self._download_counter.value += 1
            self._download_counter.release()

            # Download file, use "write_into" to write to a file, not a directory
            await client.download(source=path, destination=file_information.raw_file_path, write_into=True)

        await client.quit()


if __name__ == "__main__":
    print("Use the proteomics_preprocess.py file")
