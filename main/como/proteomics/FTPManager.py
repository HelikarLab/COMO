# ruff: noqa

"""This file is responsible downloading data found at FTP links

TODO: Find a way to mark a file as "downloaded"
    - Keep a list of file names in a ".completd" hidden folder?
"""

import asyncio
import multiprocessing
import time
import typing
from multiprocessing.sharedctypes import Synchronized
from urllib.parse import urlparse

import aioftp

from .FileInformation import FileInformation, clear_print


async def aioftp_client(
    host: str, username: str = "anonymous", password: str = "guest", port: int = 21, max_attempts: int = 3
) -> aioftp.Client:
    """This class is responsible for creating a "client" connection"""
    connection_successful: bool = False
    attempt_num: int = 1

    # Attempt to connect, throw error if unable to do so
    while not connection_successful and attempt_num <= max_attempts:
        try:
            client: aioftp.Client = aioftp.Client()
            await client.connect(host, port)
            await client.login(user=username, password=password)
            connection_successful = True
        except ConnectionResetError:
            # Make sure this print statement is on a new line on the first error
            if attempt_num == 1:
                print()

            # Line clean: https://stackoverflow.com/a/5419488/13885200
            clear_print(f"Attempt {attempt_num} of {max_attempts} failed to connect")
            attempt_num += 1
            time.sleep(5)
    if not connection_successful:
        print()
        raise ConnectionResetError("Could not connect to FTP server")

    return client


class Reader:
    def __init__(
        self,
        root_link: str,
        file_extensions: list[str],
        max_attempts: int = 3,
        port: int = 21,
        user: str = "anonymous",
        passwd: str = "guest",
    ) -> None:
        """This class is responsible for reading data about root FTP links"""
        self._root_link: str = root_link
        self._extensions: list[str] = file_extensions
        self._max_attempts: int = max_attempts
        self._port: int = port
        self._user: str = user
        self._passwd: str = passwd

        self._files: list[str] = []
        self._file_sizes: list[int] = []

        self._get_info_wrapper()

    def _get_info_wrapper(self) -> None:
        event_loop = asyncio.new_event_loop()
        asyncio.set_event_loop(event_loop)
        async_tasks = [self._get_info()]

        event_loop.run_until_complete(asyncio.wait(async_tasks))
        event_loop.close()

    async def _get_info(self) -> None:
        """This function is responsible for getting all files under the root_link"""
        url_parse = urlparse(self._root_link)

        scheme: str
        host: str
        folder: str
        if url_parse.scheme != "":
            scheme = url_parse.scheme
        else:
            scheme = ""
        if url_parse.hostname is not None:
            host = url_parse.hostname
        else:
            raise ValueError(f"Unable to identify hostname from url: {self._root_link}")
        if url_parse.path != "":
            folder = url_parse.path
        else:
            raise ValueError(f"Unable to identify folder or path from url: {self._root_link}")

        client = await aioftp_client(host=host)
        for path, info in await client.list(folder, recursive=True):
            if str(path).endswith(tuple(self._extensions)):
                download_url: str = f"{scheme}://{host}{path}"
                self._files.append(download_url)
                self._file_sizes.append(int(info["size"]))

    @property
    def files(self) -> typing.Iterator[str]:
        for file in self._files:
            yield file
        return self._files

    @property
    def file_names(self) -> typing.Iterator[str]:
        for file in self._files:
            yield file.split("/")[-1]
        # yield from [file.split("/")[-1] for file in self._files]

    @property
    def file_sizes(self) -> typing.Iterable[int]:
        for file in self._file_sizes:
            yield file
        # yield self._file_sizes


class Download:
    def __init__(
        self,
        file_information: list[FileInformation],
        core_count: int = 1,
    ) -> None:
        """This function is responsible for downloading items from the FTP urls passed into it"""
        self._file_information: list[FileInformation] = file_information
        self._core_count: int = min(core_count, 2)  # Limit to 2 downloads at a time
        self._download_counter: Synchronized = Synchronized(multiprocessing.Value("i", 1))
        self._semaphore = asyncio.Semaphore(self._core_count)

        # Find files to download
        self._download_data_wrapper()

    def _download_data_wrapper(self) -> None:
        """This function is responsible for "kicking off" asynchronous data downloading"""
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
        url_parse = urlparse(file_information.download_url)

        scheme: str
        host: str
        folder: str
        if url_parse.scheme != "":
            scheme = url_parse.scheme
        else:
            scheme = ""
        if url_parse.hostname is not None:
            host = url_parse.hostname
        else:
            raise ValueError(f"Unable to identify hostname from url: {file_information.download_url}")
        if url_parse.path != "":
            folder = url_parse.path
        else:
            raise ValueError(f"Unable to identify folder or path from url: {file_information.download_url}")

        # Convert file size from byte to MB
        size_mb: int = round(file_information.file_size / (1024**2))

        # Use a semaphore so only N number of tasks can be started at once
        async with semaphore:
            client = await aioftp_client(host)
            self._download_counter.acquire()
            clear_print(
                f"Started download {self._download_counter.value:02d} / {len(self._file_information):02d} ({size_mb} MB) - {file_information.raw_file_name}"
            )
            self._download_counter.value += 1
            self._download_counter.release()

            # Download file, use "write_into" to write to a file, not a directory
            await client.download(source=folder, destination=file_information.raw_file_path, write_into=True)

        await client.quit()


if __name__ == "__main__":
    print("Use the proteomics_preprocess.py file")
