"""
This file is responsible for viewing and downloading information found at an FTP link
"""
import aioftp
import asyncio
import os
from pathlib import Path
import sys
from urllib.parse import urlparse
from urllib.parse import ParseResult

from pprint import pprint

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
        self._extensions = preferred_extensions

        # If no preferred extensions are passed in, we must convert it to a string before using it
        # Without doing so, str().endswith() returns False. We want it to return True.
        if not self._extensions:
            self._extensions: str = ""
        else:
            # If extensions exist, convert it to a tuple of strings
            self._extensions: tuple[str] = tuple(self._extensions)

        asyncio.run(self.get_ftp_links_wrapper())

    async def get_ftp_links_wrapper(self):
        """
        This function is responsible for asynchronously downloading FTP files from the self._urls list
        """
        tasks = []
        for url in self._ftp_links:
            # Get URL components
            # From: https://stackoverflow.com/questions/9626535
            parsed_url: ParseResult = urlparse(url)

            # We are using "action" as a way to not repeat the asyncio.create_task call
            if self._print_only:
                action = self.list_ftp_files
            else:
                action = self.download_ftp_data
            tasks.append(asyncio.create_task(action(host=parsed_url.netloc, folder=parsed_url.path)))  # fmt: skip

        # Execute the tasks
        await asyncio.wait(tasks)

    async def download_ftp_data(
        self,
        host: str,
        folder: str,
        port: int = 21,
        username: str = "anonymous",
        password: str = "guest",
    ):
        """
        This function is responsible for asynchronously downloading all files under the given path

        :param host: The host to connect to
        :param folder: The path to search under
        :param port: The port to connect to
        :param username: The username to use during login
        :param password: The password to use during login
        """
        # Two clients are required to download multiple files at once
        # From: https://aioftp.readthedocs.io/client_tutorial.html#listing-paths
        list_client: aioftp.Client = aioftp.Client()
        download_client: aioftp.Client = aioftp.Client()

        # Connect "listing" client
        await list_client.connect(host=host, port=port)
        await list_client.login(user=username, password=password)
        await list_client.get_passive_connection()

        # Connect "downloading" client
        await download_client.connect(host=host, port=port)
        await download_client.login(user=username, password=password)
        await download_client.get_passive_connection()

        total_files: list[str] = [
            path if str(path).endswith(self._extensions) else None
            for path, info in (await list_client.list(folder, recursive=True))
        ]

        counter: int = 1
        async for path, info in list_client.list(folder, recursive=True):
            # Test if we are accessing a "raw" file
            if str(path).endswith(self._extensions):

                # Define download paths
                file_name: str = os.path.basename(path)
                save_path: Path = Path(self._output_dir, file_name)
                print(f"Starting download of {file_name} ({counter} / {total_files})")  # fmt: skip
                await download_client.download(source=path, destination=save_path, write_into=True)  # fmt: skip
                counter += 1

        await list_client.quit()
        await download_client.quit()

    async def list_ftp_files(
        self,
        host: str,
        folder: str,
        port: int = 21,
        username: str = "anonymous",
        password: str = "guest",
    ):
        """
        This function is responsible for printing files to the command line

        :param host: The host to connect to
        :param folder: The path of the file
        :param port: The port to connect to
        :param username: The username to use during login
        :param password: The password to use during login
        """
        client: aioftp.Client = aioftp.Client()
        await client.connect(host=host, port=port)
        await client.login(user=username, password=password)

        async for path, info in client.list(folder, recursive=True):
            # Only looking for files
            if info["type"] == "file":

                curr_path: str = str(path)

                # If no preferred extensions have been passed in OR the current path does end with a preferred extension, print it
                if (self._extensions is None) or (curr_path.endswith(self._extensions)):
                    print(path)


def sum_nums(nums: list[int]) -> int:
    """
    This function is responsible for adding an infinite number of numbers and returning the results. It should accept a single number or a list of numbers
    """
    if isinstance(nums, int):
        return nums
    else:
        return sum(nums)


if __name__ == "__main__":
    print("Use the main_proteomics.py file!")
