"""
This file is responsible for collecting information about the URL passed into it

This includes the information as depicted in the "Database.py" file
"""
import aioftp
import argparse
import asyncio
from enum import Enum
import requests
import os
from pathlib import Path
import sys
import typing
from urllib.parse import urlparse
from urllib.parse import ParseResult

from pprint import pprint

# Add parent directory to path, allows us to import the "project.py" file from the parent directory
# From: https://stackoverflow.com/a/30536516/13885200
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import project


class Action(Enum):
    """
    This class is responsible for holding the types of action the aioftp module is able to perform
    """

    download = "download"
    LIST = "list"


class Download:
    def __init__(self, urls: list[str], args: argparse.Namespace):
        """
        This function is responsible for downloading items from the FTP urls passed into it
        """
        self._output_directory: str = args.output_dir
        # self._urls: list[str] = urls
        self._urls: list[str] = [
            "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2022/02/PXD026142",
            "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2022/03/PXD009340",
            "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2022/03/PXD010319",
        ]

        # This object will be modified by get_ftp_links()
        self._files_to_download: dict[str, list[str]] = {
            "raw_file_url": [],
            "protein_groups_url": [],
        }

        asyncio.run(self.get_ftp_links_wrapper())

    async def get_ftp_links_wrapper(self, **kwargs):
        """
        This function is responsible for asynchronously downloading FTP files from the self._urls list

        :param action: The 'action' to perform. This should be a function name, such as list_ftp_files
        """
        tasks = []
        print("Starting task creation")
        for url in self._urls:
            # Get URL components
            # From: https://stackoverflow.com/questions/9626535
            parsed_url: ParseResult = urlparse(url)

            # Append tasks to complete to the "tasks" list
            tasks.append(
                asyncio.create_task(
                    self.download_ftp_data(host=parsed_url.netloc, path=parsed_url.path, **kwargs)  # fmt: skip
                )
            )

        # Execute the tasks
        print(f"{len(tasks)} task(s) created, starting await. . .")
        await asyncio.wait(tasks)

    async def download_ftp_data(
        self,
        host: str,
        path: str,
        port: int = 21,
        username: str = "anonymous",
        password: str = "guest",
        list_only: bool = False,
    ):
        client = aioftp.Client()
        await client.connect(host=host, port=port)
        await client.login(user=username, password=password)

        async for path, info in client.list(path, recursive=True):
            if list_only:
                if info["type"] == "file":
                    print(path)

            # Test if "raw" or "proteinGroups.txt" in the file path; we want these files
            # From: https://stackoverflow.com/a/8122096/13885200
            if (
                any(map(str(path).__contains__, ["raw", "proteinGroups.txt"]))
                and not list_only
            ):
                # Define paths
                download_path: str = f"ftp://{host}{path}"
                file_name = os.path.basename(path)
                save_path: Path = Path(self._output_directory, file_name)

                # Download file
                try:
                    await client.download(source=download_path, destination=save_path)
                except IndexError:
                    print(f"Index error on file {download_path}")

    async def list_ftp_files(
        self,
        host: str,
        path: str,
        port: int = 21,
        username: str = "anonymous",
        password: str = "guest",
    ):
        client = aioftp.Client()
        await client.connect(host=host, port=port)
        await client.login(user=username, password=password)

        async for path, info in client.list(path, recursive=True):
            if info["type"] == "file":
                print(path)


if __name__ == "__main__":
    print("Use the main_proteomics.py file!")
