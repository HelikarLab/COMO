"""
This file is responsible for collecting information about the URL passed into it

This includes the information as depicted in the "Database.py" file
"""
import aioftp
import argparse
import asyncio
import requests
import os
import sys
from urllib.parse import urlparse
from urllib.parse import ParseResult

from pprint import pprint

# Add parent directory to path, allows us to import the "project.py" file from the parent directory
# From: https://stackoverflow.com/a/30536516/13885200
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import project


class Download:
    def __init__(self, urls: list[str], args: argparse.Namespace):
        """
        This function is responsible for downloading items from the FTP urls passed into it
        """
        # self._urls: list[str] = urls
        self._urls: list[str] = [
            "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2022/02/PXD026142",
            # "ftp://massive.ucsd.edu/MSV000087556/",
            # "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2022/03/PXD029911",
            # "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2022/03/PXD030463",
            # "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2022/03/PXD008681",
            # "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2022/03/PXD029500",
            # "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2022/03/PXD008994",
            # "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2022/03/PXD006541",
            "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2022/03/PXD009340",
            "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2022/03/PXD010319",
        ]

        # self._invalid_urls: list[str] = self.url_validator()
        # if len(self._invalid_urls) != 0:
        #     print(f"The following FTP URLs are invalid. Re-enter them into the '{args.input_file}' file")  # fmt: skip
        #     for url in self._invalid_urls:
        #         print(url)

        asyncio.run(self.download_wrapper())

    def url_validator(self) -> list[str]:
        """
        This function is responsible for checking all FTP URLs passed into the Download() class are valid

        It will return a list of invalid URLs, if any exist
        """
        invalid_url: list[str] = []
        for url in self._urls:
            if not requests.head(url).ok:
                invalid_url.append(url)

        return invalid_url

    async def download_wrapper(self):
        """
        This function is responsible for asynchronously downloading FTP files from the self._urls list
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
                    self.get_ftp(host=parsed_url.netloc, path=parsed_url.path)
                )
            )

        # Execute the tasks
        print(f"{len(tasks)} task(s) created, starting await. . .")
        data = await asyncio.wait(tasks)
        pprint(data)

    async def get_ftp(
        self,
        host: str,
        path: str,
        port: int = 21,
        username: str = "anonymous",
        password: str = "guest",
    ):
        client = aioftp.Client()
        await client.connect(host, port=port)
        await client.login(username, password)

        async for path, info in client.list(path, recursive=True):
            if "proteinGroups.txt" in str(path):
                print(path)


if __name__ == "__main__":
    print("Use the main_proteomics.py file!")
