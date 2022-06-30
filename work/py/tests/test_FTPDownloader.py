import os
import sys
from urllib.parse import urlparse

# Add parent directory to path, allows us to import the "project.py" file from the parent directory
# From: https://stackoverflow.com/a/30536516/13885200
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import proteomics.FTPManager as FTPManager


class TestDownloader:
    """
    Test if we are able to find the files we are looking for unde
    """
    root_url: str = "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2022/02/PXD026172"
    host = urlparse(root_url).hostname
    client = FTPManager.ftp_client(host=host)
