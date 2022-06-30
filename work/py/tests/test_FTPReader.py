import os
import sys
from urllib.parse import urlparse
# Add parent directory to path, allows us to import the "project.py" file from the parent directory
# From: https://stackoverflow.com/a/30536516/13885200
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import proteomics.FTPManager as FTPManager


class TestFTPClient:
    root_url: str = "ftp://ftp.pride.ebi.ac.uk/"
    
    def test_ftp_client(self):
        host = urlparse(self.root_url).hostname
        client = FTPManager.ftp_client(host=host)
        assert client is not None

    def test_ftp_reading_folders(self):
        reader = FTPManager.Reader(self.root_url, file_extensions=["pride", "pride-archive"])
        
        # Only looking for "pride" and "pride-archive" folders in the list
        assert f"{self.root_url}pride" in reader.files
        assert f"{self.root_url}pride-archive" in reader.files
    
    def test_ftp_reading_sizes(self):
        reader = FTPManager.Reader(self.root_url, file_extensions=["pride", "pride-archive"])
        for size in reader.file_sizes:
            assert size == 0
