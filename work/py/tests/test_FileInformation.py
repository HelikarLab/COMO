import os
import sys

# Add parent directory to path, allows us to import the "project.py" file from the parent directory
# From: https://stackoverflow.com/a/30536516/13885200
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from proteomics.FileInformation import FileInformation


class TestFileInformation:
    """
    cell_type: str,
    download_url: str = None,
    study: int | str = None,
    raw_path: Path = None,
    intensity_csv: Path = None,
    mzml_path: Path = None,
    sqt_path: Path = None,
    file_size: int = None,
    """
    information = FileInformation(
        cell_type="naiveB",
        download_url="ftp://ftp.pride.ebi.ac.uk/pride-archive/2022/02/PXD026172",
    )
