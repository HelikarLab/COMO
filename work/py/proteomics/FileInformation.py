"""
This will hold all relevant information about a single file to download

This should be implemented as a list of objects
"""
import os
import sys
from pathlib import Path

import pandas as pd

# Add parent directory to path, allows us to import the "project.py" file from the parent directory
# From: https://stackoverflow.com/a/30536516/13885200
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import project


def clear_print(message: str, end: str = "\033[K\r", flush: bool = True):
    """
    Pass in your message exactly as you would like it printed, and this function will clear the screen and print it.
    """
    print(message, end=end, flush=flush)


class FileInformation:
    def __init__(
        self,
        cell_type: str,
        download_url: str = None,
        replicate: str = None,
        raw_path: Path = None,
        intensity_csv: Path = None,
        mzml_path: Path = None,
        sqt_path: Path = None,
        file_size: int = None,
    ):
        # File information
        self.cell_type: str = cell_type
        self.download_url: str = download_url
        self.replicate: str = replicate
        self.file_size: int = file_size
        self.prefix: str = f"{cell_type}_{replicate}"
        
        # Base file save paths
        if raw_path is None:
            self.raw_base_path: Path = Path(project.configs.datadir, "results", cell_type, "proteomics", "raw")
        else:
            self.raw_base_path: Path = raw_path.parent
        
        if mzml_path is None:
            self.mzml_base_path: Path = Path(project.configs.datadir, "results", cell_type, "proteomics", "mzml")
        else:
            self.mzml_base_path: Path = mzml_path.parent
        
        if sqt_path is None:
            self.sqt_base_path: Path = Path(project.configs.datadir, "results", cell_type, "proteomics", "sqt")
        else:
            self.sqt_base_path: Path = sqt_path.parent

        # File names
        self.intensity_csv: Path = Path(project.configs.datadir, "data_matrices", cell_type, f"protein_abundance_matrix_{cell_type}.csv",)
        self.base_name: str = f"{cell_type}_{replicate}_{Path(download_url).stem}"
        self.raw_file_name: str = f"{self.base_name}.raw"
        self.mzml_file_name: str = f"{self.base_name}.mzml"
        self.sqt_file_name: str = f"{self.base_name}.target.sqt"
        
        # Full file paths
        self.raw_file_path: Path = Path(self.raw_base_path, self.raw_file_name)
        self.mzml_file_path: Path = Path(self.mzml_base_path, self.mzml_file_name)
        self.sqt_file_path: Path = Path(self.sqt_base_path, self.sqt_file_name)
        
        # Intensity dataframe
        self.base_columns: list[str] = ["uniprot"]
        self.df_columns: list[str] = self.base_columns + [self.prefix]
        self.intensity_df: pd.DataFrame = pd.DataFrame(columns=self.df_columns)

    @staticmethod
    def intensity_file_path(cell_type: str) -> Path:
        """
        This function creates a single instance of the FileInformation class and returns the intensity_csv file location
        This is useful because each cell type has a specific location all data gets written to
        If all unique cell types are known, it is then possible to get their intensity csv file location
        """
        information: FileInformation = FileInformation(
            cell_type=cell_type,
            download_url="",
            replicate="",
            raw_path=Path(""),
            intensity_csv=Path(""),
            mzml_path=Path(""),
            sqt_path=Path(""),
            file_size=0,
        )
        return information.intensity_csv
