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


def clear_print(message: str, end: str = "\033[K\r", flush: bool = True) -> None:
    """
    Pass in your message exactly as you would like it printed, and this function will clear the screen and print it.
    """
    print(message, end=end, flush=flush)


class FileInformation:
    # Create an "all_instances" variable of type list[FileInformation]
    # This allows us to search through ALL instances of every FileInformation with functions declared here
    # From: https://stackoverflow.com/a/17253634
    instances: list["FileInformation"] = []
    
    def __init__(
        self,
        cell_type: str,
        download_url: str | None = None,
        study: int | str | None = None,
        raw_path: Path | None = None,
        intensity_csv: Path | None = None,
        mzml_path: Path | None = None,
        sqt_path: Path | None = None,
        file_size: int | None = None,
    ) -> None:
        # File information
        self.cell_type: str = cell_type
        
        if download_url is not None:
            self.download_url: str = download_url
        
        if file_size is not None:
            self.file_size: int = file_size
        
        # Must check for "None", as we are unable to do study[0] on a None object
        self.study: str = ""
        if isinstance(study, str):
            self.study = study
        elif isinstance(study, int):
            self.study = f"S{study}"
        self.replicate: str = ""
        self.batch: str = f"{study}"

        # Base file save paths
        self.raw_base_path: Path
        if raw_path is None:
            self.raw_base_path = Path(project.configs.datadir, "results", cell_type, "proteomics", "raw")
        else:
            self.raw_base_path = raw_path.parent

        self.mzml_base_path: Path
        if mzml_path is None:
            self.mzml_base_path = Path(project.configs.datadir, "results", cell_type, "proteomics", "mzml")
        else:
            self.mzml_base_path = mzml_path.parent
        
        self.sqt_base_path: Path
        if sqt_path is None:
            self.sqt_base_path = Path(project.configs.datadir, "results", cell_type, "proteomics", "sqt")
        else:
            self.sqt_base_path = sqt_path.parent
            
        self.intensity_csv: Path
        if intensity_csv is None:
            self.intensity_csv = Path(project.configs.datadir, "data_matrices", cell_type, f"protein_abundance_matrix_{cell_type}.csv")
        else:
            self.intensity_csv = intensity_csv
            
        # The following variables have inital values set based only on an S# batch, not a replicate
        # The set_replicate function must be called to set the values for a specific replicate, in which these variables will be reset
        # File names
        self.base_name: str = f"{self.cell_type}_{self.batch}_{Path(self.download_url).stem}"
        self.raw_file_name: str = f"{self.base_name}.raw"
        self.mzml_file_name: str = f"{self.base_name}.mzml"
        self.sqt_file_name: str = f"{self.base_name}.target.sqt"
        
        # Full file paths
        self.raw_file_path: Path = Path(self.raw_base_path, self.raw_file_name)
        self.mzml_file_path: Path = Path(self.mzml_base_path, self.mzml_file_name)
        self.sqt_file_path: Path = Path(self.sqt_base_path, self.sqt_file_name)
        
        # Intensity dataframe
        self.base_columns: list[str] = ["uniprot"]
        self.df_columns: list[str] = self.base_columns + [self.batch]
        self.intensity_df: pd.DataFrame = pd.DataFrame(columns=self.df_columns)
        
        FileInformation.instances.append(self)

    def set_replicate(self, replicate: str | int) -> None:
        """
        This function sets self.replicate, and also resets values that use the "replicate" value before it is used
        """
        # Set the initial replicate value
        if isinstance(replicate, str):
            self.replicate = replicate
        else:
            self.replicate = f"R{replicate}"
        
        # "Reset" additional values
        self.batch = f"{self.study}{self.replicate}"
        # File names
        self.base_name = f"{self.cell_type}_{self.batch}_{Path(self.download_url).stem}"
        self.raw_file_name = f"{self.base_name}.raw"
        self.mzml_file_name = f"{self.base_name}.mzml"
        self.sqt_file_name = f"{self.base_name}.target.sqt"

        # Full file paths
        self.raw_file_path = Path(self.raw_base_path, self.raw_file_name)
        self.mzml_file_path = Path(self.mzml_base_path, self.mzml_file_name)
        self.sqt_file_path = Path(self.sqt_base_path, self.sqt_file_name)

        # Intensity dataframe
        self.base_columns = ["uniprot"]
        self.df_columns = self.base_columns + [self.batch]
        self.intensity_df = pd.DataFrame(columns=self.df_columns)

    @classmethod
    def filter_instances(cls, cell_type: str) -> list["FileInformation"]:
        """
        This function finds all FileInformation objects that have the given cell type
        """
        sorted_instances: list[FileInformation] = sorted(cls.instances, key=lambda x: x.study)
        return [instance for instance in sorted_instances if instance.cell_type == cell_type]

    @staticmethod
    def intensity_file_path(cell_type: str) -> Path:
        """
        This function creates a single instance of the FileInformation class and returns the intensity_csv file location
        This is useful because each cell type has a specific location all data gets written to
        If all unique cell types are known, it is then possible to get their intensity csv file location
        """
        information: FileInformation = FileInformation(cell_type=cell_type)
        return information.intensity_csv
