import os
from pathlib import Path


class Configs:
    project_dir: str = None
    
    def __init__(self) -> None:
        """
        This class will hold all of the project configurations.
        """
        if Configs.project_dir is None:
            Configs.project_dir = self._find_project_dir()
        self.root_dir: str = str(Configs.project_dir)
        
        self.data_dir: str = os.path.join(self.root_dir, "data")
        self.config_dir: str = os.path.join(self.data_dir, "config_sheets")
        self.results_dir: str = os.path.join(self.data_dir, "results")
        self.src_dir: str = os.path.join(self.root_dir, "src")
    
    def _find_project_dir(self) -> Path:
        # Determine if "main" is in the current directory (i.e., `ls .`)
        if "main" in os.listdir(os.getcwd()):
            return Path(os.getcwd(), "main")
        
        # If "main" is not in the current directory, see if we are in a sub-directory of "main" (i.e., `ls ..`)
        current_directory = Path.cwd()
        directory_parts = current_directory.parts
        for i, directory in enumerate(directory_parts):
            if directory == "main":
                return Path(*directory_parts[0:i + 1])
        
        # If we are not in a sub-directory of "main", then we are not in the project directory
        raise FileNotFoundError(
            "Could not find project directory. "
            "Please ensure you are in the project directory, such as 'COMO' or 'main'."
        )
