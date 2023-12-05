import os
import toml
import subprocess
from pathlib import Path


# find project root dir
class Configs:
    def __init__(self, projectdir: str) -> None:
        self.root_dir = projectdir
        
        self.data_dir = os.path.join(self.root_dir, "data")
        self.config_dir = os.path.join(self.root_dir, "data", "config_sheets")
        self.results_dir = os.path.join(self.data_dir, "results")
        self.src_dir = os.path.join(self.root_dir, "src")
        
        self._configuration_file = os.path.join(self.root_dir, "configuration.toml")
        self.configuration = self._read_configuration()
    
    def _read_configuration(self) -> dict:
        configuration = toml.load(self._configuration_file)
        
        como_version = configuration["general"]["COMO_version"]
        configuration["general"]["COMO_version"] = self._get_como_version(como_version)
        
        return configuration
    
    def update_configuration(
        self,
        write_directory: str | Path,
        config_filename: str = "configuration.toml"
    ) -> None:
        with open(os.path.join(write_directory, config_filename), "w") as f:
            toml.dump(self.configuration, f)
    
    def _get_como_version(self, current_version: str):
        if current_version == "":
            try:
                git_tag_command = "git describe --tags --abbrev=0"
                current_version = subprocess.check_output(git_tag_command, shell=True).decode("utf-8").strip()
            
            except subprocess.CalledProcessError as e:
                raise ValueError("Could not find latest tag") from e
        
        return current_version


current_dir = os.getcwd()
directory_list = current_dir.split("/")

# Find the "main" directory
split_index = 1
for directory in directory_list:
    if directory == "main":
        break  # Exit the loop when we find the "main" directory
    split_index += 1  # Otherwise increment the index

# Unpack items in directory_list
# From: https://stackoverflow.com/questions/14826888
work_dir = os.path.join(*directory_list[0:split_index])

# Add leading "/", as it will not exist right now
work_dir = os.path.join("/", work_dir)
configs: Configs = Configs(work_dir)

if __name__ == '__main__':
    config = Configs("/Users/joshl/PycharmProjects/COMO/main")
    print(config.configuration)
