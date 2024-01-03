import os
from pathlib import Path


# find project root dir
class Configs:
    project_dir: str = None
    
    def __init__(self) -> None:
        if Configs.project_dir is None:
            Configs.project_dir = self._find_project_dir()
        self.rootdir: str = str(Configs.project_dir)
        
        self.datadir = os.path.join(self.rootdir, "data")
        self.configdir = os.path.join(self.rootdir, "data", "config_sheets")
        self.outputdir = os.path.join(self.rootdir, "output")
        self.pydir = os.path.join(self.rootdir, "src")
    
    def _find_project_dir(self) -> Path:
        current_dir: Path = Path.cwd()
        directory_list: tuple[str, ...] = current_dir.parts
        
        split_index: int = 0
        for directory in directory_list:
            if directory == "main":
                split_index += 1
                break
            split_index += 1
        work_dir: Path = Path(*directory_list[0:split_index])
        return work_dir


# _current_dir = os.getcwd()
# _directory_list: list[str] = _current_dir.split("/")
#
# # Find the "main" directory
# _split_index = 0
# for directory in _directory_list:
#     _split_index += 1  # Increment the index
#     if directory == "main":
#         _split_index += 1
#         break  # Exit the loop when we find the "main" directory
#
# # Unpack items in directory_list
# # From: https://stackoverflow.com/questions/14826888
# _work_dir: str = os.path.join(*_directory_list[0:_split_index])
# print(f"Providing `project.py` with work_dir: {_work_dir}")
# # Add leading "/", as it will not exist right now
# _work_dir = os.path.join("/", _work_dir)
# # configs: _Configs = _Configs(_work_dir)

if __name__ == '__main__':
    configs = Configs()
    print(configs.__dict__)
