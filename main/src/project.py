import os


# find project root dir
class Configs:
    def __init__(self, projectdir: str) -> None:
        self.rootdir = projectdir
        self.datadir = os.path.join(projectdir, "data")
        self.configdir = os.path.join(projectdir, "data", "config_sheets")
        self.outputdir = os.path.join(projectdir, "output")
        self.pydir = os.path.join(projectdir, "src")


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
