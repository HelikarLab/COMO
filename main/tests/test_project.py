import os

from src import project

# Add parent directory to path, allows us to import the "project.py" file from the parent directory
# From: https://stackoverflow.com/a/30536516/13885200



def test_config():
    configs = project.Configs()
    work_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    
    assert configs.root_dir == work_dir
    assert configs.data_dir == os.path.join(work_dir, "data")
    assert configs.config_dir == os.path.join(work_dir, "data", "config_sheets")
    assert configs.results_dir == os.path.join(work_dir, "data", "results")
    assert configs.src_dir == os.path.join(work_dir, "src")
