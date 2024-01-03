import os
import sys

# Add parent directory to path, allows us to import the "project.py" file from the parent directory
# From: https://stackoverflow.com/a/30536516/13885200

from src import project


def test_config():
    configs = project.configs
    work_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    
    assert configs.rootdir == work_dir
    assert configs.datadir == os.path.join(work_dir, "data")
    assert configs.configdir == os.path.join(work_dir, "data", "config_sheets")
    assert configs.outputdir == os.path.join(work_dir, "output")
    assert configs.pydir == os.path.join(work_dir, "src")
