import os
from pathlib import Path

# Add parent directory to path, allows us to import the "project.py" file from the parent directory
# From: https://stackoverflow.com/a/30536516/13885200

from como import project


def test_config():
    configs = project.Config()
    root_dir = Path(__file__).parent.parent

    assert configs.data_dir == root_dir / "data"
    assert configs.config_dir == root_dir / "data" / "config_sheets"
    assert configs.result_dir == root_dir / "data" / "results"
    assert configs.code_dir == root_dir / "como"
