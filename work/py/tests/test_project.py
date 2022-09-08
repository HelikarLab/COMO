import os
import sys

# Add parent directory to path, allows us to import the "project.py" file from the parent directory
# From: https://stackoverflow.com/a/30536516/13885200
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import project


def test_config():
    configs = project.configs
    madrid_work_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    assert configs.rootdir == madrid_work_dir
    assert configs.datadir == os.path.join(madrid_work_dir, "data")
    assert configs.configdir == os.path.join(madrid_work_dir, "data", "config_sheets")
    assert configs.outputdir == os.path.join(madrid_work_dir, "output")
    assert configs.pydir == os.path.join(madrid_work_dir, "py")
