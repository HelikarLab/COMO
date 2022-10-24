import os
import sys

# Add parent directory to path, allows us to import the "project.py" file from the parent directory
# From: https://stackoverflow.com/a/30536516/13885200
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import GSEpipelineFast
from fixtures.fixture_main_dir import main_dir
import project

configs = project.configs
madrid_work_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))


def test_download_gsm_id_maps(main_dir):
    config = project.configs
    data_dir = config.datadir
    gse = ""
    gpls = []
    vendor = ""
