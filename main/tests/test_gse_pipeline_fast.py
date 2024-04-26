import os
import sys

import pytest
from src import GSEpipelineFast, project


@pytest.mark.skip(reason="Skipping microarray tests")
def test_download_gsm_id_maps():
    config = project.Configs()
    data_dir = config.data_dir
    gse = ""
    gpls = []
    vendor = ""
