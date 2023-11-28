import os
import sys

from src import GSEpipelineFast, project

configs = project.configs


def test_download_gsm_id_maps():
    config = project.configs
    data_dir = config.datadir
    gse = ""
    gpls = []
    vendor = ""
