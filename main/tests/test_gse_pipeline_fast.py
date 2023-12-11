import os
import sys

from src import GSEpipelineFast, project

configs = project.configs


def test_download_gsm_id_maps():
    config = project.configs
    data_dir = config.data_dir
    gse = ""
    gpls = []
    vendor = ""
