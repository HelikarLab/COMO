import os
import sys

from src import GSEpipelineFast, project


def test_download_gsm_id_maps():
    config = project.Configs()
    data_dir = config.datadir
    gse = ""
    gpls = []
    vendor = ""
