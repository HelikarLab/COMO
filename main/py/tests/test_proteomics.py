"""
This file tests the proteomics module
"""

import os
import sys
from pathlib import Path

# Add parent directory to path, allows us to import the "project.py" file from the parent directory
# From: https://stackoverflow.com/a/30536516/13885200
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from proteomics.Crux import MZMLtoSQT, RAWtoMZML, SQTtoCSV
from proteomics.FileInformation import FileInformation
from proteomics.FTPManager import Download, Reader
from proteomics.proteomics_preprocess import ParseCSVInput, PopulateInformation


class TestCrux:
    pass


class TestFileInformation:
    def test_creating_instances(self, tmp_path):
        cell_type: str = "cd8NaiveT"
        download_url: str = "ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2015/02/PXD001065/Aus-141_2.raw"
        study: str = "S1"
        replicate: str = "R1"
        raw_path: Path = tmp_path / f"cd8NaiveT_{study}{replicate}_Aus-141_2.raw"
        mzml_path: Path = tmp_path / f"cd8NaiveT_{study}{replicate}_Aus-141_2.mzML"
        sqt_path: Path = tmp_path / f"cd8NaiveT_{study}{replicate}_Aus-141_2.sqt"
        intensity_path: Path = tmp_path / f"protein_abundance_matrix_{cell_type}.csv"

        instance: FileInformation = FileInformation(
            cell_type=cell_type,
            download_url=download_url,
            study=study,
            raw_path=raw_path,
            intensity_csv=intensity_path,
            mzml_path=mzml_path,
            sqt_path=sqt_path,
            file_size=0
        )
        instance.set_replicate(replicate)

        assert instance.cell_type == cell_type
        assert instance.download_url == download_url
        assert instance.study == study
        assert instance.replicate == replicate
        assert instance.batch == study + replicate
        assert instance.raw_base_path == raw_path.parent
        assert instance.mzml_base_path == mzml_path.parent
        assert instance.sqt_base_path == sqt_path.parent
        assert instance.intensity_csv == intensity_path
        assert instance.file_size == 0


class TestFTPManager:
    pass


class TestProteomicsPreprocess:
    pass
