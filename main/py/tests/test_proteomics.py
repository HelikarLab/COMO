"""
This file tests the proteomics module
"""

import os
import sys
from pathlib import Path
from pytest_mock import mocker

# Add parent directory to path, allows us to import the "project.py" file from the parent directory
# From: https://stackoverflow.com/a/30536516/13885200
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from proteomics.Crux import MZMLtoSQT, RAWtoMZML, SQTtoCSV
from proteomics.FileInformation import FileInformation
from proteomics.FTPManager import Download, Reader, ftp_client
from proteomics.proteomics_preprocess import ParseCSVInput, PopulateInformation


class TestCrux:
    pass


class TestFileInformation:
    def test_creating_instances(self, tmp_path):
        """
        This is an example instance. If anything is incorrectly changed in the FileInformation class, this test will fail.
        """

        # Define several constants
        cell_type: str = "t_gondii"
        study: str = "S1"
        replicate: str = "R1"
        file_name: str = "2DE_2DE_10-2D"

        # Define paths; this file was chosen because it is only 14MB in size
        download_url: str = f"ftp://ftp.pride.ebi.ac.uk/pride/data/archive/2015/02/PXD000297/{file_name}.raw"
        raw_path: Path = Path(tmp_path, f"{cell_type}_{study}{replicate}_{file_name}.raw")
        mzml_path: Path = Path(tmp_path, f"{cell_type}_{study}{replicate}_{file_name}.mzml")
        sqt_path: Path = Path(tmp_path, f"{cell_type}_{study}{replicate}_{file_name}.target.sqt")
        intensity_path: Path = Path(tmp_path, f"protein_abundance_matrix_{cell_type}.csv")

        # Create an instance of FileInformation
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

        # Set the replicate (`R1`)
        instance.set_replicate(replicate)

        # Check that the instance is correct
        assert instance.cell_type == cell_type
        assert instance.download_url == download_url
        assert instance.study == study
        assert instance.replicate == replicate
        assert instance.batch == study + replicate
        assert str(instance.raw_file_path) == str(raw_path)
        assert str(instance.mzml_file_path) == str(mzml_path)
        assert str(instance.sqt_file_path) == str(sqt_path)
        assert str(instance.intensity_csv) == str(intensity_path)
        assert instance.file_size == 0


class TestFTPManager:
    def test_ftp_client(self):
        """
        This test checks that the ftp_client function works as expected
        """
        host: str = "ftp.pride.ebi.ac.uk"
        port: int = 21
        client = ftp_client(host=host, port=port, user="anonymous", passwd="guest")

        assert client.host == host
        assert client.port == port


class TestProteomicsPreprocess:
    pass
