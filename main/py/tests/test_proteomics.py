"""
This file tests the proteomics module
"""

from pathlib import Path

import aioftp
import pytest
from proteomics.Crux import MZMLtoSQT, RAWtoMZML, SQTtoCSV
from proteomics.FileInformation import FileInformation
from proteomics.FTPManager import Download, Reader, aioftp_client
from proteomics.proteomics_preprocess import ParseCSVInput, PopulateInformation

from fixtures.fixture_ftp_server import fixture_ftp_server
from fixtures.fixture_ftp_server import ftp_file_names


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
    @pytest.mark.asyncio
    async def test_ftp_client(self):
        """
        This test checks that the ftp_client function works as expected
        """
        host: str = "ftp.pride.ebi.ac.uk"
        port: int = 21
        client: aioftp.Client = await aioftp_client(host=host, port=port, username="anonymous", password="guest")

        assert client.server_host == host
        assert client.server_port == port
        assert await client.get_current_directory() == Path("/")
        assert await client.quit() is None

    @pytest.mark.skip(reason="pyftpdlib is broken, no way to test this")
    def test_reader(self, ftpserver, fixture_ftp_server, ftp_file_names):
        # Use pytest_localftpserver and fixtures.fixture_ftp_server.fix
        # Now we can get login information for our local FTP server
        file_extensions: list[str] = ["raw"]
        login_data: str = ftpserver.get_login_data(style="url", anon=True)
        dict_login_data: dict = ftpserver.get_login_data(style="dict", anon=False)
        port: int = dict_login_data["port"]
        user: str = dict_login_data["user"]
        passwd: str = dict_login_data["passwd"]

        # Now we can use our reader class and collect files from the local FTP server
        reader = Reader(
            file_extensions=file_extensions,
            root_link=login_data,
            port=port,
            user=user,
            passwd=passwd
        )

        # Assert that we have files ending in "raw", and that the file name found is in the ftp file names list
        for file in reader.file_names:
            assert file.endswith(tuple(file_extensions))
            assert file in ftp_file_names

    @pytest.mark.skip(reason="pyftpdlib is broken, no way to test this")
    def test_download(self, mocker, tmp_path):
        """
        This checks that the Download class works as expected
        """


class TestProteomicsPreprocess:
    pass
