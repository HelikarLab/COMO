from pathlib import Path

import pytest


@pytest.fixture(scope="module")
def ftp_file_names() -> list[str]:
    return ["file_one.raw", "file_two.raw", "file_three.raw", "config.txt", "README.md"]


@pytest.fixture(scope="module")
def fixture_ftp_server(ftpserver, ftp_file_names):
    # We must "upload" files to the local FTP server
    # This will not upload real files, it will only create files that we can collect with the Reader class
    # From: https://github.com/oz123/pytest-localftpserver/blob/c1ebb0068f4197fad5e14e9c1ef3488de5cf877f/tests/test_pytest_localftpserver.py#L362

    ftpserver.reset_tmp_dirs()

    # Get the base file path to "upload" files
    base_path = Path(ftpserver.get_local_base_path())

    # Define a list of files to create on the local ftp server
    files_on_server: list[Path] = []
    for filename in ftp_file_names:
        abs_file_path: Path = base_path / filename
        files_on_server.append(abs_file_path)
        with abs_file_path.open("w") as f:
            f.write(filename)
