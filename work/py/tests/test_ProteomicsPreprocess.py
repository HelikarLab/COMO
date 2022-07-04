import os
import pytest
import sys

# Add parent directory to path, allows us to import the "project.py" file from the parent directory
# From: https://stackoverflow.com/a/30536516/13885200
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from proteomics.proteomics_preprocess import ParseCSVInput
from proteomics.proteomics_preprocess import PopulateInformation
from fixtures_proteomics_preprocess import CSVConstants
from fixtures_proteomics_preprocess
from fixtures_proteomics_preprocess import create_fasta_database




# ----------------------------- #
# Test the ParseCSVInput class  #
# ----------------------------- #
@pytest.fixture
def csv_data() -> CSVConstants:
    """
    Create a fixture that returns data from the above class
    """
    return CSVConstants()
def test_parse_ftp_urls(csv_data):
    parsed_csv = ParseCSVInput(csv_data.file_path)
    assert csv_data.ftp_url == parsed_csv.ftp_urls
def test_parse_cell_types(csv_data):
    parsed_csv = ParseCSVInput(csv_data.file_path)
    assert csv_data.cell_type == parsed_csv.input_cell_types
def test_parse_studies(csv_data):
    parsed_csv = ParseCSVInput(csv_data.file_path)
    assert csv_data.study == parsed_csv.studies


# ---------------------------------- #
# Test the PopulateInformation class #
# ---------------------------------- #
def test_populate_information():
    pass


# ---------------------------------- #
# Test the parse_args() class        #
# ---------------------------------- #
def test_parse_args(csv_data, create_fasta_database, create_argparse_object):
    input_file = csv_data
    parsed_args = create_argparse_object
    
    assert str(parsed_args.input_csv) == str(input_file)
    assert str(parsed_args.database) == str(create_fasta_database)
    assert parsed_args.extensions == ["raw", "fasta"]
    assert parsed_args.skip_download == False
    assert parsed_args.skip_mzml == False
    assert parsed_args.skip_sqt == False
