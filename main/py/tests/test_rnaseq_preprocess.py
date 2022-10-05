import pytest
import os
import sys

# Add parent directory to path, allows us to import the "project.py" file from the parent directory
# From: https://stackoverflow.com/a/30536516/13885200
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import rnaseq_preprocess


# Define a list of arguments to test
@pytest.mark.parametrize(
    "args",
    [
        # Test using data in MADRID_input data
        [
            "--context-names", "naiveB,immNK",
            "--gene-format", "Ensembl",
            "--taxon-id", 9606,
            "--create-matrix"  # Create a matrix with provided data
        ],
    ]
)
def test_arg_input(args):
    """
    This function asserts that the arguments passed into the function are correct
    """
    print(args)
    # --context-names -> naiveB,immNK
    # --gene-format -> Ensembl
    # --taxon-id
    # --provide-matrix
    # --create-matrix
    # --matrix

    parsed_args = rnaseq_preprocess.parse_args(args)
    print(parsed_args)
    assert 0

