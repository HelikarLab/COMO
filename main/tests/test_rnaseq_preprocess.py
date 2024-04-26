import os
import sys

import pytest

# Add parent directory to path, allows us to import the "project.py" file from the parent directory
# From: https://stackoverflow.com/a/30536516/13885200
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import rnaseq_preprocess


# Define a list of arguments to test
@pytest.mark.parametrize(
    "args",
    [
        # Test using data in COMO_input data
        [
            "--context-names",
            "naiveB immNK",
            "--gene-format",
            "Ensembl",
            "--taxon-id",
            "9606",
            "--create-matrix",
        ],
        [
            "--context-names",
            "dimNK brightNK",
            "--gene-format",
            "SYMBOL",
            "--taxon-id",
            "human",
            "--provide-matrix",
            "--matrix",
            "COMO_input/counts_matrix.tsv",
        ],
    ],
)
def test_arg_input(args: list[str]):
    """
    This function asserts that the arguments passed into the function are correct
    """
    print(args)
    context_names = args[1]
    gene_format = args[3]
    taxon_id = args[5]
    matrix_mode = args[6]

    parsed = rnaseq_preprocess.parse_args(args)

    assert [
        context_name in parsed.context_names for context_name in context_names.split()
    ]
    assert parsed.gene_format == gene_format
    assert parsed.taxon_id == taxon_id

    if matrix_mode == "--create-matrix":
        assert parsed.make_matrix is True
    elif matrix_mode == "--provide-matrix":
        assert parsed.make_matrix is False
        assert parsed.provided_matrix_fname == args[8]
