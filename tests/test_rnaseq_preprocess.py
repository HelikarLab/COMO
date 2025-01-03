from __future__ import annotations

import sys

import pytest
from como import rnaseq_preprocess


# Define a list of arguments to test
@pytest.mark.skip("Not testing args because the input 'argv' has been removed")
@pytest.mark.parametrize(
    "args",
    [
        # Test using data in COMO_input data
        ["--context-names", "naiveB immNK", "--gene-format", "Ensembl", "--taxon-id", "9606", "--create-matrix"],
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
    """Asserts that the arguments passed into the function are correct."""
    context_names = args[1]
    gene_format = args[3]
    taxon_id = args[5]
    matrix_mode = args[6]

    sys.argv = args
    parsed = rnaseq_preprocess._parse_args()

    assert [context_name in parsed.context_names for context_name in context_names.split()]
    assert parsed.gene_format == gene_format
    assert parsed.taxon_id == taxon_id

    if matrix_mode == "--create-matrix":
        assert parsed.make_matrix is True
    elif matrix_mode == "--provide-matrix":
        assert parsed.make_matrix is False
        assert parsed.provided_matrix_fname == args[8]
