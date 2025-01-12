from como.data_types import SourceTypes


def test_source_types():
    """Validate that source types always go in the order of 'trna', 'mrna', 'scrna', 'protemics'."""
    expected_order = ["trna", "mrna", "scrna", "proteomics"]
    for expected, source in zip(expected_order, SourceTypes):
        expected: str
        source: SourceTypes
        assert expected == source.value
