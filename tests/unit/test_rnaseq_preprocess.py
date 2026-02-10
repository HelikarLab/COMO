from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
from como.rnaseq_preprocess import (
    _organize_gene_counts_files,
    _process_first_multirun_sample,
    _QuantInformation,
    _sample_name_from_filepath,
    _StudyMetrics,
)

from tests.fixtures.collect_files import (
    PackedFilepaths,
    all_gene_count_filepaths,
    any_como_input_filepath,
    como_input_data_directory,
    packed_filepaths,
    strand_filepath,
)


class TestQuantInformation:
    valid_data = Path("main/data/COMO_input/naiveB/quantification/S1/naiveB_S1R1_quant.genes.sf").resolve()
    invalid_data = Path("main/data/COMO_input/naiveB/strandedness/S1/naiveB_S1R1_strandedness.txt").resolve()

    def test_build_from_sf_valid_file(self) -> None:
        quant: _QuantInformation = _QuantInformation.build_from_sf(TestQuantInformation.valid_data)
        assert len(quant.gene_names) == len(quant.count_matrix) == 78899
        assert quant.sample_name == "naiveB_S1R1"
        assert quant.filepath.as_posix().endswith(
            "/COMO/main/data/COMO_input/naiveB/quantification/S1/naiveB_S1R1_quant.genes.sf"
        )

    def test_build_from_sf_invalid_file(self):
        with pytest.raises(ValueError, match=r"Building quantification information requires a '.sf' file; received: "):
            _QuantInformation.build_from_sf(TestQuantInformation.invalid_data)

    def test_build_from_missing_file(self):
        with pytest.raises(FileNotFoundError, match=r"Unable to find the .sf file"):
            _QuantInformation.build_from_sf(Path("missing_file.sf"))


def test_sample_name_from_filepath(any_como_input_filepath: Path):
    expected = "_".join(any_como_input_filepath.stem.split("_")[:2])
    assert _sample_name_from_filepath(any_como_input_filepath) == expected


def test_organize_gene_counts_files(como_input_data_directory: Path):
    metric: _StudyMetrics
    for metric in _organize_gene_counts_files(como_input_data_directory):
        assert len(metric.sample_names) == metric.num_samples == len(metric.quant_files) == len(metric.strand_files)

        for file in metric.quant_files:
            assert f"/{metric.study_name}/" in file.as_posix()
            assert file.as_posix().endswith("_quant.genes.sf")
            assert file.suffix == ".sf"

        for file in metric.strand_files:
            assert f"/{metric.study_name}/" in file.as_posix()
            assert "strandedness" in file.as_posix()
            assert file.suffix == ".txt"


def test_process_first_multirun_sample(strand_filepath: Path, all_gene_count_filepaths: list[Path]):
    result: pd.DataFrame = _process_first_multirun_sample(strand_filepath, all_gene_count_filepaths)
    assert result.columns[0] == "ensembl_gene_id"
    assert len(result.columns) == 2
    assert result.columns.tolist()[1] in strand_filepath.as_posix()


def test_pack_filepaths(packed_filepaths: PackedFilepaths):
    assert packed_filepaths.sample_name in packed_filepaths.fragment_size.as_posix()
    assert packed_filepaths.sample_name in packed_filepaths.gene_count.as_posix()
    assert packed_filepaths.sample_name in packed_filepaths.insert_size.as_posix()
    assert packed_filepaths.sample_name in packed_filepaths.layout.as_posix()
    assert packed_filepaths.sample_name in packed_filepaths.preparation_method.as_posix()
    assert packed_filepaths.sample_name in packed_filepaths.strandedness.as_posix()
