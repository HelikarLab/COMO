from __future__ import annotations

from pathlib import Path

import pandas as pd
import pytest
from como.rnaseq_preprocess import (
    _organize_gene_counts_files,
    _process_first_multirun_sample,
    _process_standard_replicate,
    _sample_name_from_filepath,
    _STARinformation,
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


class TestSTARInformation:
    valid_data = Path("main/data/COMO_input/naiveB/geneCounts/S1/naiveB_S1R1.tab").resolve()
    invalid_data = Path("main/data/COMO_input/naiveB/fragmentSizes/S1/naiveB_S1R1_fragment_size.txt").resolve()

    @pytest.mark.asyncio
    async def test_build_from_tab_valid_file(self):
        """Validate building STAR information object."""
        star: _STARinformation = await _STARinformation.build_from_tab(TestSTARInformation.valid_data)
        assert len(star.gene_names) == len(star.count_matrix) == 61541
        assert len(star.num_unmapped) == 3
        assert len(star.num_multimapping) == 3
        assert len(star.num_no_feature) == 3
        assert len(star.num_ambiguous) == 3

    @pytest.mark.asyncio
    async def test_build_from_tab_invalid_file(self):
        """Validate error on invalid file."""
        with pytest.raises(ValueError, match="Building STAR information requires a '.tab' file"):
            await _STARinformation.build_from_tab(TestSTARInformation.invalid_data)


def test_sample_name_from_filepath(any_como_input_filepath: Path):
    expected = "_".join(any_como_input_filepath.stem.split("_")[:2])
    assert _sample_name_from_filepath(any_como_input_filepath) == expected


def test_organize_gene_counts_files(como_input_data_directory: Path):
    metrics: list[_StudyMetrics] = _organize_gene_counts_files(como_input_data_directory)
    for metric in metrics:
        assert len(metric.sample_names) == metric.num_samples == len(metric.count_files) == len(metric.strand_files)

        for file in metric.count_files:
            assert f"/{metric.study_name}/" in file.as_posix()
            assert "geneCounts" in file.as_posix()
            assert file.suffix == ".tab"

        for file in metric.strand_files:
            assert f"/{metric.study_name}/" in file.as_posix()
            assert "strandedness" in file.as_posix()
            assert file.suffix == ".txt"


@pytest.mark.asyncio
async def test_process_first_multirun_sample(strand_filepath: Path, all_gene_count_filepaths: list[Path]):
    result: pd.DataFrame = await _process_first_multirun_sample(strand_filepath, all_gene_count_filepaths)
    assert result.columns[0] == "ensembl_gene_id"
    assert len(result.columns) == 2
    assert result.columns[1] in strand_filepath.as_posix()


def test_pack_filepaths(packed_filepaths: PackedFilepaths):
    assert packed_filepaths.sample_name in packed_filepaths.fragment_size.as_posix()
    assert packed_filepaths.sample_name in packed_filepaths.gene_count.as_posix()
    assert packed_filepaths.sample_name in packed_filepaths.insert_size.as_posix()
    assert packed_filepaths.sample_name in packed_filepaths.layout.as_posix()
    assert packed_filepaths.sample_name in packed_filepaths.preparation_method.as_posix()
    assert packed_filepaths.sample_name in packed_filepaths.strandedness.as_posix()
