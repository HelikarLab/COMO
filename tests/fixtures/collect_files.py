from pathlib import Path
from typing import NamedTuple

import pytest
from _pytest.fixtures import SubRequest

_fragment_size_filepaths = list(Path("main/data/COMO_input").absolute().rglob("*fragment_size*.txt"))
_gene_count_filepaths = list(Path("main/data/COMO_input").absolute().rglob("*.tab"))
_insert_size_filepaths = list(Path("main/data/COMO_input").absolute().rglob("*_insert_size.txt"))
_layout_filepaths = list(Path("main/data/COMO_input").absolute().rglob("*_layout.txt"))
_preparation_method_filepaths = list(Path("main/data/COMO_input").absolute().rglob("*_prep_method.txt"))
_strandedness_filepaths = list(Path("main/data/COMO_input").absolute().rglob("*_strandedness.txt"))


class PackedFilepaths(NamedTuple):
    sample_name: str
    fragment_size: Path
    gene_count: Path
    insert_size: Path
    layout: Path
    preparation_method: Path
    strandedness: Path


@pytest.fixture(params=_fragment_size_filepaths)
def fragment_size_filepath(request: SubRequest) -> Path:
    return request.param


@pytest.fixture(params=_gene_count_filepaths)
def gene_count_filepath(request: SubRequest) -> Path:
    return request.param


@pytest.fixture
def all_gene_count_filepaths() -> list[Path]:
    return _gene_count_filepaths


@pytest.fixture(params=_insert_size_filepaths)
def insert_size_filepath(request: SubRequest) -> Path:
    return request.param


@pytest.fixture(params=_layout_filepaths)
def layout_filepath(request: SubRequest) -> Path:
    return request.param


@pytest.fixture(params=_preparation_method_filepaths)
def prep_method_filepath(request: SubRequest) -> Path:
    return request.param


@pytest.fixture(params=_strandedness_filepaths)
def strand_filepath(request: SubRequest) -> Path:
    return request.param


@pytest.fixture(
    params=[
        file
        for filepaths in [
            _fragment_size_filepaths,
            _gene_count_filepaths,
            _insert_size_filepaths,
            _layout_filepaths,
            _preparation_method_filepaths,
            _strandedness_filepaths,
        ]
        for file in filepaths
    ]
)
def any_como_input_filepath(request: SubRequest) -> Path:
    return request.param


@pytest.fixture(params=[Path("main/data/COMO_input/naiveB").absolute(), Path("main/data/COMO_input/smB").absolute()])
def como_input_data_directory(request: SubRequest) -> Path:
    return request.param


@pytest.fixture(params=["naiveB", "smB"])
def packed_filepaths(request: SubRequest) -> PackedFilepaths:
    return PackedFilepaths(
        sample_name=request.param,
        fragment_size=Path(f"main/data/COMO_input/{request.param}/fragmentSizes/{request.param}_fragment_size.txt"),
        gene_count=Path(f"main/data/COMO_input/{request.param}/geneCounts/{request.param}.tab"),
        insert_size=Path(f"main/data/COMO_input/{request.param}/insertSizes/{request.param}_insert_size.txt"),
        layout=Path(f"main/data/COMO_input/{request.param}/layouts/{request.param}_layout.txt"),
        preparation_method=Path(f"main/data/COMO_input/{request.param}/prepMethods/{request.param}_prep_method.txt"),
        strandedness=Path(f"main/data/COMO_input/{request.param}/strandedness/{request.param}_strandedness.txt"),
    )
