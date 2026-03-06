from __future__ import annotations

import asyncio
import csv
import functools
import io
import json
import re
import sys
from dataclasses import dataclass, field
from itertools import chain
from pathlib import Path
from typing import Final, Literal, TextIO, cast

import numpy as np
import numpy.typing as npt
import pandas as pd
from loguru import logger

from como.data_types import LogLevel, RNAType
from como.pipelines.identifier import build_gene_info, get_remaining_identifiers
from como.utils import read_file, set_up_logging


@dataclass
class _QuantInformation:
    gene_names: list[str]
    count_matrix: pd.DataFrame
    sample_name: str
    filepath: Path

    @classmethod
    def build_from_sf(cls, filepath: Path) -> _QuantInformation:
        if filepath.suffix != ".sf":
            raise ValueError(f"Building quantification information requires a '.sf' file; received: '{filepath}'")
        if not filepath.exists():
            raise FileNotFoundError(f"Unable to find the .sf file: {filepath}")

        sample_name = filepath.stem.removesuffix("_quant.genes")
        df = read_file(
            filepath,
            sep="\t",
            names=["ensembl_gene_id", "length", "effective_length", "tpm", sample_name],
            header=0,
        )
        return cls(
            gene_names=df["ensembl_gene_id"].to_list(),
            count_matrix=df,
            sample_name=sample_name,
            filepath=filepath,
        )


@dataclass
class _StudyMetrics:
    study_name: str
    quant_files: list[Path]
    strand_files: list[Path]
    __sample_names: list[str] = field(default_factory=list)
    __num_samples: int = 0

    @property
    def sample_names(self) -> list[str]:
        return self.__sample_names

    @property
    def num_samples(self):
        return self.__num_samples

    def __post_init__(self):
        self.__num_samples = len(self.quant_files)
        self.__sample_names = [f.stem for f in self.quant_files]

        if len(self.quant_files) != len(self.strand_files):
            raise ValueError(
                f"Unequal number of count files and strand files for study '{self.study_name}'. "
                f"Found {len(self.quant_files)} count files and {len(self.strand_files)} strand files."
            )

        if self.num_samples != len(self.quant_files):
            raise ValueError(
                f"Unequal number of samples and count files for study '{self.study_name}'. "
                f"Found {self.num_samples} samples and {len(self.quant_files)} count files."
            )

        if self.num_samples != len(self.strand_files):
            raise ValueError(
                f"Unequal number of samples and strand files for study '{self.study_name}'. "
                f"Found {self.num_samples} samples and {len(self.strand_files)} strand files."
            )

        if self.__num_samples == 1:
            raise ValueError(f"Only one sample exists for study {self.study_name}. Provide at least two samples")

        self.quant_files.sort()
        self.strand_files.sort()
        self.__sample_names.sort()


@dataclass(slots=True)
class SampleConfiguration:
    sample_name: str
    effective_lengths: pd.DataFrame
    mean_effective_length: float
    layout: str
    strand: str
    study: str
    library_prep: str

    def __post_init__(self):
        """Validate the effective lengths dataframe to ensure it has the expected structure and content."""
        if len(self.effective_lengths.columns) > 2:
            raise ValueError(
                f"Effective lengths dataframe for sample '{self.sample_name}' has more than 2 columns, "
                f"expected 'name' and 'effective_length'"
            )

        if "name" not in self.effective_lengths.columns:
            raise ValueError(f"Effective lengths dataframe for sample '{self.sample_name}' is missing 'name' column")

        if "effective_length" not in self.effective_lengths.columns:
            raise ValueError(f"Sample '{self.sample_name}' is missing 'effective_length' column")

    @classmethod
    def to_dataframe(cls, samples: list[SampleConfiguration]) -> tuple[pd.DataFrame, pd.DataFrame]:
        """Convert a list of SampleConfiguration to a dataframe.

        :param samples: The list of SampleConfiguration objects to convert.
        :return: A tuple of dataframes:
            [0]: The sample configuration as a dataframe
            [1]: The effective lengths as a separate data frame with `same_name` as columns
        """
        config = pd.DataFrame(
            columns=["sample_name", "mean_effective_length", "layout", "strand", "study", "library_prep"]
        )

        genes = set()
        for s in samples:
            genes.update(s.effective_lengths["name"].to_list())

        lengths = pd.DataFrame(data=np.float64(0.0), columns=[s.sample_name for s in samples], index=list(genes))
        for sample in samples:
            ids: list[str] = sample.effective_lengths["name"].to_list()
            data: npt.NDArray[np.floating] = sample.effective_lengths["effective_length"].to_numpy(dtype=np.float64)
            lengths.loc[ids, sample.sample_name] = data

        return config, lengths

        return config, lengths


def _sample_name_from_filepath(file: Path) -> str:
    group = re.search(r".+_S\d+R\d+(r\d+)?", file.stem)
    if not group:
        raise ValueError(
            "Filename does not match expected pattern 'contextName_SXRYrZ' where "
            "X is the study number, Y is the replicate number, and Z is the optional run number"
        )
    return group.group()


def _require_one(
    paths: list[Path | None],
    kind: Literal["layout", "strand", "preparation", "fragment"],
    label: str,
) -> Path:
    if len(paths) == 1 and isinstance(paths[0], Path):
        return paths[0]
    if len(paths) > 1:
        message = (
            f"Multiple matching {kind} files for {label}, "
            f"make sure there is only one copy for each replicate in COMO_input"
        )
    elif len(paths) > 1 and paths[0] is None:
        message = f"No {kind} file found for {label}, make sure there is one copy for each replicate in COMO_input"
    else:
        message = f"No {kind} file found for {label}, make sure there is one copy for each replicate in COMO_input"

    raise ValueError(message)


def _organize_gene_counts_files(data_dir: Path) -> list[_StudyMetrics]:
    quant_dir = Path(data_dir, "quantification")
    strand_dir = Path(data_dir, "strandedness")

    if not quant_dir.exists():
        raise FileNotFoundError(f"Quantification directory not found: {quant_dir}")

    if not strand_dir.exists():
        raise FileNotFoundError(f"Strandedness directory not found: {strand_dir}")

    quantification_directories: list[Path] = sorted([p for p in quant_dir.glob("*") if not p.name.startswith(".")])
    strandedness_directories: list[Path] = sorted([p for p in strand_dir.glob("*") if not p.name.startswith(".")])

    if len(quantification_directories) != len(strandedness_directories):
        raise ValueError(
            f"Unequal number of quantification directories and strandedness directories. "
            f"Found {len(quantification_directories)} quantification directories and "
            f"{len(strandedness_directories)} strandedness directories."
            f"\nQuantification directory: {quant_dir}\nStrandedness directory: {strand_dir}"
        )

    # For each study, collect gene count files, fragment files, insert size files, layouts, and strandedness information
    study_metrics: list[_StudyMetrics] = []
    for quant, strand_dir in zip(quantification_directories, strandedness_directories, strict=True):
        quant_files = list(quant.glob("*_quant.genes.sf"))
        strand_files = list(strand_dir.glob("*.txt"))
        if len(quant_files) == 0:
            raise ValueError(f"No quant found for study '{quant.stem}'.")
        if len(strand_files) == 0:
            raise ValueError(f"No strandedness files found for study '{quant.stem}'.")

        study_metrics.append(
            _StudyMetrics(
                study_name=quant.stem,
                quant_files=quant_files,
                strand_files=strand_files,
            )
        )
    return study_metrics


def _process_first_multirun_sample(strand_file: Path, all_quant_files: list[Path]):
    quant_information: list[_QuantInformation] = [_QuantInformation.build_from_sf(f) for f in all_quant_files]

    counts: list[pd.DataFrame] = []
    for info in quant_information:
        count = info.count_matrix[["ensembl_gene_id", info.sample_name]]
        count.columns = ["ensembl_gene_id", "counts"]
        counts.append(count)
    sample_counts = pd.concat(counts, axis=0, ignore_index=True)

    # Set na values to 0
    sample_counts = sample_counts.fillna(value=0)
    sample_counts["counts"] = sample_counts["counts"].astype(float)

    count_avg = cast(  # type checkers think `.groupby(...).mean()` returns a pd.Series, force pd.DataFrame
        pd.DataFrame, cast(object, sample_counts.groupby("ensembl_gene_id", as_index=False)["counts"].mean())
    )
    count_avg["counts"] = np.ceil(count_avg["counts"].astype(int))
    count_avg.columns = ["ensembl_gene_id", _sample_name_from_filepath(strand_file)]
    return count_avg


def _process_standard_replicate(counts_file: Path, strand_file: Path, sample_name: str):
    quant_information = _QuantInformation.build_from_sf(counts_file)
    return quant_information.count_matrix


def _prepare_sample_counts(
    sample_name: str,
    counts_file: Path,
    strand_file: Path,
    all_quant_files: list[Path],
) -> pd.DataFrame | None:
    # Test if the counts_file is the first run in a multi-run smaple
    if re.search(r"R\d+r1", counts_file.as_posix()):
        return _process_first_multirun_sample(strand_file=strand_file, all_quant_files=all_quant_files)
    elif re.search(r"R\d+r[2-9]+", counts_file.as_posix()):
        return None
    else:
        return _process_standard_replicate(counts_file, strand_file, sample_name)


def _create_sample_counts_matrix(metrics: _StudyMetrics) -> pd.DataFrame:
    adjusted_index = 0
    counts: pd.DataFrame | None = _prepare_sample_counts(
        sample_name=metrics.sample_names[0],
        counts_file=metrics.quant_files[0],
        strand_file=metrics.strand_files[0],
        all_quant_files=metrics.quant_files,
    )

    for i in range(1, metrics.num_samples):
        new_counts = _prepare_sample_counts(
            sample_name=metrics.sample_names[i],
            counts_file=metrics.quant_files[i],
            strand_file=metrics.strand_files[i],
            all_quant_files=metrics.quant_files,
        )
        if new_counts is None:
            adjusted_index += 1
            continue

        assert isinstance(counts, pd.DataFrame)  # noqa: S101
        counts = counts.merge(new_counts, on="ensembl_gene_id", how="outer")
        counts = counts.fillna(value=0)

        # Remove run number "r\d+" from multi-run names
        if re.search(r"R\d+r1", metrics.sample_names[i]):
            new_sample_name = re.sub(r"r\d+", "", metrics.sample_names[i])
            old_col_name = counts.columns[i + 1 - adjusted_index]
            counts.rename(columns={old_col_name: new_sample_name}, inplace=True)

    if counts is None:
        raise ValueError(f"No valid counts were processed for study '{metrics.study_name}'")

    return counts


def _write_counts_matrix(
    *,
    config_df: pd.DataFrame,
    fragment_lengths: pd.DataFrame,
    como_context_dir: Path,
    output_counts_matrix_filepath: Path,
    output_fragment_lengths_filepath: Path,
    rna: RNAType,
) -> pd.DataFrame:
    """Create a counts matrix file by reading gene counts table(s).

    :param config_df: Configuration DataFrame containing sample information.
    :param fragment_lengths: DataFrame containing effective lengths for each gene and sample,
        used for zFPKM normalization.
    :param como_context_dir: Path to the COMO_input directory containing gene count files.
    :param output_counts_matrix_filepath: Path where the output counts matrix CSV will be saved.
    :param output_fragment_lengths_filepath: Path where the output fragment lengths CSV will be saved.
    :param rna: RNAType enum indicating whether to process 'trna' or 'mrna' samples.
    :return: A pandas DataFrame representing the final counts matrix.
    """
    study_metrics = _organize_gene_counts_files(data_dir=como_context_dir)
    counts: list[pd.DataFrame] = [_create_sample_counts_matrix(metric) for metric in study_metrics]
    rna_specific_sample_names = set(
        config_df.loc[config_df["library_prep"].str.lower() == rna.value.lower(), "sample_name"].tolist()
    )

    final_matrix: pd.DataFrame = functools.reduce(
        lambda left, right: pd.merge(left, right, on="ensembl_gene_id", how="outer"), counts
    )
    final_matrix.fillna(value=0, inplace=True)
    final_matrix.iloc[:, 1:] = final_matrix.iloc[:, 1:].astype(int)
    final_matrix = final_matrix[["ensembl_gene_id", *rna_specific_sample_names]]

    output_counts_matrix_filepath.parent.mkdir(parents=True, exist_ok=True)
    output_fragment_lengths_filepath.parent.mkdir(parents=True, exist_ok=True)

    final_matrix.to_csv(output_counts_matrix_filepath, index=False)
    fragment_lengths[rna_specific_sample_names].to_csv(output_fragment_lengths_filepath, index=True)

    logger.success(f"Wrote gene count matrix for '{rna.value}' RNA at '{output_counts_matrix_filepath}'")

    return final_matrix


def _create_config_df(  # noqa: C901
    context_name: str,
    /,
    como_context_dir: Path,
    gene_count_dirname: str = "geneCounts",
    layout_dirname: str = "layouts",
    strandedness_dirname: str = "strandedness",
    quantification_dir: str = "quantification",
    prep_method_dirname: str = "prepMethods",
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Create configuration sheet.

    The configuration file created is based on the gene counts matrix.
    If using zFPKM normalization technique, mean fragment lengths will be fetched

    Args:
        context_name: Name of the context, used as a prefix for sample names.
        como_context_dir: Path to the COMO_input directory containing subdirectories for
            gene counts, layouts, strandedness, fragment sizes, and prep methods.
        gene_count_dirname: Name of the subdirectory containing gene count files.
        layout_dirname: Name of the subdirectory containing layout files.
        strandedness_dirname: Name of the subdirectory containing strandedness files.
        quantification_dir: Name of the subdirectory containing Salmon's quantification files.
        prep_method_dirname: Name of the subdirectory containing library preparation method files.

    Returns:
        [0]: A pandas DataFrame representing the configuration sheet.
        [1]: Fragment lengths for downstream calculations
    """
    label_regex: Final = re.compile(r"(?P<study>S\d{1,3})(?P<rep>R\d{1,3})(?P<run>r\d{1,3})?")
    quant_files: list[Path] = list((como_context_dir / quantification_dir).rglob("*.genes.sf"))
    # gene_counts: list[Path] = list((como_context_dir / gene_count_dirname).rglob("*.tab"))
    if not quant_files:
        raise FileNotFoundError(f"No gene count files found in '{gene_count_dirname}'")

    auxillary_directories = {
        "layout": como_context_dir / layout_dirname,
        "strand": como_context_dir / strandedness_dirname,
        "quantification": como_context_dir / quantification_dir,
        "prep": como_context_dir / prep_method_dirname,
    }
    aux_lookup: dict[str, dict[str, Path]] = {kind: {} for kind in auxillary_directories}
    for kind, root in auxillary_directories.items():
        kind: str
        root: Path
        for p in root.rglob("*"):
            if p.is_file():
                m = label_regex.search(p.stem)
                if m:
                    aux_lookup[kind][m.group(0)] = p
    if "layout" not in aux_lookup:
        raise ValueError

    rows: list[SampleConfiguration] = []
    for quant_file in sorted(quant_files):
        m = label_regex.search(quant_file.as_posix())
        if m is None:
            raise ValueError(f"Filename '{quant_file.name}' does not match contextName_SXRYrZ.tab pattern")

        label = m.group()
        study_number = m["study"]
        rep_number = m["rep"]
        sample_id = f"{context_name}_{study_number}{rep_number}"

        layout_path = _require_one([aux_lookup["layout"].get(label)], "layout", label)
        strand_path = _require_one([aux_lookup["strand"].get(label)], "strand", label)
        prep_path = _require_one([aux_lookup["prep"].get(label)], "preparation", label)

        layout = layout_path.read_text().rstrip()
        strand = strand_path.read_text().rstrip()
        prep = prep_path.read_text().rstrip()
        if prep not in {"total", "mrna"}:
            raise ValueError(f"Prep method must be 'total' or 'mrna' (got '{prep}') for {label}")
        if layout == "":
            raise FileNotFoundError(message=f"No layout file found for '{label}'.")

        quant_paths = [p for p in aux_lookup["quantification"].values() if p.name == f"{sample_id}_quant.genes.sf"]
        if (
            not quant_paths
            and layout in ["paired-end", "", None]
            and prep.lower() in [RNAType.TRNA.value.lower(), RNAType.MRNA.value.lower()]
        ):
            log_and_raise_error(
                message=f"No quantification file found for '{label}'; defaulting to 100 bp (needed for zFPKM).",
                error=FileNotFoundError,
                level=LogLevel.WARNING,
            )
        elif len(quant_paths) == 1 and layout == "single-end":
            effective_len = pd.DataFrame({"Name": [], "EffectiveLength": []})
            mean_effective_len = 0.0  # cannot compute FPKM for single-ended data
        else:
            df = read_file(quant_file)
            df.columns = [c.lower() for c in df.columns]
            df = df.rename(columns={"effectivelength": "effective_length"})

            effective_len = df[["name", "effective_length"]]
            effective_len["effective_length"] = effective_len["effective_length"].astype(np.float64)
            mean_effective_len: float = effective_len["effective_length"].sum() / len(df)

        rows.append(
            SampleConfiguration(
                sample_name=sample_id,
                effective_lengths=effective_len,
                mean_effective_length=mean_effective_len,
                layout=layout,
                strand=strand,
                study=study_number,
                library_prep=prep,
            )
        )

    return SampleConfiguration.to_dataframe(rows)


async def _create_gene_info_file(  # noqa: C901
    *,
    counts_matrix_filepaths: list[Path],
    output_filepath: Path,
    taxon: int,
    cache: bool,
):
    """Create a gene information file context.

    The gene information file will be created by reading each matrix filepath in the provided list
    """

    async def read_ensembl_gene_ids(file: Path) -> list[str]:
        data_ = read_file(file, h5ad_as_df=False)
        if isinstance(data_, pd.DataFrame):
            return data_["ensembl_gene_id"].tolist()
        try:
            conversion = get_remaining_identifiers(ids=data_.var_names.tolist(), taxon=taxon)
        except json.JSONDecodeError as e:
            raise ValueError(f"Got a JSON decode error for file '{counts_matrix_filepaths}' ({e})")

        # Remove NA values from entrez_gene_id dataframe column
        conversion = conversion[~conversion["ensembl_gene_id"].isna()]
        return conversion["ensembl_gene_id"].tolist()

    logger.info(
        "Fetching gene info - this can take up to 5 minutes "
        "depending on the number of genes and your internet connection"
    )

    ensembl_ids: set[str] = set(
        chain.from_iterable(await asyncio.gather(*[read_ensembl_gene_ids(f) for f in counts_matrix_filepaths]))
    )
    gene_data: list[dict[str, str | int | list[str] | list[int] | None]] = await MyGene(cache=cache).query(
        items=list(ensembl_ids),
        taxon=taxon,
        scopes="ensemblgene",
    )

    n = len(gene_data)
    all_gene_symbols: list[str] = ["-"] * n
    all_entrez_ids: list[str | int] = ["-"] * n
    all_ensembl_ids: list[str] = ["-"] * n
    all_sizes: list[int] = [-1] * n

    def _avg_pos(value: int | list[int] | None) -> int:
        if value is None:
            return 0
        if isinstance(value, list):
            return int(sum(value) / len(value)) if value else 0
        return int(value)

    for i, data in enumerate(gene_data):
        data: dict[str, str | int | list[str] | list[int] | None]
        if "genomic_pos.start" not in data:
            log_and_raise_error(
                message="Unexpectedly missing key 'genomic_pos.start'", error=KeyError, level=LogLevel.WARNING
            )
        if "genomic_pos.end" not in data:
            log_and_raise_error(
                message="Unexpectedly missing key 'genomic_pos.end'", error=KeyError, level=LogLevel.WARNING
            )
        if "ensembl.gene" not in data:
            log_and_raise_error(
                message="Unexpectedly missing key 'ensembl.gene'", error=KeyError, level=LogLevel.WARNING
            )

        start = data["genomic_pos.start"]
        end = data["genomic_pos.end"]
        ensembl_id = data["ensembl.gene"]

        if not isinstance(start, int):
            log_and_raise_error(
                message=f"Unexpected type for 'genomic_pos.start': expected int, got {type(start)}",
                error=TypeError,
                level=LogLevel.WARNING,
            )
        if not isinstance(end, int):
            log_and_raise_error(
                message=f"Unexpected type for 'genomic_pos.end': expected int, got {type(start)}",
                error=TypeError,
                level=LogLevel.WARNING,
            )
        if not isinstance(ensembl_id, str):
            log_and_raise_error(
                message=f"Unexpected type for 'ensembl.gene': expected str, got {type(ensembl_id)}",
                error=ValueError,
                level=LogLevel.WARNING,
            )

        size = end - start
        all_ensembl_ids[i] = ",".join(map(str, ensembl_id)) if isinstance(ensembl_id, list) else ensembl_id
        all_gene_symbols[i] = str(data.get("symbol", "-"))
        all_entrez_ids[i] = str(data.get("entrezgene", "-"))
        all_sizes[i] = max(size, -1)  # use `size` otherwise -1

    gene_info: pd.DataFrame = pd.DataFrame(
        {
            "ensembl_gene_id": all_ensembl_ids,
            "gene_symbol": all_gene_symbols,
            "entrez_gene_id": all_entrez_ids,
            "size": all_sizes,
        }
    )

    # remove rows where every gene size value is -1 (not available)
    gene_info = gene_info[~(gene_info == -1).all(axis=1)]

    gene_info["ensembl_gene_id"] = gene_info["ensembl_gene_id"].str.split(",")  # extend lists into multiple rows
    gene_info = gene_info.explode(column=["ensembl_gene_id"])
    # we would set `entrez_gene_id` to int here as well, but not all ensembl ids are mapped to entrez ids,
    #   and as a result, there are still "-" values in the entrez id column that cannot be converted to an integer

    gene_info = gene_info.sort_values(by="ensembl_gene_id")
    output_filepath.parent.mkdir(parents=True, exist_ok=True)
    gene_info.to_csv(output_filepath, index=False)
    logger.success(f"Gene Info file written at '{output_filepath}'")


def _process_como_input(
    context_name: str,
    output_config_filepath: Path,
    como_context_dir: Path,
    output_counts_matrix_filepath: Path,
    output_fragment_lengths_filepath: Path,
    rna: RNAType,
) -> None:
    config_df, fragment_lengths = _create_config_df(context_name, como_context_dir=como_context_dir)

    _write_counts_matrix(
        config_df=config_df,
        fragment_lengths=fragment_lengths,
        como_context_dir=como_context_dir,
        output_counts_matrix_filepath=output_counts_matrix_filepath,
        output_fragment_lengths_filepath=output_fragment_lengths_filepath,
        rna=rna,
    )
    with pd.ExcelWriter(output_config_filepath) as writer:
        subset_config = config_df[config_df["library_prep"].str.lower() == rna.value.lower()]
        subset_config.to_excel(writer, sheet_name=context_name, header=True, index=False)


async def _process(
    context_name: str,
    taxon: int,
    output_gene_info_filepath: Path,
    como_context_dir: Path | None,
    input_matrix_filepath: list[Path] | None,
    output_trna_fragment_lengths_filepath: Path | None,
    output_mrna_fragment_lengths_filepath: Path | None,
    output_trna_config_filepath: Path | None,
    output_mrna_config_filepath: Path | None,
    output_trna_matrix_filepath: Path | None,
    output_mrna_matrix_filepath: Path | None,
    *,
    cache: bool,
    create_gene_info_only: bool,
):
    rna_types: list[tuple[RNAType, Path, Path, Path]] = []
    if output_trna_config_filepath and output_trna_matrix_filepath and output_trna_fragment_lengths_filepath:
        rna_types.append(
            (
                RNAType.TRNA,
                output_trna_config_filepath,
                output_trna_matrix_filepath,
                output_trna_fragment_lengths_filepath,
            )
        )
    if output_mrna_config_filepath and output_mrna_matrix_filepath and output_mrna_fragment_lengths_filepath:
        rna_types.append(
            (
                RNAType.MRNA,
                output_mrna_config_filepath,
                output_mrna_matrix_filepath,
                output_mrna_fragment_lengths_filepath,
            )
        )

    # if provided, iterate through como-input specific directories
    if not create_gene_info_only:
        if como_context_dir is None:
            raise ValueError("como_context_dir must be provided if create_gene_info_only is False")
        if output_trna_fragment_lengths_filepath is None:
            raise ValueError("output_fragment_lengths_filepath must be provided if create_gene_info_only is False")

        for rna, out_config, out_matrix, out_frag_len in rna_types:
            _process_como_input(
                context_name=context_name,
                output_config_filepath=out_config,
                como_context_dir=como_context_dir,
                output_counts_matrix_filepath=out_matrix,
                output_fragment_lengths_filepath=out_frag_len,
                rna=rna,
            )

    # create the gene info filepath based on provided data
    input_files = []
    if input_matrix_filepath:
        input_files.extend(input_matrix_filepath)
    if output_trna_matrix_filepath:
        input_files.append(output_trna_matrix_filepath)
    if output_mrna_matrix_filepath:
        input_files.append(output_mrna_matrix_filepath)

    await _create_gene_info_file(
        counts_matrix_filepaths=input_files,
        output_filepath=output_gene_info_filepath,
        taxon=taxon,
        cache=cache,
    )


async def rnaseq_preprocess(  # noqa: C901
    context_name: str,
    taxon: int,
    output_gene_info_filepath: str | Path,
    como_context_dir: str | Path | None = None,
    input_matrix_filepath: str | Path | list[str] | list[Path] | list[str | Path] | None = None,
    output_trna_fragment_lengths_filepath: str | Path | None = None,
    output_mrna_fragment_lengths_filepath: str | Path | None = None,
    output_trna_metadata_filepath: str | Path | None = None,
    output_mrna_metadata_filepath: str | Path | None = None,
    output_trna_count_matrix_filepath: str | Path | None = None,
    output_mrna_count_matrix_filepath: str | Path | None = None,
    cache: bool = True,
    log_level: LogLevel | str = LogLevel.INFO,
    log_location: str | TextIO = sys.stderr,
    *,
    create_gene_info_only: bool = False,
) -> None:
    """Preprocesses RNA-seq data for downstream analysis.

    Fetches additional gene information from a provided matrix or gene counts,
        or optionally creates this matrix using gene count files obtained using STAR aligner

    :param context_name: The context/cell type being processed
    :param taxon: The NCBI taxonomy ID
    :param output_gene_info_filepath: Path to the output gene information CSV file
    :param output_trna_fragment_lengths_filepath: Path to the output tRNA fragment lengths CSV file
    :param output_mrna_fragment_lengths_filepath: Path to the output mRNA fragment lengths CSV file
    :param output_trna_metadata_filepath: Path to the output tRNA config file (if in "create" mode)
    :param output_mrna_metadata_filepath: Path to the output mRNA config file (if in "create" mode)
    :param output_trna_count_matrix_filepath: The path to write total RNA count matrices
    :param output_mrna_count_matrix_filepath: The path to write messenger RNA count matrices
    :param como_context_dir: If in "create" mode, the input path(s) to the COMO_input directory of the current context
        i.e., the directory containing "fragmentSizes", "geneCounts", "insertSizeMetrics", etc. directories
    :param input_matrix_filepath: If in "provide" mode, the path(s) to the count matrices to be processed~
    :param cache: Should HTTP requests be cached
    :param log_level: The logging level
    :param log_location: The logging location
    :param create_gene_info_only: If True, only create the gene info file and skip general preprocessing steps
    """
    set_up_logging(level=log_level, location=log_location)

    # ruff: disable[ASYNC240]
    if not output_gene_info_filepath:
        raise ValueError("output_gene_info_filepath must be provided")

    output_gene_info_filepath = Path(output_gene_info_filepath).resolve()

    context_dir = None
    if como_context_dir is not None:
        context_dir = Path(como_context_dir).resolve()

    in_matrix = None
    if isinstance(input_matrix_filepath, list):
        in_matrix = [Path(i).resolve() for i in input_matrix_filepath]
    else:
        if isinstance(input_matrix_filepath, (str, Path)):
            in_matrix = [Path(input_matrix_filepath).resolve()]

    if output_trna_metadata_filepath is not None:
        output_trna_metadata_filepath = Path(output_trna_metadata_filepath).resolve()
    if output_mrna_metadata_filepath is not None:
        output_mrna_metadata_filepath = Path(output_mrna_metadata_filepath).resolve()
    if output_trna_count_matrix_filepath is not None:
        output_trna_count_matrix_filepath = Path(output_trna_count_matrix_filepath).resolve()
    if output_mrna_count_matrix_filepath is not None:
        output_mrna_count_matrix_filepath = Path(output_mrna_count_matrix_filepath).resolve()
    if output_trna_fragment_lengths_filepath is not None:
        output_trna_fragment_lengths_filepath = Path(output_trna_fragment_lengths_filepath).resolve()
    if output_mrna_fragment_lengths_filepath is not None:
        output_mrna_fragment_lengths_filepath = Path(output_mrna_fragment_lengths_filepath).resolve()

    await _process(
        context_name=context_name,
        taxon=taxon,
        como_context_dir=context_dir,
        input_matrix_filepath=in_matrix,
        output_gene_info_filepath=output_gene_info_filepath,
        output_trna_config_filepath=output_trna_metadata_filepath,
        output_trna_matrix_filepath=output_trna_count_matrix_filepath,
        output_trna_fragment_lengths_filepath=output_trna_fragment_lengths_filepath,
        output_mrna_config_filepath=output_mrna_metadata_filepath,
        output_mrna_matrix_filepath=output_mrna_count_matrix_filepath,
        output_mrna_fragment_lengths_filepath=output_mrna_fragment_lengths_filepath,
        cache=cache,
        create_gene_info_only=create_gene_info_only,
    )
    # ruff: enable[ASYNC240]
