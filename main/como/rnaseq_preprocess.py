from __future__ import annotations

import functools
import json
import re
import sys
from dataclasses import dataclass, field
from itertools import chain
from pathlib import Path
from typing import Final, Literal, TextIO

import numpy as np
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
    tpm: pd.DataFrame
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
    def to_dataframe(cls, samples: list[SampleConfiguration]) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """Convert a list of SampleConfiguration to a dataframe.

        :param samples: The list of SampleConfiguration objects to convert.
        :return: A tuple of dataframes:
            [0]: The sample configuration as a dataframe
            [1]: The effective lengths as a separate data frame with `sample_name` as columns
            [2]: The salmon TPM values as a separate data frame with `sample_name` as columns
        """
        config: pd.DataFrame = pd.DataFrame.from_records(
            [(s.sample_name, s.mean_effective_length, s.layout, s.strand, s.study, s.library_prep) for s in samples],
            columns=["sample_name", "mean_effective_length", "layout", "strand", "study", "library_prep"],
        )
        lengths = pd.concat(
            [
                s.effective_lengths.set_index("name")["effective_length"]
                .astype(np.float64, copy=False)
                .rename(s.sample_name)
                for s in samples
            ],
            axis=1,
        ).fillna(0.0)
        lengths.index.name = "ensembl_gene_id"
        tpms = pd.concat(
            [s.tpm.set_index("name")["tpm"].astype(np.float64, copy=False).rename(s.sample_name) for s in samples],
            axis=1,
        )
        tpms.index.name = "ensembl_gene_id"

        return config, lengths, tpms


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

    count_avg = sample_counts.groupby("ensembl_gene_id", as_index=False).mean()
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
    if counts is None:
        raise ValueError(f"Unable to process the first sample for study '{metrics.study_name}'")
    counts = counts.drop(columns=["length", "effective_length", "tpm"])

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
        new_counts = new_counts.drop(columns=["length", "effective_length", "tpm"])
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


def _write_matrices(
    *,
    config_df: pd.DataFrame,
    fragment_lengths: pd.DataFrame,
    tpm_df: pd.DataFrame,
    como_context_dir: Path,
    out_counts_filepath: Path,
    out_frag_length_filepath: Path,
    out_tpm_filepath: Path,
    rna: RNAType,
) -> pd.DataFrame:
    """Create a counts matrix file by reading gene counts table(s).

    :param config_df: Configuration DataFrame containing sample information.
    :param fragment_lengths: DataFrame containing effective lengths for each gene and sample,
        used for zFPKM normalization.
    :param tpm_df: DataFrame containing TPM values for each gene and sample, used for zFPKM normalization.
    :param como_context_dir: Path to the COMO_input directory containing gene count files.
    :param out_counts_filepath: Path where the output counts matrix CSV will be saved.
    :param out_frag_length_filepath: Path where the output fragment lengths CSV will be saved.
    :param out_tpm_filepath: Path where the output TPM values CSV will be saved.
    :param rna: RNAType enum indicating whether to process 'trna' or 'mrna' samples.
    :return: A pandas DataFrame representing the final counts matrix.
    """
    study_metrics = _organize_gene_counts_files(data_dir=como_context_dir)
    counts: list[pd.DataFrame] = [_create_sample_counts_matrix(metric) for metric in study_metrics]
    rna_specific_sample_names = sorted(
        set(config_df.loc[config_df["library_prep"].str.lower() == rna.value.lower(), "sample_name"].tolist())
    )

    final_matrix: pd.DataFrame = functools.reduce(
        lambda left, right: pd.merge(left, right, on="ensembl_gene_id", how="outer"), counts
    )
    final_matrix.fillna(value=0, inplace=True)
    final_matrix = final_matrix[["ensembl_gene_id", *rna_specific_sample_names]]
    final_matrix = final_matrix.reindex(columns=["ensembl_gene_id", *rna_specific_sample_names])
    fragment_lengths = fragment_lengths.reindex(columns=rna_specific_sample_names).sort_index()
    tpm_df = tpm_df.reindex(columns=rna_specific_sample_names).sort_index()

    out_counts_filepath.parent.mkdir(parents=True, exist_ok=True)
    out_frag_length_filepath.parent.mkdir(parents=True, exist_ok=True)
    final_matrix.to_csv(out_counts_filepath, index=False, float_format="%.15g")
    fragment_lengths[rna_specific_sample_names].to_csv(out_frag_length_filepath, index=True, float_format="%.15g")
    tpm_df[rna_specific_sample_names].to_csv(out_tpm_filepath, index=True, float_format="%.15g")

    logger.success(f"Wrote gene count matrix for '{rna.value}' RNA at '{out_counts_filepath}'")

    return final_matrix


def _create_config_df(
    context_name: str,
    /,
    como_context_dir: Path,
    layout_dirname: str = "layouts",
    strandedness_dirname: str = "strandedness",
    quantification_dir: str = "quantification",
    prep_method_dirname: str = "prepMethods",
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Create configuration sheet.

    The configuration file created is based on the gene counts matrix.
    If using zFPKM normalization technique, mean fragment lengths will be fetched

    Args:
        context_name: Name of the context, used as a prefix for sample names.
        como_context_dir: Path to the COMO_input directory containing subdirectories for
            gene counts, layouts, strandedness, fragment sizes, and prep methods.
        layout_dirname: Name of the subdirectory containing layout files.
        strandedness_dirname: Name of the subdirectory containing strandedness files.
        quantification_dir: Name of the subdirectory containing Salmon's quantification files.
        prep_method_dirname: Name of the subdirectory containing library preparation method files.

    Returns:
        [0]: A pandas DataFrame representing the configuration sheet.
        [1]: Fragment lengths for downstream calculations
        [2]: Salmon TPM values for downstream calculations
    """
    label_regex: Final = re.compile(r"(?P<study>S\d{1,3})(?P<rep>R\d{1,3})(?P<run>r\d{1,3})?")
    quant_files: list[Path] = sorted((como_context_dir / quantification_dir).rglob("*.genes.sf"))
    # gene_counts: list[Path] = list((como_context_dir / gene_count_dirname).rglob("*.tab"))
    if not quant_files:
        raise FileNotFoundError(f"No gene count files found in '{quantification_dir}'")

    auxillary_directories = {
        "layout": como_context_dir / layout_dirname,
        "strand": como_context_dir / strandedness_dirname,
        "quantification": como_context_dir / quantification_dir,
        "prep": como_context_dir / prep_method_dirname,
    }
    aux_lookup: dict[str, dict[str, Path]] = {kind: {} for kind in auxillary_directories}
    for kind, root in auxillary_directories.items():
        files = (p for p in root.rglob("*") if p.is_file())
        for file in files:
            m = label_regex.search(file.stem)
            if m:
                aux_lookup[kind][m.group(0)] = file

    rows: list[SampleConfiguration] = []
    for quant_file in quant_files:
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
            raise FileNotFoundError(f"No layout file found for '{label}'.")

        df = read_file(quant_file)
        df.columns = [c.lower() for c in df.columns]
        df = df.rename(columns={"effectivelength": "effective_length", "TPM": "tpm"})
        df["effective_length"] = df["effective_length"].astype(float)
        df["tpm"] = df["tpm"].astype(float)

        tpm = df[["name", "tpm"]].copy()
        effective_len = df[["name", "effective_length"]].copy()
        mean_effective_len: float = effective_len["effective_length"].sum() / df.shape[1]

        rows.append(
            SampleConfiguration(
                sample_name=sample_id,
                effective_lengths=effective_len,
                tpm=tpm,
                mean_effective_length=mean_effective_len,
                layout=layout,
                strand=strand,
                study=study_number,
                library_prep=prep,
            )
        )

    return SampleConfiguration.to_dataframe(rows)


def _create_gene_info_file(
    *,
    counts_matrix_filepaths: list[Path],
    output_filepath: Path,
    taxon: int,
    cache: bool,
):
    """Create a gene information file context.

    The gene information file will be created by reading each matrix filepath in the provided list
    """

    def read_ensembl_gene_ids(file: Path) -> list[str]:
        data_ = read_file(file, h5ad_as_df=False)
        if isinstance(data_, pd.DataFrame):
            return data_["ensembl_gene_id"].tolist()
        try:
            conversion = get_remaining_identifiers(ids=data_.var_names.tolist(), taxon=taxon)
        except json.JSONDecodeError as e:
            raise ValueError(f"Got a JSON decode error for file '{counts_matrix_filepaths}' ({e})") from e

        # Remove NA values from entrez_gene_id dataframe column
        conversion = conversion[~conversion["ensembl_gene_id"].isna()]
        conversion = conversion.explode(column=["ensembl_gene_id"])
        return conversion["ensembl_gene_id"].tolist()

    logger.info(
        "Fetching gene info - this can take up to 5 minutes "
        "depending on the number of genes and your internet connection"
    )

    ensembl_ids: set[str] = set(chain.from_iterable(read_ensembl_gene_ids(f) for f in counts_matrix_filepaths))
    gene_info = build_gene_info(ids=ensembl_ids, taxon=taxon, cache=cache)
    output_filepath.parent.mkdir(parents=True, exist_ok=True)
    gene_info.to_csv(output_filepath, index=False)
    logger.success(f"Gene Info file written at '{output_filepath}'")


def _process_como_input(
    context_name: str,
    como_context_dir: Path,
    out_mrna_config: Path | None,
    out_trna_config: Path | None,
    out_mrna_counts: Path | None,
    out_trna_counts: Path | None,
    out_mrna_fragments: Path | None,
    out_trna_fragments: Path | None,
    out_mrna_tpm: Path | None,
    out_trna_tpm: Path | None,
) -> None:
    config_df, fragment_lengths, tpm_df = _create_config_df(context_name, como_context_dir=como_context_dir)

    if out_trna_config and out_trna_counts and out_trna_fragments and out_trna_tpm:
        _write_matrices(
            config_df=config_df,
            fragment_lengths=fragment_lengths,
            tpm_df=tpm_df,
            como_context_dir=como_context_dir,
            out_counts_filepath=out_trna_counts,
            out_frag_length_filepath=out_trna_fragments,
            out_tpm_filepath=out_trna_tpm,
            rna=RNAType.TRNA,
        )
        with pd.ExcelWriter(out_trna_config) as writer:
            subset_config = config_df[config_df["library_prep"].str.lower() == RNAType.TRNA.value.lower()]
            subset_config.to_excel(writer, sheet_name=context_name, header=True, index=False)

    if out_mrna_config and out_mrna_counts and out_mrna_fragments and out_mrna_tpm:
        _write_matrices(
            config_df=config_df,
            fragment_lengths=fragment_lengths,
            tpm_df=tpm_df,
            como_context_dir=como_context_dir,
            out_counts_filepath=out_mrna_counts,
            out_frag_length_filepath=out_mrna_fragments,
            out_tpm_filepath=out_mrna_tpm,
            rna=RNAType.MRNA,
        )
        with pd.ExcelWriter(out_mrna_config) as writer:
            subset_config = config_df[config_df["library_prep"].str.lower() == RNAType.MRNA.value.lower()]
            subset_config.to_excel(writer, sheet_name=context_name, header=True, index=False)


def _process(
    context_name: str,
    taxon: int,
    output_gene_info_filepath: Path,
    como_context_dir: Path | None,
    in_matrix_fp: list[Path] | None,
    out_trna_frag_lengths_fp: Path | None,
    out_mrna_frag_lengths_fp: Path | None,
    out_trna_config_fp: Path | None,
    out_mrna_config_fp: Path | None,
    out_trna_matrix_fp: Path | None,
    out_mrna_matrix_fp: Path | None,
    out_mrna_tpm_fp: Path | None,
    out_trna_tpm_fp: Path | None,
    *,
    cache: bool,
    create_gene_info_only: bool,
):
    if not create_gene_info_only:
        if como_context_dir is None:
            raise ValueError("como_context_dir must be provided if create_gene_info_only is False")
        if out_trna_frag_lengths_fp is None:
            raise ValueError("output_fragment_lengths_filepath must be provided if create_gene_info_only is False")

        _process_como_input(
            context_name=context_name,
            como_context_dir=como_context_dir,
            out_mrna_config=out_mrna_config_fp,
            out_trna_config=out_trna_config_fp,
            out_mrna_counts=out_mrna_matrix_fp,
            out_trna_counts=out_trna_matrix_fp,
            out_mrna_fragments=out_mrna_frag_lengths_fp,
            out_trna_fragments=out_trna_frag_lengths_fp,
            out_mrna_tpm=out_mrna_tpm_fp,
            out_trna_tpm=out_trna_tpm_fp,
        )

    # create the gene info filepath based on provided data
    input_files = []
    if in_matrix_fp:
        input_files.extend(in_matrix_fp)
    if out_trna_matrix_fp:
        input_files.append(out_trna_matrix_fp)
    if out_mrna_matrix_fp:
        input_files.append(out_mrna_matrix_fp)

    _create_gene_info_file(
        counts_matrix_filepaths=input_files,
        output_filepath=output_gene_info_filepath,
        taxon=taxon,
        cache=cache,
    )


def rnaseq_preprocess(  # noqa: C901
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
    output_trna_tpm_filepath: str | Path | None = None,
    output_mrna_tpm_filepath: str | Path | None = None,
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
    :param output_trna_tpm_filepath: The path to write total RNA TPM matrices
    :param output_mrna_tpm_filepath: The path to write messenger RNA TPM matrices
    :param como_context_dir: If in "create" mode, the input path(s) to the COMO_input directory of the current context
        i.e., the directory containing "fragmentSizes", "geneCounts", "insertSizeMetrics", etc. directories
    :param input_matrix_filepath: If in "provide" mode, the path(s) to the count matrices to be processed~
    :param cache: Should HTTP requests be cached
    :param log_level: The logging level
    :param log_location: The logging location
    :param create_gene_info_only: If True, only create the gene info file and skip general preprocessing steps
    """
    set_up_logging(level=log_level, location=log_location)

    if not output_gene_info_filepath:
        raise ValueError("output_gene_info_filepath must be provided")

    if any(
        (
            output_trna_fragment_lengths_filepath,
            output_trna_metadata_filepath,
            output_trna_count_matrix_filepath,
            output_trna_tpm_filepath,
        ),
    ) and not all(
        (
            output_trna_fragment_lengths_filepath,
            output_trna_metadata_filepath,
            output_trna_count_matrix_filepath,
            output_trna_tpm_filepath,
        ),
    ):
        missing = [
            name
            for name, value in {
                "output_trna_fragment_lengths_filepath": output_trna_fragment_lengths_filepath,
                "output_trna_metadata_filepath": output_trna_metadata_filepath,
                "output_trna_count_matrix_filepath": output_trna_count_matrix_filepath,
                "output_trna_tpm_filepath": output_trna_tpm_filepath,
            }.items()
            if value is None
        ]
        raise ValueError(
            f"If any tRNA output filepaths are provided, all must be provided. Missing: {', '.join(missing)}"
        )

    if any(
        (
            output_mrna_fragment_lengths_filepath,
            output_mrna_metadata_filepath,
            output_mrna_count_matrix_filepath,
            output_mrna_tpm_filepath,
        ),
    ) and not all(
        (
            output_mrna_fragment_lengths_filepath,
            output_mrna_metadata_filepath,
            output_mrna_count_matrix_filepath,
            output_mrna_tpm_filepath,
        ),
    ):
        missing = [
            name
            for name, value in {
                "output_mrna_fragment_lengths_filepath": output_mrna_fragment_lengths_filepath,
                "output_mrna_metadata_filepath": output_mrna_metadata_filepath,
                "output_mrna_count_matrix_filepath": output_mrna_count_matrix_filepath,
                "output_mrna_tpm_filepath": output_mrna_tpm_filepath,
            }.items()
            if value is None
        ]
        raise ValueError(
            f"If any mRNA output filepaths are provided, all must be provided. Missing: {', '.join(missing)}"
        )

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
    if output_trna_tpm_filepath is not None:
        output_trna_tpm_filepath = Path(output_trna_tpm_filepath).resolve()
    if output_mrna_tpm_filepath is not None:
        output_mrna_tpm_filepath = Path(output_mrna_tpm_filepath).resolve()

    _process(
        context_name=context_name,
        taxon=taxon,
        como_context_dir=context_dir,
        in_matrix_fp=in_matrix,
        output_gene_info_filepath=output_gene_info_filepath,
        out_trna_config_fp=output_trna_metadata_filepath,
        out_trna_matrix_fp=output_trna_count_matrix_filepath,
        out_trna_frag_lengths_fp=output_trna_fragment_lengths_filepath,
        out_mrna_config_fp=output_mrna_metadata_filepath,
        out_mrna_matrix_fp=output_mrna_count_matrix_filepath,
        out_mrna_frag_lengths_fp=output_mrna_fragment_lengths_filepath,
        out_trna_tpm_fp=output_trna_tpm_filepath,
        out_mrna_tpm_fp=output_mrna_tpm_filepath,
        cache=cache,
        create_gene_info_only=create_gene_info_only,
    )
