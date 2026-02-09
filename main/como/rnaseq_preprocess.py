from __future__ import annotations

import asyncio
import csv
import functools
import io
import json
import re
import sys
from collections.abc import Sequence
from dataclasses import asdict, dataclass, field
from itertools import chain
from pathlib import Path
from typing import Final, Literal, cast

import numpy as np
import numpy.typing as npt
import pandas as pd
import pandera.pandas as pa
import pandera.typing.pandas as pat
from fast_bioservices.biothings.mygene import MyGene
from fast_bioservices.pipeline import gene_symbol_to_ensembl_and_gene_id
from loguru import logger

from como.data_types import PATH_TYPE, LogLevel, RNAType
from como.utils import _listify, _log_and_raise_error, _read_file, _set_up_logging


@dataclass
class _QuantInformation:
    gene_names: list[str]
    count_matrix: pd.DataFrame
    sample_name: str
    filepath: Path

    @classmethod
    def build_from_sf(cls, filepath: Path) -> _QuantInformation:
        if filepath.suffix != ".sf":
            _log_and_raise_error(
                f"Building quantification information requires a '.sf' file; received: '{filepath}'",
                error=ValueError,
                level=LogLevel.ERROR,
            )
        if not filepath.exists():
            _log_and_raise_error(
                f"Unable to find the .sf file: {filepath}",
                error=FileNotFoundError,
                level=LogLevel.ERROR,
            )

        sample_name = filepath.stem.removesuffix("_quant.genes.sf")
        df = pd.read_csv(
            io.StringIO(filepath.read_text()),
            sep="\t",
            names=["ensembl_gene_id", "length", "effective_length", "tpm", sample_name],
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
            _log_and_raise_error(
                (
                    f"Unequal number of count files and strand files for study '{self.study_name}'. "
                    f"Found {len(self.quant_files)} count files and {len(self.strand_files)} strand files."
                ),
                error=ValueError,
                level=LogLevel.ERROR,
            )

        if self.num_samples != len(self.quant_files):
            _log_and_raise_error(
                (
                    f"Unequal number of samples and count files for study '{self.study_name}'. "
                    f"Found {self.num_samples} samples and {len(self.quant_files)} count files."
                ),
                error=ValueError,
                level=LogLevel.ERROR,
            )

        if self.num_samples != len(self.strand_files):
            _log_and_raise_error(
                (
                    f"Unequal number of samples and strand files for study '{self.study_name}'. "
                    f"Found {self.num_samples} samples and {len(self.strand_files)} strand files."
                ),
                error=ValueError,
                level=LogLevel.ERROR,
            )

        if self.__num_samples == 1:
            _log_and_raise_error(
                f"Only one sample exists for study {self.study_name}. Provide at least two samples",
                error=ValueError,
                level=LogLevel.ERROR,
            )

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
        if len(self.effective_lengths.columns) > 2:
            _log_and_raise_error(
                message=f"Effective lengths dataframe for sample '{self.sample_name}' has more than 2 columns, expected 'name' and 'effective_length'",
                error=ValueError,
                level=LogLevel.ERROR,
            )
        if "name" not in self.effective_lengths.columns:
            _log_and_raise_error(
                message=f"Effective lengths dataframe for sample '{self.sample_name}' is missing 'name' column",
                error=ValueError,
                level=LogLevel.ERROR,
            )
        if "effective_length" not in self.effective_lengths.columns:
            _log_and_raise_error(
                message=f"Effective lengths dataframe for sample '{self.sample_name}' is missing 'effective_length' column",
                error=ValueError,
                level=LogLevel.ERROR,
            )

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


def _sample_name_from_filepath(file: Path) -> str:
    return re.search(r".+_S\d+R\d+(r\d+)?", file.stem).group()


def _require_one(
    paths: list[Path],
    kind: Literal["layout", "strand", "preparation", "fragment"],
    label: str,
) -> Path | None:
    if len(paths) == 1:
        return paths[0]
    if len(paths) == 0:
        return None
    _log_and_raise_error(
        f"Multiple matching {kind} files for {label}, make sure there is only one copy for each replicate in COMO_input",
        error=ValueError,
        level=LogLevel.ERROR,
    )


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
        _log_and_raise_error(
            (
                f"Unequal number of quantification directories and strandedness directories. "
                f"Found {len(quantification_directories)} quantification directories and "
                f"{len(strandedness_directories)} strandedness directories."
                f"\nQuantification directory: {quant_dir}\nStrandedness directory: {strand_dir}"
            ),
            error=ValueError,
            level=LogLevel.ERROR,
        )

    # For each study, collect gene count files, fragment files, insert size files, layouts, and strandedness information
    study_metrics: list[_StudyMetrics] = []
    for quant, strand_dir in zip(quantification_directories, strandedness_directories, strict=True):
        quant_files = list(quant.glob("*_quant.genes.sf"))
        strand_files = list(strand_dir.glob("*.txt"))
        if len(quant_files) == 0:
            _log_and_raise_error(f"No quant found for study '{quant.stem}'.", error=ValueError, level=LogLevel.ERROR)
        if len(strand_files) == 0:
            _log_and_raise_error(
                f"No strandedness files found for study '{quant.stem}'.",
                error=ValueError,
                level=LogLevel.ERROR,
            )

        study_metrics.append(
            _StudyMetrics(
                study_name=quant.stem,
                quant_files=quant_files,
                strand_files=strand_files,
            )
        )
    return study_metrics


def _process_first_multirun_sample(strand_file: Path, all_quant_files: list[Path]):
    sample_count = pd.DataFrame()
    quant_information: list[_QuantInformation] = [_QuantInformation.build_from_sf(f) for f in all_quant_files]

    for info in quant_information:
        run_counts = info.count_matrix[["ensembl_gene_id", info.sample_name]]
        run_counts.columns = ["ensembl_gene_id", "counts"]
        sample_count = (
            run_counts if sample_count.empty else sample_count.join(run_counts, on=["ensembl_gene_id"], how="outer")
        )

    # Set na values to 0
    sample_count = sample_count.fillna(value=0)
    sample_count["counts"] = sample_count["counts"].astype(float)

    count_avg = sample_count.groupby("ensembl_gene_id", as_index=False)["counts"].mean()
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
        counts: pd.DataFrame = counts.merge(new_counts, on="ensembl_gene_id", how="outer")
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
    :param fragment_lengths: DataFrame containing effective lengths for each gene and sample, used for zFPKM normalization.
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

    final_matrix: pd.DataFrame = functools.reduce(lambda left, right: pd.merge(left, right, on="ensembl_gene_id", how="outer"), counts)
    final_matrix.fillna(value=0, inplace=True)
    final_matrix.iloc[:, 1:] = final_matrix.iloc[:, 1:].astype(int)
    final_matrix = cast(pd.DataFrame, final_matrix[["ensembl_gene_id", *rna_specific_sample_names]])

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
        _log_and_raise_error(
            f"No gene count files found in '{gene_count_dirname}'",
            error=FileNotFoundError,
            level=LogLevel.ERROR,
        )

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

    rows: list[SampleConfiguration] = []
    for quant_file in sorted(quant_files):
        m = label_regex.search(quant_file.as_posix())
        if m is None:
            _log_and_raise_error(
                f"Filename '{quant_file.name}' does not match contextName_SXRYrZ.tab pattern",
                error=ValueError,
                level=LogLevel.ERROR,
            )
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
            _log_and_raise_error(
                f"Prep method must be 'total' or 'mrna' (got '{prep}') for {label}",
                error=ValueError,
                level=LogLevel.ERROR,
            )
        if layout == "":
            _log_and_raise_error(
                message=f"No layout file found for '{label}'.",
                error=FileNotFoundError,
                level=LogLevel.WARNING,
            )

        quant_paths = [p for p in aux_lookup["quantification"].values() if p.name == f"{sample_id}_quant.genes.sf"]
        if (
            not quant_paths
            and layout in ["paired-end", "", None]
            and prep.lower() in [RNAType.TRNA.value.lower(), RNAType.MRNA.value.lower()]
        ):
            _log_and_raise_error(
                message=f"No quantification file found for '{label}'; defaulting to 100 bp (needed for zFPKM).",
                error=FileNotFoundError,
                level=LogLevel.WARNING,
            )
        elif len(quant_paths) == 1 and layout == "single-end":
            effective_len = pd.DataFrame({"Name": [], "EffectiveLength": []})
            mean_effective_len = 0.0  # cannot compute FPKM for single-ended data
        else:
            df = _read_file(quant_file)
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

    # 6-3-25: Intentionally left commented-out code to test its replacement
    # gene_counts_dir = como_context_dir / gene_count_dirname
    # layout_dir = como_context_dir / layout_dirname
    # strandedness_dir = como_context_dir / strandedness_dirname
    # fragment_sizes_dir = como_context_dir / fragment_sizes_dirname
    # prep_method_dir = como_context_dir / prep_method_dirname
    #
    # gene_counts_files = list(gene_counts_dir.rglob("*.tab"))
    # sample_names: list[str] = []
    # fragment_lengths: list[int | float] = []
    # layouts: list[str] = []
    # strands: list[str] = []
    # groups: list[str] = []
    # preparation_method: list[str] = []
    #
    # if len(gene_counts_files) == 0:
    #     _log_and_raise_error(f"No gene count files found in '{gene_counts_dir}'.", error=FileNotFoundError, level=LogLevel.ERROR)
    #
    # for gene_count_filename in sorted(gene_counts_files):
    #     # Match S___R___r___
    #     # \d{1,3} matches 1-3 digits
    #     # (?:r\d{1,3})? optionally matches a "r" followed by three digits
    #     label = re.findall(r"S\d{1,3}R\d{1,3}(?:r\d{1,3})?", gene_count_filename.as_posix())[0]
    #     if not label:
    #         _log_and_raise_error(
    #             (
    #                 f"\n\nFilename of '{gene_count_filename}' is not valid. "
    #                 f"Should be 'contextName_SXRYrZ.tab', "
    #                 f"where X is the study/batch number, Y is the replicate number, "
    #                 f"and Z is the run number."
    #                 "\n\nIf not a multi-run sample, exclude 'rZ' from the filename."
    #             ),
    #             error=ValueError,
    #             level=LogLevel.ERROR,
    #         )
    #
    #     study_number = re.findall(r"S\d{1,3}", label)[0]
    #     rep_number = re.findall(r"R\d{1,3}", label)[0]
    #     run_number = re.findall(r"r\d{1,3}", label)
    #
    #     multi_flag = 0
    #     if len(run_number) > 0:
    #         if run_number[0] != "r1":
    #             continue
    #         label_glob = f"{study_number}{rep_number}r*"  # S__R__r*
    #         runs = [run for run in gene_counts_files if re.search(label_glob, run.as_posix())]
    #         multi_flag = 1
    #         frag_files = []
    #
    #         for run in runs:
    #             run_number = re.findall(r"R\d{1,3}", run.as_posix())[0]
    #             replicate = re.findall(r"r\d{1,3}", run.as_posix())[0]
    #             frag_filename = "".join([context_name, "_", study_number, run_number, replicate, "_fragment_size.txt"])
    #             frag_files.append(como_context_dir / fragment_sizes_dirname / study_number / frag_filename)
    #
    #     layout_files: list[Path] = list(layout_dir.rglob(f"{context_name}_{label}_layout.txt"))
    #     strand_files: list[Path] = list(strandedness_dir.rglob(f"{context_name}_{label}_strandedness.txt"))
    #     frag_files: list[Path] = list(fragment_sizes_dir.rglob(f"{context_name}_{label}_fragment_size.txt"))
    #     prep_files: list[Path] = list(prep_method_dir.rglob(f"{context_name}_{label}_prep_method.txt"))
    #
    #     layout = "UNKNOWN"
    #     if len(layout_files) == 0:
    #         logger.warning(
    #             f"No layout file found for {label}, writing as 'UNKNOWN', "
    #             f"this should be defined if you are using zFPKM or downstream 'rnaseq_gen.py' will not run"
    #         )
    #     elif len(layout_files) == 1:
    #         with layout_files[0].open("r") as file:
    #             layout = file.read().strip()
    #     elif len(layout_files) > 1:
    #         _log_and_raise_error(
    #             f"Multiple matching layout files for {label}, make sure there is only one copy for each replicate in COMO_input",
    #             error=ValueError,
    #             level=LogLevel.ERROR,
    #         )
    #
    #     strand = "UNKNOWN"
    #     if len(strand_files) == 0:
    #         logger.warning(
    #             f"No strandedness file found for {label}, writing as 'UNKNOWN'. "
    #             f"This will not interfere with the analysis since you have already set rnaseq_preprocess.py to "
    #             f"infer the strandedness when writing the counts matrix"
    #         )
    #     elif len(strand_files) == 1:
    #         with strand_files[0].open("r") as file:
    #             strand = file.read().strip()
    #     elif len(strand_files) > 1:
    #         _log_and_raise_error(
    #             f"Multiple matching strandedness files for {label}, make sure there is only one copy for each replicate in COMO_input",
    #             error=ValueError,
    #             level=LogLevel.ERROR,
    #         )
    #
    #     prep = "total"
    #     if len(prep_files) == 0:
    #         logger.warning(f"No prep file found for {label}, assuming 'total', as in 'Total RNA' library preparation")
    #     elif len(prep_files) == 1:
    #         with prep_files[0].open("r") as file:
    #             prep = file.read().strip().lower()
    #             if prep not in ["total", "mrna"]:
    #                 _log_and_raise_error(
    #                     f"Prep method must be either 'total' or 'mrna' for {label}",
    #                     error=ValueError,
    #                     level=LogLevel.ERROR,
    #                 )
    #     elif len(prep_files) > 1:
    #         _log_and_raise_error(
    #             f"Multiple matching prep files for {label}, make sure there is only one copy for each replicate in COMO_input",
    #             error=ValueError,
    #             level=LogLevel.ERROR,
    #         )
    #
    #     mean_fragment_size = 100
    #     if len(frag_files) == 0 and prep != RNAType.TRNA.value:
    #         logger.warning(
    #             f"No fragment file found for {label}, using '100'. You should define this if you are going to use downstream zFPKM normalization"
    #         )
    #     elif len(frag_files) == 1:
    #         if layout == "single-end":
    #             mean_fragment_size = 0
    #         else:
    #             if not multi_flag:
    #                 frag_df = pd.read_table(frag_files[0], low_memory=False)
    #                 frag_df["meanxcount"] = frag_df["frag_mean"] * frag_df["frag_count"]
    #                 mean_fragment_size = sum(frag_df["meanxcount"] / sum(frag_df["frag_count"]))
    #
    #             else:
    #                 mean_fragment_sizes = np.array([])
    #                 library_sizes = np.array([])
    #                 for ff in frag_files:
    #                     frag_df = pd.read_table(ff, low_memory=False, sep="\t", on_bad_lines="skip")
    #                     frag_df["meanxcount"] = frag_df["frag_mean"] * frag_df["frag_count"]
    #                     mean_fragment_size = sum(frag_df["meanxcount"] / sum(frag_df["frag_count"]))
    #                     mean_fragment_sizes = np.append(mean_fragment_sizes, mean_fragment_size)
    #                     library_sizes = np.append(library_sizes, sum(frag_df["frag_count"]))
    #
    #                 mean_fragment_size = sum(mean_fragment_sizes * library_sizes) / sum(library_sizes)
    #     elif len(frag_files) > 1:
    #         _log_and_raise_error(
    #             f"Multiple matching fragment files for {label}, make sure there is only one copy for each replicate in COMO_input",
    #             error=ValueError,
    #             level=LogLevel.ERROR,
    #         )
    #
    #     sample_names.append(f"{context_name}_{study_number}{rep_number}")
    #     fragment_lengths.append(mean_fragment_size)
    #     layouts.append(layout)
    #     strands.append(strand)
    #     groups.append(study_number)
    #     preparation_method.append(prep)
    #
    # out_df = pd.DataFrame(
    #     {
    #         "sample_name": sample_names,
    #         "fragment_length": fragment_lengths,
    #         "layout": layouts,
    #         "strand": strands,
    #         "study": groups,
    #         "library_prep": preparation_method,
    #     }
    # ).sort_values("sample_name")
    # return out_df


async def _create_gene_info_file(
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
        data = _read_file(file, h5ad_as_df=False)
        if isinstance(data, pd.DataFrame):
            data: pd.DataFrame
            return data["ensembl_gene_id"].tolist()
        try:
            conversion = await gene_symbol_to_ensembl_and_gene_id(symbols=data.var_names.tolist(), taxon=taxon)
        except json.JSONDecodeError as e:
            _log_and_raise_error(
                f"Got a JSON decode error for file '{counts_matrix_filepaths}' ({e})",
                error=ValueError,
                level=LogLevel.CRITICAL,
            )

        # Remove NA values from entrez_gene_id dataframe column
        return conversion["ensembl_gene_id"].tolist()

    logger.info("Fetching gene info - this can take up to 5 minutes depending on the number of genes and your internet connection")

    ensembl_ids: set[str] = set(chain.from_iterable(await asyncio.gather(*[read_ensembl_gene_ids(f) for f in counts_matrix_filepaths])))
    gene_data: list[dict[str, str | int | list[str] | list[int] | None]] = await MyGene(cache=cache).query(
        items=list(ensembl_ids),
        taxon=taxon,
        scopes="ensemblgene",
    )
    gene_info: pd.DataFrame = pd.DataFrame(
        data=None,
        columns=pd.Index(data=["ensembl_gene_id", "gene_symbol", "entrez_gene_id", "size"]),
        index=pd.Index(data=list(range(len(ensembl_ids)))),
    )

    for i, data in enumerate(gene_data):
        data: dict[str, str | int | list[str] | list[int] | None]
        ensembl_genes: str | list[str] = cast(str | list[str], data.get("ensembl.gene", "-"))
        start_pos: int | list[int] = cast(int | list[int], data.get("genomic_pos.start", 0))
        end_pos: int | list[int] = cast(int | list[int], data.get("genomic_pos.end", 0))

        avg_start: int | float = sum(start_pos) / len(start_pos) if isinstance(start_pos, list) else start_pos
        avg_end: int | float = sum(end_pos) / len(end_pos) if isinstance(end_pos, list) else end_pos
        size: int = int(avg_end - avg_start)

        gene_info.at[i, "gene_symbol"] = data.get("symbol", "-")
        gene_info.at[i, "entrez_gene_id"] = data.get("entrezgene", "-")
        gene_info.at[i, "ensembl_gene_id"] = ",".join(ensembl_genes) if isinstance(ensembl_genes, list) else ensembl_genes
        gene_info.at[i, "size"] = size if size > 0 else -1

    gene_info["size"] = gene_info["size"].astype(str)  # replace no-length values with "-" to match rows where every value is "-"
    gene_info["size"] = gene_info["size"].replace("-1", "-")
    gene_info = cast(pd.DataFrame, gene_info[~(gene_info == "-").all(axis=1)])  # remove rows where every value is "-"

    gene_info["ensembl_gene_id"] = gene_info["ensembl_gene_id"].str.split(",")  # extend lists into multiple rows
    gene_info = gene_info.explode(column=["ensembl_gene_id"])
    gene_info["size"] = gene_info["size"].astype(int)
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
    if output_trna_config_filepath is not None and output_trna_fragment_lengths_filepath is not None:
        rna_types.append(
            (
                RNAType.TRNA,
                output_trna_config_filepath,
                output_trna_matrix_filepath,
                output_trna_fragment_lengths_filepath,
            )
        )
    if output_mrna_config_filepath is not None and output_mrna_fragment_lengths_filepath is not None:
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
            _log_and_raise_error(
                message="como_context_dir must be provided if create_gene_info_only is False",
                error=ValueError,
                level=LogLevel.ERROR,
            )
        if output_trna_fragment_lengths_filepath is None:
            _log_and_raise_error(
                message="output_fragment_lengths_filepath must be provided if create_gene_info_only is False",
                error=ValueError,
                level=LogLevel.ERROR,
            )

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


async def rnaseq_preprocess(
    context_name: str,
    taxon: int,
    output_gene_info_filepath: Path,
    como_context_dir: Path | None = None,
    input_matrix_filepath: Path | list[Path] | None = None,
    output_trna_fragment_lengths_filepath: Path | None = None,
    output_mrna_fragment_lengths_filepath: Path | None = None,
    output_trna_metadata_filepath: Path | None = None,
    output_mrna_metadata_filepath: Path | None = None,
    output_trna_count_matrix_filepath: Path | None = None,
    output_mrna_count_matrix_filepath: Path | None = None,
    cache: bool = True,
    log_level: LogLevel | str = LogLevel.INFO,
    log_location: str | io.TextIOWrapper = sys.stderr,
    *,
    create_gene_info_only: bool = False,
) -> None:
    """Preprocesses RNA-seq data for downstream analysis.

    Fetches additional gene information from a provided matrix or gene counts,
        or optionally creates this matrix using gene count files obtained using STAR aligner

    :param context_name: The context/cell type being processed
    :param taxon: The NCBI taxonomy ID
    :param output_gene_info_filepath: Path to the output gene information CSV file
    :param output_trna_fragment_lengths_filepath: Path to the output tRNA fragment lengths CSV file (if in "create" mode)
    :param output_mrna_fragment_lengths_filepath: Path to the output mRNA fragment lengths CSV file (if in "create" mode)
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
    _set_up_logging(level=log_level, location=log_location)

    output_gene_info_filepath = output_gene_info_filepath.resolve()

    if como_context_dir:
        como_context_dir = como_context_dir.resolve()
    input_matrix_filepath = [i.resolve() for i in _listify(input_matrix_filepath)] if input_matrix_filepath else None
    output_trna_metadata_filepath = output_trna_metadata_filepath.resolve() if output_trna_metadata_filepath else None
    output_mrna_metadata_filepath = output_mrna_metadata_filepath.resolve() if output_mrna_metadata_filepath else None
    output_trna_count_matrix_filepath = (
        output_trna_count_matrix_filepath.resolve() if output_trna_count_matrix_filepath else None
    )
    output_mrna_count_matrix_filepath = (
        output_mrna_count_matrix_filepath.resolve() if output_mrna_count_matrix_filepath else None
    )

    await _process(
        context_name=context_name,
        taxon=taxon,
        como_context_dir=como_context_dir,
        input_matrix_filepath=input_matrix_filepath,
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
