from __future__ import annotations

import asyncio
import contextlib
import re
import sys
from dataclasses import dataclass, field
from io import StringIO, TextIOWrapper
from itertools import chain
from pathlib import Path
from typing import Literal

import aiofiles
import numpy as np
import pandas as pd
import scanpy as sc
from fast_bioservices.biothings.mygene import MyGene
from fast_bioservices.pipeline import ensembl_to_gene_id_and_symbol, gene_symbol_to_ensembl_and_gene_id
from loguru import logger

from como.types import RNAPrepMethod, type_path, type_rna
from como.utils import _listify


@dataclass
class _STARinformation:
    num_unmapped: list[int]
    num_multimapping: list[int]
    num_no_feature: list[int]
    num_ambiguous: list[int]
    gene_names: list[str]
    count_matrix: pd.DataFrame

    @property
    def num_genes(self) -> int:
        return len(self.count_matrix)

    @classmethod
    async def build_from_tab(cls, filepath: Path) -> _STARinformation:
        if filepath.suffix != ".tab":
            raise ValueError(f"Building STAR information requires a '.tab' file; received: '{filepath}'")

        async with aiofiles.open(filepath) as i_stream:
            unmapped, multimapping, no_feature, ambiguous = await asyncio.gather(
                *[i_stream.readline(), i_stream.readline(), i_stream.readline(), i_stream.readline()]
            )
            num_unmapped = [int(i) for i in unmapped.rstrip("\n").split("\t")[1:]]
            num_multimapping = [int(i) for i in multimapping.rstrip("\n").split("\t")[1:]]
            num_no_feature = [int(i) for i in no_feature.rstrip("\n").split("\t")[1:]]
            num_ambiguous = [int(i) for i in ambiguous.rstrip("\n").split("\t")[1:]]
            remainder = await i_stream.read()

        string_io = StringIO(remainder)
        df = pd.read_csv(string_io, sep="\t", header=None)
        df.columns = [
            "ensembl_gene_id",
            "unstranded_rna_counts",
            "first_read_transcription_strand",
            "second_read_transcription_strand",
        ]
        df = df[~df["ensembl_gene_id"].isna()]
        return _STARinformation(
            num_unmapped=num_unmapped,
            num_multimapping=num_multimapping,
            num_no_feature=num_no_feature,
            num_ambiguous=num_ambiguous,
            gene_names=df["ensembl_gene_id"].values.tolist(),
            count_matrix=df,
        )


@dataclass
class _StudyMetrics:
    study_name: str
    count_files: list[Path]
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
        self.__num_samples = len(self.count_files)
        self.__sample_names = [f.stem for f in self.count_files]

        if len(self.count_files) != len(self.strand_files):
            raise ValueError(
                f"Unequal number of count files and strand files for study '{self.study_name}'. "
                f"Found {len(self.count_files)} count files and {len(self.strand_files)} strand files."
            )

        if self.num_samples != len(self.count_files):
            raise ValueError(
                f"Unequal number of samples and count files for study '{self.study_name}'. "
                f"Found {self.num_samples} samples and {len(self.count_files)} count files."
            )

        if self.num_samples != len(self.strand_files):
            raise ValueError(
                f"Unequal number of samples and strand files for study '{self.study_name}'. "
                f"Found {self.num_samples} samples and {len(self.strand_files)} strand files."
            )

        if self.__num_samples == 1:
            raise ValueError(f"Only one sample exists for study {self.study_name}. Provide at least two samples")

        self.count_files.sort()
        self.strand_files.sort()
        self.__sample_names.sort()


def _sample_name_from_filepath(file: Path) -> str:
    return re.search(r".+_S\d+R\d+", file.stem).group()


def _organize_gene_counts_files(data_dir: Path) -> list[_StudyMetrics]:
    gene_count_dir = Path(data_dir, "geneCounts")
    strand_dir = Path(data_dir, "strandedness")

    gene_counts_directories: list[Path] = sorted([p for p in gene_count_dir.glob("*") if not p.name.startswith(".")])
    strandedness_directories: list[Path] = sorted([p for p in strand_dir.glob("*") if not p.name.startswith(".")])

    if len(gene_counts_directories) != len(strandedness_directories):
        raise ValueError(
            f"Unequal number of gene count directories and strandedness directories. "
            f"Found {len(gene_counts_directories)} gene count directories and {len(strandedness_directories)} strandedness directories."  # noqa: E501
            f"\nGene count directory: {gene_count_dir}\nStrandedness directory: {strand_dir}"
        )

    # For each study, collect gene count files, fragment files, insert size files, layouts, and strandedness information
    study_metrics: list[_StudyMetrics] = []
    for gene_dir, strand_dir in zip(gene_counts_directories, strandedness_directories):
        if gene_dir.stem != strand_dir.stem:
            raise ValueError(
                f"Gene directory name of '{gene_dir.stem}' does not match stranded directory name of '{strand_dir.stem}'"  # noqa: E501
            )

        study_metrics.append(
            _StudyMetrics(
                study_name=gene_dir.stem,
                count_files=list(gene_dir.glob("*.tab")),
                strand_files=list(strand_dir.glob("*.txt")),
            )
        )
    return study_metrics


async def _process_first_multirun_sample(strand_file: Path, all_counts_files: list[Path]):
    sample_count = pd.DataFrame()
    for file in all_counts_files:
        star_information = await _STARinformation.build_from_tab(file)
        strand_information = strand_file.read_text().rstrip("\n").lower()

        if strand_information not in ("none", "first_read_transcription_strand", "second_read_transcription_strand"):
            raise ValueError(
                f"Unrecognized Strand Information: {strand_information}; "
                f"expected 'none', 'first_read_transcription_strand', or 'second_read_transcription_strand'"
            )

        if strand_information == "none":
            strand_information = "unstranded_rna_counts"

        run_counts = star_information.count_matrix[["ensembl_gene_id", strand_information]]
        run_counts.columns = pd.Index(["ensembl_gene_id", "counts"])
        sample_count = (
            run_counts
            if sample_count.empty
            else sample_count.merge(run_counts, on=["ensembl_gene_id", "counts"], how="outer")
        )

    # Set na values to 0
    sample_count = sample_count.fillna(value="0")
    sample_count.iloc[:, 1:] = sample_count.iloc[:, 1:].apply(pd.to_numeric)

    count_sums: pd.DataFrame = pd.DataFrame(sample_count.sum(axis=1, numeric_only=True))
    count_sums.insert(0, "ensembl_gene_id", sample_count["ensembl_gene_id"])
    count_sums.columns = pd.Index(["ensembl_gene_id", _sample_name_from_filepath(strand_file)])
    return count_sums


async def _process_standard_replicate(counts_file: Path, strand_file: Path, sample_name: str):
    star_information = await _STARinformation.build_from_tab(counts_file)
    strand_information = strand_file.read_text().rstrip("\n").lower()

    if strand_information not in ("none", "first_read_transcription_strand", "second_read_transcription_strand"):
        raise ValueError(
            f"Unrecognized Strand Information: {strand_information}; "
            f"expected 'none', 'first_read_transcription_strand', or 'second_read_transcription_strand'"
        )

    if strand_information == "none":
        strand_information = "unstranded_rna_counts"

    sample_count = star_information.count_matrix[["ensembl_gene_id", strand_information]]
    sample_count.columns = pd.Index(["ensembl_gene_id", sample_name])
    return sample_count


async def _prepare_sample_counts(
    sample_name: str,
    counts_file: Path,
    strand_file: Path,
    all_counts_files: list[Path],
) -> pd.DataFrame | Literal["SKIP"]:
    # Test if the counts_file is the first run in a multi-run smaple
    if re.search(r"R\d+r1", counts_file.as_posix()):
        return await _process_first_multirun_sample(strand_file=strand_file, all_counts_files=all_counts_files)
    elif re.search(r"R\d+r\d+", counts_file.as_posix()):
        return "SKIP"
    else:
        return await _process_standard_replicate(counts_file, strand_file, sample_name)


async def _create_sample_counts_matrix(metrics: _StudyMetrics) -> pd.DataFrame:
    adjusted_index = 0
    counts: pd.DataFrame | Literal["SKIP"] = await _prepare_sample_counts(
        sample_name=metrics.sample_names[0],
        counts_file=metrics.count_files[0],
        strand_file=metrics.strand_files[0],
        all_counts_files=metrics.count_files,
    )

    for i in range(1, metrics.num_samples):
        new_counts = await _prepare_sample_counts(
            sample_name=metrics.sample_names[i],
            counts_file=metrics.count_files[i],
            strand_file=metrics.strand_files[i],
            all_counts_files=metrics.count_files,
        )
        if isinstance(new_counts, str) and new_counts == "SKIP":
            adjusted_index += 1
            continue

        counts: pd.DataFrame = counts.merge(new_counts, on="ensembl_gene_id", how="outer")
        counts = counts.fillna(value=0)

        # Remove run number "r\d+" from multi-run names
        if re.search(r"R\d+r1", metrics.sample_names[i]):
            new_sample_name = re.sub(r"r\d+", "", metrics.sample_names[i])
            counts.columns[i + 1 - adjusted_index] = new_sample_name

    return counts


async def _write_counts_matrix(
    *,
    config_df: pd.DataFrame,
    como_context_dir: Path,
    output_counts_matrix_filepath: Path,
    rna: type_rna,
) -> pd.DataFrame:
    """Create a counts matrix file by reading gene counts table(s)."""
    study_metrics = _organize_gene_counts_files(data_dir=como_context_dir)
    counts: list[pd.DataFrame] = await asyncio.gather(
        *[_create_sample_counts_matrix(metric) for metric in study_metrics]
    )
    final_matrix = pd.DataFrame()
    for count in counts:
        final_matrix = count if final_matrix.empty else pd.merge(final_matrix, count, on="ensembl_gene_id", how="outer")

    rna_specific_sample_names = config_df.loc[config_df["library_prep"] == rna, "sample_name"].tolist()
    final_matrix = final_matrix[["ensembl_gene_id", *rna_specific_sample_names]]

    output_counts_matrix_filepath.parent.mkdir(parents=True, exist_ok=True)
    final_matrix.to_csv(output_counts_matrix_filepath, index=False)
    logger.success(f"Wrote gene count matrix for '{rna}' RNA at '{output_counts_matrix_filepath}'")
    return final_matrix


async def _create_config_df(context_name: str, /, como_input_dir: Path) -> pd.DataFrame:  # noqa: C901
    """Create configuration sheet.

    The configuration file created is based on the gene counts matrix.
     If using zFPKM normalization technique, mean fragment lengths will be fetched
    """
    gene_counts_files = list(Path(como_input_dir, context_name, "geneCounts").rglob("*.tab"))
    sample_names: list[str] = []
    fragment_lengths: list[int | float] = []
    layouts: list[str] = []
    strands: list[str] = []
    groups: list[str] = []
    preparation_method: list[str] = []

    for gene_count_filename in sorted(gene_counts_files):
        try:
            # Match S___R___r___
            # \d{1,3} matches 1-3 digits
            # (?:r\d{1,3})? matches an option "r" followed by three digits
            label = re.findall(r"S\d{1,3}R\d{1,3}(?:r\d{1,3})?", gene_count_filename.as_posix())[0]

        except IndexError as e:
            raise IndexError(
                f"\n\nFilename of '{gene_count_filename}' is not valid. "
                f"Should be 'contextName_SXRYrZ.tab', where X is the study/batch number, Y is the replicate number, "
                f"and Z is the run number."
                "\n\nIf not a multi-run sample, exclude 'rZ' from the filename."
            ) from e

        study_number = re.findall(r"S\d{1,3}", label)[0]
        rep_number = re.findall(r"R\d{1,3}", label)[0]
        run = re.findall(r"r\d{1,3}", label)

        multi_flag = 0
        if len(run) > 0:
            if run[0] != "r1":
                continue
            else:
                label_glob = study_number + rep_number + "r*"
                runs = [run for run in gene_counts_files if re.search(label_glob, run.as_posix())]
                multi_flag = 1
                frag_files = []

                for r in runs:
                    r_label = re.findall(r"r\d{1,3}", r.as_posix())[0]
                    R_label = re.findall(r"R\d{1,3}", r.as_posix())[0]  # noqa: N806
                    frag_filename = "".join([context_name, "_", study_number, R_label, r_label, "_fragment_size.txt"])
                    frag_files.append(como_input_dir / context_name / "fragmentSizes" / study_number / frag_filename)

        context_path = como_input_dir / context_name
        layout_files: list[Path] = list((context_path / "layouts").rglob(f"{context_name}_{label}_layout.txt"))
        strand_files: list[Path] = list((context_path / "strandedness").rglob(f"{context_name}_{label}_strandedness.txt"))  # fmt: skip  # noqa: E501
        frag_files: list[Path] = list((context_path / "fragmentSizes").rglob(f"{context_name}_{label}_fragment_size.txt"))  # fmt: skip  # noqa: E501
        prep_files: list[Path] = list((context_path / "prepMethods").rglob(f"{context_name}_{label}_prep_method.txt"))

        layout = "UNKNOWN"
        if len(layout_files) == 0:
            logger.warning(
                f"No layout file found for {label}, writing as 'UNKNOWN', "
                f"this should be defined by user if using zFPKM or rnaseq_gen.py will not run"
            )
        elif len(layout_files) == 1:
            with layout_files[0].open("r") as file:
                layout = file.read().strip()
        elif len(layout_files) > 1:
            raise ValueError(
                f"Multiple matching layout files for {label}, "
                f"make sure there is only one copy for each replicate in COMO_input"
            )

        strand = "UNKNOWN"
        if len(strand_files) == 0:
            logger.warning(
                f"No strandedness file found for {label}, writing as 'UNKNOWN'. "
                f"This will not interfere with the analysis since you have already set rnaseq_preprocess.py to "
                f"infer the strandedness when writing the counts matrix"
            )
        elif len(strand_files) == 1:
            with strand_files[0].open("r") as file:
                strand = file.read().strip()
        elif len(strand_files) > 1:
            raise ValueError(
                f"Multiple matching strandedness files for {label}, "
                f"make sure there is only one copy for each replicate in COMO_input"
            )

        prep = "total"
        if len(prep_files) == 0:
            logger.warning(f"No prep file found for {label}, assuming 'total' as in Total RNA library preparation")
        elif len(prep_files) == 1:
            with prep_files[0].open("r") as file:
                prep = file.read().strip().lower()
                if prep not in ["total", "mrna"]:
                    raise ValueError(f"Prep method must be either 'total' or 'mrna' for {label}")
        elif len(prep_files) > 1:
            raise ValueError(
                f"Multiple matching prep files for {label}, "
                f"make sure there is only one copy for each replicate in COMO_input"
            )

        mean_fragment_size = 100
        if len(frag_files) == 0:
            logger.warning(
                f"No fragment file found for {label}, using '100'. "
                f"This must be defined by the user in order to use zFPKM normalization"
            )
        elif len(frag_files) == 1:
            if layout == "single-end":
                mean_fragment_size = 0
            else:
                if not multi_flag:
                    frag_df = pd.read_table(frag_files[0], low_memory=False)
                    frag_df["meanxcount"] = frag_df["frag_mean"] * frag_df["frag_count"]
                    mean_fragment_size = sum(frag_df["meanxcount"] / sum(frag_df["frag_count"]))

                else:
                    mean_fragment_sizes = np.array([])
                    library_sizes = np.array([])
                    for ff in frag_files:
                        frag_df = pd.read_table(ff, low_memory=False, sep="\t", on_bad_lines="skip")
                        frag_df["meanxcount"] = frag_df["frag_mean"] * frag_df["frag_count"]
                        mean_fragment_size = sum(frag_df["meanxcount"] / sum(frag_df["frag_count"]))
                        mean_fragment_sizes = np.append(mean_fragment_sizes, mean_fragment_size)
                        library_sizes = np.append(library_sizes, sum(frag_df["frag_count"]))

                    mean_fragment_size = sum(mean_fragment_sizes * library_sizes) / sum(library_sizes)
        elif len(frag_files) > 1:
            raise ValueError(
                f"Multiple matching fragment files for {label}, "
                f"make sure there is only one copy for each replicate in COMO_input"
            )

        sample_names.append(f"{context_name}_{study_number}{rep_number}")
        fragment_lengths.append(mean_fragment_size)
        layouts.append(layout)
        strands.append(strand)
        groups.append(study_number)
        preparation_method.append(prep)

    out_df = pd.DataFrame(
        {
            "sample_name": sample_names,
            "fragment_length": fragment_lengths,
            "layout": layouts,
            "strand": strands,
            "study": groups,
            "library_prep": preparation_method,
        }
    ).sort_values("sample_name")
    return out_df


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

    async def read_counts(file: Path) -> list[str]:
        data = await asyncio.to_thread(pd.read_csv if file.suffix == ".csv" else sc.read_h5ad, file)
        conversion = (
            await ensembl_to_gene_id_and_symbol(ids=data["ensembl_gene_id"].tolist(), taxon=taxon)
            if isinstance(data, pd.DataFrame)
            else await gene_symbol_to_ensembl_and_gene_id(symbols=data.var_names.tolist(), taxon=taxon)
        )
        return conversion["entrez_gene_id"].tolist()

    logger.info("Fetching gene info (this may take 1-5 minutes)")
    genes = set(chain.from_iterable(await asyncio.gather(*[read_counts(f) for f in counts_matrix_filepaths])))

    mygene = MyGene(cache=cache)
    gene_data = await mygene.query(items=list(genes), taxon=taxon, scopes="entrezgene")
    gene_info: pd.DataFrame = pd.DataFrame(
        data=None,
        columns=pd.Index(data=["ensembl_gene_id", "gene_symbol", "entrez_gene_id", "start_position", "end_position"]),
        index=pd.Index(data=range(len(gene_data))),
    )
    for i, data in enumerate(gene_data):
        ensembl_ids = data.get("ensembl.gene", "-")
        if isinstance(ensembl_ids, list):
            ensembl_ids = ensembl_ids[0]

        start_pos = data.get("genomic_pos.start", 0)
        start_pos = sum(start_pos) / len(start_pos) if isinstance(start_pos, list) else start_pos
        end_pos = data.get("genomic_pos.end", 0)
        end_pos = sum(end_pos) / len(end_pos) if isinstance(end_pos, list) else end_pos

        gene_info.at[i, "gene_symbol"] = data.get("symbol", "-")
        gene_info.at[i, "entrez_gene_id"] = data.get("entrezgene", "-")
        gene_info.at[i, "ensembl_gene_id"] = ensembl_ids
        gene_info.at[i, "start_position"] = start_pos
        gene_info.at[i, "end_position"] = end_pos

    gene_info = gene_info[
        (
            (gene_info["entrez_gene_id"] != "-")
            & (gene_info["ensembl_gene_id"] != "-")
            & (gene_info["gene_symbol"] != "-")
        )
    ]
    gene_info["size"] = gene_info["end_position"].astype(int) - gene_info["start_position"].astype(int)
    gene_info.drop(columns=["start_position", "end_position"], inplace=True)
    gene_info.sort_values(by="ensembl_gene_id", inplace=True)
    gene_info.to_csv(output_filepath, index=False)
    logger.success(f"Gene Info file written at '{output_filepath}'")


async def _create_matrix_file(
    context_name: str,
    output_config_filepath: Path,
    como_context_dir: type_path,
    output_counts_matrix_filepath: Path,
    rna: type_rna,
) -> None:
    config_df = await _create_config_df(context_name, como_input_dir=como_context_dir)
    await _write_counts_matrix(
        config_df=config_df,
        como_context_dir=como_context_dir,
        output_counts_matrix_filepath=output_counts_matrix_filepath,
        rna=rna,
    )
    with pd.ExcelWriter(output_config_filepath) as writer:
        subset_config = config_df[config_df["library_prep"] == rna]
        subset_config.to_excel(writer, sheet_name=context_name, header=True, index=False)


async def _process(
    context_name: str,
    taxon: int,
    output_gene_info_filepath: Path,
    como_context_dir: Path | None,
    input_matrix_filepath: list[Path] | None,
    output_trna_config_filepath: Path | None,
    output_mrna_config_filepath: Path | None,
    output_trna_matrix_filepath: Path | None,
    output_mrna_matrix_filepath: Path | None,
    cache: bool,
):
    rna_types: list[tuple[type_rna, Path, Path]] = []
    if output_trna_config_filepath:
        rna_types.append(("total", output_trna_config_filepath, output_trna_matrix_filepath))
    if output_mrna_config_filepath:
        rna_types.append(("polya", output_mrna_config_filepath, output_mrna_matrix_filepath))

    # if provided, iterate through como-input specific directories
    tasks = []
    for rna, output_config_filepath, output_matrix_filepath in rna_types:
        tasks.append(
            asyncio.create_task(
                _create_matrix_file(
                    context_name=context_name,
                    output_config_filepath=output_config_filepath,
                    como_context_dir=como_context_dir,
                    output_counts_matrix_filepath=output_matrix_filepath,
                    rna=rna,
                )
            )
        )

    await asyncio.gather(*tasks)

    # create the gene info filepath based on provided data
    await _create_gene_info_file(
        counts_matrix_filepaths=[
            f
            for f in [*input_matrix_filepath, output_trna_matrix_filepath, output_mrna_matrix_filepath]
            if f is not None
        ],
        output_filepath=output_gene_info_filepath,
        taxon=taxon,
        cache=cache,
    )




async def rnaseq_preprocess(  # noqa: C901
    context_names: str | list[str],
    mode: Literal["create", "provide"] | list[Literal["create", "provide"]],
    taxon_id: type_taxon | list[type_taxon],
    input_como_dirpath: type_path | list[type_path] | None = None,
    input_matrix_filepath: type_path | list[type_path] | None = None,
    output_gene_info_filepath: Path | None = None,
    output_trna_config_filepath: Path | None = None,
    output_mrna_config_filepath: Path | None = None,
    output_count_matrices_dirpath: list[Path] | None = None,
    output_trna_count_matrix: list[Path] | None = None,
    output_mrna_count_matrix: list[Path] | None = None,
    cache: bool = True,
) -> None:
    """Preprocesses RNA-seq data for downstream analysis.

    Fetches additional gene information from a provided matrix or gene counts,
        or optionally creates this matrix using gene count files obtained using STAR aligner

    :param context_name: The context/cell type being processed
    :param taxon: The NCBI taxonomy ID
    :param output_gene_info_filepath: Path to the output gene information CSV file
    :param output_trna_config_filepath: Path to the output tRNA config file (if in "create" mode)
    :param output_polya_config_filepath: Path to the output mRNA config file (if in "create" mode)
    :param output_trna_count_matrix_filepath: The path to write total RNA count matrices
    :param output_polya_count_matrix_filepath: The path to write messenger RNA count matrices
    :param como_context_dir: If in "create" mode, the input path(s) to the COMO_input directory of the current context
        i.e., the directory containing "fragmentSizes", "geneCounts", "insertSizeMetrics", etc. directories
    :param input_matrix_filepath: If in "provide" mode, the path(s) to the count matrices to be processed
    :param preparation_method: The preparation method
    :param cache: Should HTTP requests be cached
    :param log_level: The logging level
    :param log_location: The logging location
    """
    with contextlib.suppress(ValueError):
        logger.remove(0)
        logger.add(
            sink=log_location,
            level=log_level,
            format="<green>{time:YYYY-MM-DD HH:mm:ss}</green> | <level>{level: <8}</level> | <cyan>{name}</cyan>:<cyan>{line}</cyan> - <level>{message}</level>",  # noqa: E501
        )

    output_gene_info_filepath = output_gene_info_filepath.resolve()
    como_context_dir = como_context_dir.resolve()
    input_matrix_filepath = [i.resolve() for i in _listify(input_matrix_filepath)] if input_matrix_filepath else None
    output_trna_config_filepath = (
        output_trna_config_filepath.resolve() if output_trna_config_filepath else output_trna_config_filepath
    )
    output_polya_config_filepath = (
        output_polya_config_filepath.resolve() if output_polya_config_filepath else output_polya_config_filepath
    )
    output_trna_count_matrix_filepath = (
        output_trna_count_matrix_filepath.resolve()
        if output_trna_count_matrix_filepath
        else output_trna_count_matrix_filepath
    )
    output_polya_count_matrix_filepath = (
        output_polya_count_matrix_filepath.resolve()
        if output_polya_count_matrix_filepath
        else output_polya_count_matrix_filepath
    )

    input_matrix_filepath = _listify(input_matrix_filepath)
    preparation_method = _listify(preparation_method)

    if len(input_matrix_filepath) != len(preparation_method):
        raise ValueError(
            "input_matrix_filepath (--input-matrix-filepath) and "
            "preparation_method (--preparation-method) must be the same length."
        )
    await _process(
        context_name=context_name,
        taxon=taxon,
        como_context_dir=como_context_dir,
        input_matrix_filepath=input_matrix_filepath,
        output_gene_info_filepath=output_gene_info_filepath,
        output_trna_config_filepath=output_trna_config_filepath,
        output_mrna_config_filepath=output_polya_config_filepath,
        output_trna_matrix_filepath=output_trna_count_matrix_filepath,
        output_mrna_matrix_filepath=output_polya_count_matrix_filepath,
        cache=cache,
    )
