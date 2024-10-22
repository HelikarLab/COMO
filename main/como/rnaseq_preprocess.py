import sys
from pathlib import Path

sys.path.insert(0, Path(__file__).parent.parent.as_posix())

import argparse
import re
from dataclasses import dataclass, field
from typing import Literal, Optional, Union

import numpy as np
import pandas as pd
from fast_bioservices import BioDBNet, Input, Output, Taxon
from loguru import logger

from como import Config, stringlist_to_list


@dataclass
class _Arguments:
    context_names: list[str]
    gene_format: Input
    taxon_id: Taxon | int | str
    mode: Literal["create", "provide"]
    provided_matrix_fname: str = None

    def __post_init__(self):
        if self.mode == "provide" and self.provided_matrix_fname is None:
            raise ValueError("If provide_matrix is True, then provided_matrix_fname must be provided")

        if self.gene_format.upper() in ["ENSEMBL", "ENSEMBLE", "ENSG", "ENSMUSG", "ENSEMBL ID", "ENSEMBL GENE ID"]:
            self.gene_format: Input = Input.ENSEMBL_GENE_ID
        elif self.gene_format.upper() in ["HGNC SYMBOL", "HUGO", "HUGO SYMBOL", "SYMBOL", "HGNC", "GENE SYMBOL"]:
            self.gene_format: Input = Input.GENE_SYMBOL
        elif self.gene_format.upper() in ["ENTREZ", "ENTRES", "ENTREZ ID", "ENTREZ NUMBER" "GENE ID"]:
            self.gene_format: Input = Input.GENE_ID
        else:  # provided invalid gene format
            raise ValueError(f"Gene format (--gene_format) is invalid; accepts 'Ensembl', 'Entrez', and 'HGNC symbol'; provided: {self.gene_format}")

        # handle species alternative ids
        if isinstance(self.taxon_id, str) and not self.taxon_id.isdigit():
            if self.taxon_id.upper() in ["HUMAN", "HOMO SAPIENS"]:
                self.taxon_id = Taxon.HOMO_SAPIENS
            elif self.taxon_id.upper() in ["MOUSE", "MUS MUSCULUS"]:
                self.taxon_id = Taxon.MUS_MUSCULUS
            else:
                raise ValueError(f"Taxon id (--taxon-id) is invalid; accepts 'human', 'mouse', or an integer value; provided: {self.taxon_id}")
        else:
            try:
                # If taxon id can't be found in the Taxon Enum, do nothing
                self.taxon_id = Taxon.from_int(int(self.taxon_id))
            except ValueError:
                pass


@dataclass
class _STARinformation:
    num_unmapped: list[int]
    num_multimapping: list[int]
    num_no_feature: list[int]
    num_ambiguous: list[int]
    gene_names: list[str]
    count_matrix: pd.DataFrame

    @property
    def num_genes(self):
        return len(self.count_matrix)

    @classmethod
    def build_from_tab(cls, filepath: Path) -> "_STARinformation":
        if filepath.suffix != ".tab":
            raise ValueError(f"Building STAR information requires a '.tab' file; received: '{filepath}'")
        with open(filepath) as i_stream:
            num_unmapped = [int(i) for i in next(i_stream).rstrip("\n").split("\t")[1:]]
            num_multimapping = [int(i) for i in next(i_stream).rstrip("\n").split("\t")[1:]]
            num_no_feature = [int(i) for i in next(i_stream).rstrip("\n").split("\t")[1:]]
            num_ambiguous = [int(i) for i in next(i_stream).rstrip("\n").split("\t")[1:]]

        df = pd.read_csv(
            filepath,
            sep="\t",
            skiprows=4,
            names=["ensembl_gene_id", "unstranded_rna_counts", "first_read_transcription_strand", "second_read_transcription_strand"],
        )
        # Remove NA values
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
    strandedness_files: list[Path]
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

        if len(self.count_files) != len(self.strandedness_files):
            raise ValueError(
                f"Unequal number of count files and strand files for study '{self.study_name}'. Found {len(self.count_files)} count files and {len(self.strandedness_files)} strand files."
            )

        if self.num_samples != len(self.count_files):
            raise ValueError(
                f"Unequal number of samples and count files for study '{self.study_name}'. Found {self.num_samples} samples and {len(self.count_files)} count files."
            )

        if self.num_samples != len(self.strandedness_files):
            raise ValueError(
                f"Unequal number of samples and strand files for study '{self.study_name}'. Found {self.num_samples} samples and {len(self.strandedness_files)} strand files."
            )

        if self.__num_samples == 1:
            raise ValueError(f"Only one sample exists for study {self.study_name}. Provide at least two samples")

        self.count_files.sort()
        self.strandedness_files.sort()
        self.__sample_names.sort()


def _context_from_filepath(file: Path) -> str:
    return file.stem.split("_S")[0]


def _sample_name_from_filepath(file: Path) -> str:
    return re.search(r".+_S\d+R\d+", file.stem).group()


def _organize_gene_counts_files(data_dir: Path) -> list[_StudyMetrics]:
    root_gene_count_dir = Path(data_dir, "geneCounts")
    root_strandedness_dir = Path(data_dir, "strandedness")

    gene_counts_directories: list[Path] = sorted(list(Path(root_gene_count_dir).glob("*")))
    strandedness_directories: list[Path] = sorted(list(Path(root_strandedness_dir).glob("*")))

    if len(gene_counts_directories) != len(strandedness_directories):
        raise ValueError(
            f"Unequal number of gene count directories and strandedness directories. Found {len(gene_counts_directories)} gene count directories and {len(strandedness_directories)} strandedness directories.\nGene count directory: {root_gene_count_dir}\nStrandedness directory: {root_strandedness_dir}"
        )

    # For each study, collect gene count files, fragment files, insert size files, layouts, and strandedness information
    study_metrics: list[_StudyMetrics] = []
    for gene_dir, strand_dir in zip(gene_counts_directories, strandedness_directories):
        if gene_dir.stem != strand_dir.stem:
            raise ValueError(f"Gene directory name of '{gene_dir.stem}' does not match stranded directory name of '{strand_dir.stem}'")

        study_metrics.append(
            _StudyMetrics(
                study_name=gene_dir.stem,
                count_files=list(gene_dir.glob("*.tab")),
                strandedness_files=list(strand_dir.glob("*.txt")),
            )
        )
    return study_metrics


def _process_first_multirun_sample(strand_file: Path, all_counts_files: list[Path]):
    sample_count = pd.DataFrame()
    for file in all_counts_files:
        star_information = _STARinformation.build_from_tab(file)
        strand_information = strand_file.read_text().rstrip("\n").lower()

        if strand_information not in ("none", "first_read_transcription_strand", "second_read_transcription_strand"):
            raise ValueError(
                f"Unrecognized Strand Information: {strand_information}; expected 'none', 'first_read_transcription_strand', or 'second_read_transcription_strand'"
            )

        if strand_information == "none":
            strand_information = "unstranded_rna_counts"

        run_counts = star_information.count_matrix[["gene_id", strand_information]]
        run_counts.columns = pd.Index(["ensembl_gene_id", "counts"])
        if sample_count.empty:
            sample_count = run_counts
        else:
            # Merge to take all items from both data frames
            sample_count = sample_count.merge(run_counts, on="ensembl_gene_id", how="outer")

    # Set na values to 0
    sample_count = sample_count.fillna(value="0")
    sample_count.iloc[:, 1:] = sample_count.iloc[:, 1:].apply(pd.to_numeric)

    count_sums: pd.DataFrame = pd.DataFrame(sample_count.sum(axis=1, numeric_only=True))
    count_sums.insert(0, "ensembl_gene_id", sample_count["ensembl_gene_id"])
    count_sums.columns = pd.Index(["ensembl_gene_id", _sample_name_from_filepath(strand_file)])
    return count_sums


def _process_standard_replicate(counts_file: Path, strand_file: Path, sample_name: str):
    star_information = _STARinformation.build_from_tab(counts_file)
    strand_information = strand_file.read_text().rstrip("\n").lower()

    if strand_information not in ("none", "first_read_transcription_strand", "second_read_transcription_strand"):
        raise ValueError(
            f"Unrecognized Strand Information: {strand_information}; expected 'none', 'first_read_transcription_strand', or 'second_read_transcription_strand'"
        )

    if strand_information == "none":
        strand_information = "unstranded_rna_counts"

    sample_count = star_information.count_matrix[["ensembl_gene_id", strand_information]]
    sample_count.columns = pd.Index(["ensembl_gene_id", sample_name])
    return sample_count


def _prepare_sample_counts(sample_name: str, counts_file: Path, strand_file: Path, all_counts_files: list[Path]) -> pd.DataFrame | Literal["SKIP"]:
    # Test if the counts_file is the first run in a multi-run smaple
    if re.search(r"R\d+r1", counts_file.as_posix()):
        return _process_first_multirun_sample(strand_file=strand_file, all_counts_files=all_counts_files)
    elif re.search(r"R\d+r\d+", counts_file.as_posix()):
        return "SKIP"
    else:
        return _process_standard_replicate(counts_file, strand_file, sample_name)


def _create_sample_counts_matrix(metrics: _StudyMetrics) -> pd.DataFrame:
    adjusted_index = 0
    counts: pd.DataFrame | Literal["SKIP"] = _prepare_sample_counts(
        sample_name=metrics.sample_names[0],
        counts_file=metrics.count_files[0],
        strand_file=metrics.strandedness_files[0],
        all_counts_files=metrics.count_files,
    )

    for i in range(1, metrics.num_samples):
        new_counts = _prepare_sample_counts(
            sample_name=metrics.sample_names[i],
            counts_file=metrics.count_files[i],
            strand_file=metrics.strandedness_files[i],
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


def _create_context_counts_matrix(data_dir: Path, output_dir: Path):
    study_metrics = _organize_gene_counts_files(data_dir=data_dir)
    final_matrix: pd.DataFrame = pd.DataFrame()
    for metric in study_metrics:
        counts: pd.DataFrame = _create_sample_counts_matrix(metric)
        final_matrix = counts if final_matrix.empty else pd.merge(final_matrix, counts, on="ensembl_gene_id", how="outer")

    output_filename = output_dir / f"gene_counts_matrix_full_{data_dir.stem}.csv"
    output_filename.parent.mkdir(parents=True, exist_ok=True)
    logger.info(f"Writing context '{data_dir.stem}' gene count matrix: {output_filename}")
    final_matrix.to_csv(output_filename, index=False)


def _create_counts_matrix(context_name: str, config: Config):
    """
    Create a counts matrix by reading gene counts tables in COMO_input/<context name>/<study number>/geneCounts/
    Uses R in backend (_create_context_counts_matrix.R)
    """
    config = Config()
    input_dir = config.data_dir / "COMO_input" / context_name
    matrix_output_dir = config.data_dir / "data_matrices" / context_name
    _create_context_counts_matrix(data_dir=input_dir, output_dir=matrix_output_dir)


def _create_config_df(context_name: str) -> pd.DataFrame:
    """
    Create configuration sheet at /main/data/config_sheets/rnaseq_data_inputs_auto.xlsx
    based on the gene counts matrix. If using zFPKM normalization technique, fetch mean fragment lengths from
    /work/data/COMO_input/<context name>/<study number>/fragmentSizes/
    """
    config = Config()
    gene_counts_files = list(Path(config.data_dir, "COMO_input", context_name, "geneCounts").rglob("*.tab"))

    sample_names: list[str] = []
    fragment_lengths: list[int | float] = []
    layouts: list[str] = []
    strands: list[str] = []
    groups: list[str] = []
    preparation_method: list[str] = []

    for gcfilename in sorted(gene_counts_files):
        try:
            # Match S___R___r___
            # \d{1,3} matches 1-3 digits
            # (?:r\d{1,3})? matches an option "r" followed by three digits
            label = re.findall(r"S\d{1,3}R\d{1,3}(?:r\d{1,3})?", gcfilename.as_posix())[0]

        except IndexError:
            raise IndexError(
                f"\n\nFilename of '{gcfilename}' is not valid. Should be 'contextName_SXRYrZ.tab', where X is the "
                "study/batch number, Y is the replicate number, and Z is the run number."
                "\n\nIf not a multi-run sample, exclude 'rZ' from the filename."
            )

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
                    R_label = re.findall(r"R\d{1,3}", r.as_posix())[0]
                    frag_filename = "".join([context_name, "_", study_number, R_label, r_label, "_fragment_size.txt"])
                    frag_files.append(config.data_dir / "COMO_input" / context_name / "fragmentSizes" / study_number / frag_filename)

        context_path = config.data_dir / "COMO_input" / context_name
        layout_files = list(Path(context_path, "layouts").rglob(f"{context_name}_{label}_layou.txt"))
        strand_files = list(Path(context_path, "strandedness").rglob(f"{context_name}_{label}_strandedness.txt"))
        frag_files = list(Path(context_path, "fragmentSizes").rglob(f"{context_name}_{label}_fragment_size.txt"))
        prep_files = list(Path(context_path, "prepMethods").rglob(f"{context_name}_{label}_prep_method.txt"))

        # Get layout
        layout = "UNKNOWN"
        if len(layout_files) == 0:
            logger.warning(
                f"No layout file found for {label}, writing as 'UNKNOWN', this should be defined by user if using zFPKM or rnaseq_gen.py will not run"
            )
        elif len(layout_files) == 1:
            with open(layout_files[0]) as file:
                layout = file.read().strip()
        elif len(layout_files) > 1:
            raise ValueError(f"Multiple matching layout files for {label}, make sure there is only one copy for each replicate in COMO_input")

        # Get strandedness
        strand = "UNKNOWN"
        if len(strand_files) == 0:
            logger.warning(
                f"No strandedness file found for {label}, writing as 'UNKNOWN'. This will not interfere with the analysis since you have already set rnaseq_preprocess.py to infer the strandedness when writing the counts matrix"
            )
        elif len(strand_files) == 1:
            with open(strand_files[0]) as file:
                strand = file.read().strip()
        elif len(strand_files) > 1:
            raise ValueError(f"Multiple matching strandedness files for {label}, make sure there is only one copy for each replicate in COMO_input")

        # Get preparation method
        prep = "total"
        if len(prep_files) == 0:
            logger.warning(f"No prep file found for {label}, assuming 'total' as in Total RNA library preparation")
        elif len(prep_files) == 1:
            with open(prep_files[0]) as file:
                prep = file.read().strip().lower()
                if prep not in ["total", "mrna"]:
                    raise ValueError(f"Prep method must be either 'total' or 'mrna' for {label}")
        elif len(prep_files) > 1:
            raise ValueError(f"Multiple matching prep files for {label}, make sure there is only one copy for each replicate in COMO_input")

        # Get fragment length
        mean_fragment_size = 100
        if len(frag_files) == 0:
            logger.warning(f"\nNo fragment file found for {label}, using '100'. This must be defined by the user in order to use zFPKM normalization")
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
            raise ValueError(f"Multiple matching fragment files for {label}, make sure there is only one copy for each replicate in COMO_input")

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


def _split_config_df(df):
    """
    Split a config dataframe into two seperate ones. One for Total RNA library prep, one for mRNA
    """
    df_t = df[df["library_prep"] == "total"]
    df_m = df[df["library_prep"] == "mrna"]

    return df_t, df_m


def _split_counts_matrices(count_matrix_all, df_total, df_mrna):
    """
    Split a counts-matrix dataframe into two seperate ones. One for Total RNA library prep, one for mRNA
    """
    logger.info(f"Reading gene count matrix file at '{count_matrix_all}'")
    matrix_all = pd.read_csv(count_matrix_all)
    matrix_total = matrix_all[["ensembl_gene_id"] + [n for n in matrix_all.columns if n in df_total["sample_name"].tolist()]]
    matrix_mrna = matrix_all[["ensembl_gene_id"] + [n for n in matrix_all.columns if n in df_mrna["sample_name"].tolist()]]

    return matrix_total, matrix_mrna


def _create_gene_info_file(matrix_file_list: list[str], input_format: Input, taxon_id, config: Config):
    """
    Create gene info file for specified context by reading first column in its count matrix file at
     results/<context name>/gene_info_<context name>.csv
    """

    logger.info("Fetching gene info")
    gene_info_file = config.data_dir / "gene_info.csv"
    genes = set()
    for file in matrix_file_list:
        df = pd.read_csv(file, low_memory=False)
        genes.update(df["ensembl_gene_id"].astype(str).tolist())
    genes = list(genes)

    # Create our output database format
    # Do not include values equal to "form"
    # Remove items not equal to `form` because the input database cannot exist as an output database
    output_db: list[Output] = [
        i for i in [Output.ENSEMBL_GENE_ID, Output.GENE_SYMBOL, Output.GENE_ID, Output.CHROMOSOMAL_LOCATION] if i.value != input_format.value
    ]

    biodbnet = BioDBNet(cache=True)
    gene_info = biodbnet.db2db(
        input_values=genes,
        input_db=input_format,
        output_db=output_db,
        taxon=taxon_id,
    )

    gene_info.rename(columns={Output.ENSEMBL_GENE_ID.value: "ensembl_gene_id"}, inplace=True)
    gene_info["start_position"] = gene_info["Chromosomal Location"].str.extract(r"chr_start: (\d+)")
    gene_info["end_position"] = gene_info["Chromosomal Location"].str.extract(r"chr_end: (\d+)")
    gene_info.rename(columns={"Gene Symbol": "hgnc_symbol", "Gene ID": "entrez_gene_id"}, inplace=True)
    gene_info.drop(["Chromosomal Location"], axis=1, inplace=True)
    gene_info.to_csv(gene_info_file, index=False)
    logger.info(f"Gene Info file written at '{gene_info_file}'")


def _handle_context_batch(
    context_names: list[str],
    mode,
    input_format: Input,
    taxon_id,
    provided_matrix_file,
    config: Config,
):
    """
    Handle iteration through each context type and create files according to flag used (config, matrix, info)
    """
    trnaseq_config_filename = config.config_dir / "trnaseq_data_inputs_auto.xlsx"
    mrnaseq_config_filename = config.config_dir / "mrnaseq_data_inputs_auto.xlsx"

    tflag = False  # turn on when any total set is found to prevent writer from being init multiple times or empty
    mflag = False  # turn on when any mrna set is found to prevent writer from being init multiple times or empty

    logger.info(f"Found {len(context_names)} contexts to process: {', '.join(context_names)}")

    tmatrix_files = []
    mmatrix_files = []
    for context_name in context_names:
        context_name = context_name.strip(" ")
        logger.info(f"Preprocessing {context_name}")
        gene_output_dir = config.result_dir / context_name
        matrix_output_dir = config.data_dir / "data_matrices" / context_name

        gene_output_dir.parent.mkdir(parents=True, exist_ok=True)
        matrix_output_dir.parent.mkdir(parents=True, exist_ok=True)

        logger.info(f"Gene info output directory is '{gene_output_dir}'")

        matrix_path_all = matrix_output_dir / f"gene_counts_matrix_full_{context_name}.csv"
        matrix_path_total = matrix_output_dir / f"gene_counts_matrix_total_{context_name}.csv"
        matrix_path_mrna = matrix_output_dir / f"gene_counts_matrix_mrna_{context_name}.csv"

        if mode == "make":
            _create_counts_matrix(context_name, config=config)
            # TODO: warn user or remove samples that are all 0 to prevent density plot error in zFPKM
            df = _create_config_df(context_name)
            df_t, df_m = _split_config_df(df)

            if not df_t.empty:
                if not tflag:
                    tflag = True
                    twriter = pd.ExcelWriter(trnaseq_config_filename)

                tmatrix_files.append(matrix_path_total)
                df_t.to_excel(twriter, sheet_name=context_name, header=True, index=False)

            if not df_m.empty:
                if not mflag:
                    mflag = True
                    mwriter = pd.ExcelWriter(mrnaseq_config_filename)

                mmatrix_files.append(matrix_path_mrna)
                df_m.to_excel(mwriter, sheet_name=context_name, header=True, index=False)

            tmatrix, mmatrix = _split_counts_matrices(matrix_path_all, df_t, df_m)
            if len(tmatrix.columns) >= 1:
                tmatrix.to_csv(matrix_path_total, header=True, index=False)
            if len(mmatrix.columns) >= 1:
                mmatrix.to_csv(matrix_path_mrna, header=True, index=False)

    if mode == "make":
        if tflag:
            twriter.close()
        if mflag:
            mwriter.close()

        _create_gene_info_file(tmatrix_files + mmatrix_files, input_format, taxon_id, config=config)

    else:
        matrix_files: list[str] = stringlist_to_list(provided_matrix_file)
        _create_gene_info_file(matrix_files, input_format, taxon_id, config=config)


def rnaseq_preprocess(
    context_names: list[str],
    mode: str,
    input_format: Input,
    taxon_id: Union[int, str],
    matrix_file: Optional[str | Path] = None,
    config: Config = None,
) -> None:
    config = Config() if config is None else config
    if mode not in ["make", "provide"]:
        raise ValueError("mode must be either 'make' or 'provide'")

    if input_format not in [Input.ENSEMBL_GENE_ID, Input.GENE_SYMBOL, Input.GENE_ID]:
        raise ValueError(f"input_format must be either 'ENSEMBL_GENE_ID', 'GENE_SYMBOL', or 'GENE_ID' (provided: {input_format})")

    if not isinstance(taxon_id, int) and taxon_id not in ["human", "mouse"]:
        raise ValueError("taxon_id must be either an integer, or accepted string ('mouse', 'human')")

    _handle_context_batch(
        context_names=context_names,
        mode=mode,
        input_format=input_format,
        taxon_id=taxon_id,
        provided_matrix_file=Path(matrix_file if matrix_file is not None else "").as_posix(),
        config=config,
    )


def _parse_args():
    """
    Parse arguments to rnaseq_preprocess.py, create a gene info files for each provided context at:
    /work/data/results/<context name>/gene_info_<context name>.csv.

     If using --info-matrix or --info-matrix-config:
    create gene count matrix file at /work/data/data_matrices/<context name>/gene_counts_matrix_<context name>.csv,

    If using --info-matrix-config:
    create config file at /work/data/config_sheets/rnaseq_data_inputs_auto.xlsx
    """

    parser = argparse.ArgumentParser(
        prog="rnaseq_preprocess.py",
        description="""
            Fetches additional gene information from a provided matrix or gene counts, or optionally creates this
            matrix using gene count files obtained using STAR aligner. Creation of counts matrix from STAR aligner
            output requires that the 'COMO_input' folder exists and is correctly structured according to the
            normalization technique being used. A correctly structured folder can be made using our Snakemake-based
            alignment pipeline at:
            https://github.com/HelikarLab/FastqToGeneCounts""",
        epilog="""
            For additional help, please post questions/issues in the MADRID GitHub repo at
            https://github.com/HelikarLab/MADRID or email babessell@gmail.com""",
    )

    parser.add_argument(
        "-n",
        "--context-names",
        type=str,
        nargs="+",
        required=True,
        dest="context_names",
        help="""Tissue/cell name of models to generate. These names should correspond to the folders
                             in 'COMO_input/' if creating count matrix files, or to
                             'work/data/data_matrices/<context name>/gene_counts_matrix_<context name>.csv' if supplying
                             the count matrix as an imported .csv file. If making multiple models in a batch, then
                             use the format: "context1 context2 context3". """,
    )

    parser.add_argument(
        "-f",
        "--gene-format",
        type=str,
        required=False,
        default="Ensembl Gene ID",
        dest="gene_format",
        help="Format of Genes, accepts 'Ensembl', 'Entrez', or'HGNC symbol'",
    )

    parser.add_argument(
        "-i",
        "--taxon-id",
        required=False,
        default=9606,
        dest="taxon_id",
        help="BioDbNet taxon ID number, also accepts 'human', or 'mouse'",
    )

    parser.add_argument(
        "--mode",
        type=str,
        required=True,
        dest="mode",
        choices=["make", "provide"],
        help="Mode of rnaseq_preprocess.py, either 'make' or 'provide'",
    )
    parser.add_argument(
        "--matrix",
        required=False,
        dest="provided_matrix_fname",
        default=None,
        help="Name of provided counts matrix in /work/data/data_matrices/<context name>/<NAME OF FILE>.csv",
    )

    parsed = parser.parse_args()
    parsed.context_names = stringlist_to_list(parsed.context_names)
    args = _Arguments(**vars(parsed))

    return args


if __name__ == "__main__":
    args: _Arguments = _parse_args()
    taxon_id_value = args.taxon_id.value if isinstance(args.taxon_id, Taxon) else args.taxon_id

    rnaseq_preprocess(
        context_names=args.context_names,
        mode=args.mode,
        input_format=args.gene_format,
        taxon_id=taxon_id_value,
        matrix_file=args.provided_matrix_fname,
    )
