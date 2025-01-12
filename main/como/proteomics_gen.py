from __future__ import annotations

import sys
from io import TextIOWrapper
from pathlib import Path

import numpy as np
import pandas as pd
from fast_bioservices.biodbnet import BioDBNet, Input, Output
from loguru import logger

from como.data_types import LOG_FORMAT, LogLevel
from como.project import Config
from como.proteomics_preprocessing import protein_transform_main
from como.utils import _log_and_raise_error, _set_up_logging, return_placeholder_data


# Load Proteomics
def process_proteomics_data(path: Path) -> pd.DataFrame:
    """Load proteomics data from a given context and filename."""
    # Preprocess data, drop na, duplicate ';' in symbol,
    matrix: pd.DataFrame = pd.read_csv(path)
    if "gene_symbol" not in matrix.columns:
        _log_and_raise_error(
            "No gene_symbol column found in proteomics data.",
            error=ValueError,
            level=LogLevel.ERROR,
        )

    matrix["gene_symbol"] = matrix["gene_symbol"].astype(str)
    matrix.dropna(subset=["gene_symbol"], inplace=True)
    matrix = matrix.assign(gene_symbol=matrix["gene_symbol"].str.split(";")).explode("gene_symbol")
    return matrix


# read map to convert to entrez
async def load_gene_symbol_map(gene_symbols: list[str], entrez_map: Path | None = None):
    """Add descirption...."""
    if entrez_map and entrez_map.exists():
        df = pd.read_csv(entrez_map, index_col="gene_symbol")
    else:
        biodbnet = BioDBNet()
        df = await biodbnet.async_db2db(
            values=gene_symbols,
            input_db=Input.GENE_SYMBOL,
            output_db=[Output.GENE_ID, Output.ENSEMBL_GENE_ID],
        )
        df.loc[df["gene_id"] == "-", ["gene_id"]] = np.nan
        df.to_csv(entrez_map, index_label="gene_symbol")

    return df[~df.index.duplicated()]


def abundance_to_bool_group(
    context_name,
    abundance_filepath: Path,
    output_gaussian_png_filepath: Path,
    output_gaussian_html_filepath: Path,
    output_z_score_matrix_filepath: Path,
    abundance_matrix: pd.DataFrame,
    replicate_ratio: float,
    high_confidence_replicate_ratio: float,
    quantile: float,
    output_boolean_filepath: Path,
):
    """Convert proteomic data to boolean expression."""
    abundance_matrix.to_csv(abundance_filepath, index_label="entrez_gene_id")
    protein_transform_main(
        abundance_df=abundance_matrix,
        output_gaussian_png_filepath=output_gaussian_png_filepath,
        output_gaussian_html_filepath=output_gaussian_html_filepath,
        output_z_score_matrix_filepath=output_z_score_matrix_filepath,
    )

    # Logical Calculation
    abundance_matrix_nozero = abundance_matrix.replace(0, np.nan)
    thresholds = abundance_matrix_nozero.quantile(quantile, axis=0)
    testbool = pd.DataFrame(0, columns=abundance_matrix.columns, index=abundance_matrix.index)

    for col in abundance_matrix.columns:
        testbool.loc[abundance_matrix[col] > thresholds[col], [col]] = 1

    abundance_matrix["expressed"] = 0
    abundance_matrix["high"] = 0
    abundance_matrix["pos"] = abundance_matrix[abundance_matrix > 0].sum(axis=1) / abundance_matrix.count(axis=1)
    abundance_matrix.loc[(abundance_matrix["pos"] >= replicate_ratio), ["expressed"]] = 1
    abundance_matrix.loc[(abundance_matrix["pos"] >= high_confidence_replicate_ratio), ["high"]] = 1

    abundance_matrix.to_csv(output_boolean_filepath, index_label="entrez_gene_id")


def to_bool_context(context_name, group_ratio, hi_group_ratio, group_names):
    """Convert proteomic data to boolean expression."""
    config = Config()
    output_dir = config.result_dir / context_name / "proteomics"
    merged_df = pd.DataFrame(columns=["entrez_gene_id", "expressed", "high"])
    merged_df.set_index(["entrez_gene_id"], inplace=True)
    merged_hi_df = merged_df

    for group in group_names:
        read_filepath = output_dir / f"bool_prot_Matrix_{context_name}_{group}.csv"
        read_df = pd.read_csv(read_filepath)
        read_df.set_index("entrez_gene_id", inplace=True)
        read_df = read_df[["expressed", "high"]]

        if not merged_df.empty:
            merged_df = pd.merge(merged_df, read_df["expressed"], right_index=True, left_index=True)
            merged_hi_df = pd.merge(merged_hi_df, read_df["high"], right_index=True, left_index=True)

        else:
            merged_df = read_df["expressed"].to_frame()
            merged_hi_df = read_df["high"].to_frame()

    if len(merged_df.columns) > 1:
        merged_df.apply(lambda x: sum(x) / len(merged_df.columns) >= group_ratio, axis=1, result_type="reduce")
        merged_hi_df.apply(lambda x: sum(x) / len(merged_hi_df.columns) >= hi_group_ratio, axis=1, result_type="reduce")

    out_df = pd.merge(merged_df, merged_hi_df, right_index=True, left_index=True)
    out_filepath = output_dir / f"Proteomics_{context_name}.csv"
    out_df.to_csv(out_filepath, index_label="entrez_gene_id")
    logger.success(f"Test Data Saved to {out_filepath}")


# read data from csv files
def load_proteomics_tests(filename, context_name):
    """Load statistical test results."""
    config = Config()

    def load_empty_dict():
        return "dummy", return_placeholder_data()

    if not filename or filename == "None":  # if not using proteomics load empty dummy data matrix
        return load_empty_dict()

    inquiry_full_path = config.data_dir / "config_sheets" / filename
    if not inquiry_full_path.exists():
        _log_and_raise_error(
            f"Error: file not found {inquiry_full_path}", error=FileNotFoundError, level=LogLevel.ERROR
        )

    filename = f"Proteomics_{context_name}.csv"
    full_save_filepath = config.result_dir / context_name / "proteomics" / filename
    if full_save_filepath.exists():
        data = pd.read_csv(full_save_filepath, index_col="entrez_gene_id")
        logger.success(f"Read from {full_save_filepath}")
        return context_name, data

    else:
        logger.warning(
            f"Proteomics gene expression file for {context_name} was not found at {full_save_filepath}. "
            f"Is this intentional?"
        )
        return load_empty_dict()


async def proteomics_gen(
    context_name: str,
    config_filepath: Path,
    matrix_filepath: Path,
    output_boolean_filepath: Path,
    output_z_score_matrix_filepath: Path,
    output_gaussian_png_filepath: Path | None = None,
    output_gaussian_html_filepath: Path | None = None,
    input_entrez_map: Path | None = None,
    replicate_ratio: float = 0.5,
    batch_ratio: float = 0.5,
    high_confidence_replicate_ratio: float = 0.7,
    high_confidence_batch_ratio: float = 0.7,
    quantile: int = 25,
):
    """Generate proteomics data."""
    if not config_filepath.exists():
        _log_and_raise_error(
            f"Config file not found at {config_filepath}",
            error=FileNotFoundError,
            level=LogLevel.ERROR,
        )
    if config_filepath.suffix not in (".xlsx", ".xls"):
        _log_and_raise_error(
            f"Config file must be an xlsx or xls file at {config_filepath}",
            error=FileNotFoundError,
            level=LogLevel.ERROR,
        )

    if not matrix_filepath.exists():
        _log_and_raise_error(
            f"Matrix file not found at {matrix_filepath}",
            error=FileNotFoundError,
            level=LogLevel.ERROR,
        )
    if matrix_filepath.suffix not in {".csv"}:
        _log_and_raise_error(
            f"Matrix file must be a csv file at {matrix_filepath}",
            error=FileNotFoundError,
            level=LogLevel.ERROR,
        )

    if quantile < 0 or quantile > 100:
        _log_and_raise_error(
            "Quantile must be an integer from 0 to 100",
            error=ValueError,
            level=LogLevel.ERROR,
        )
    quantile /= 100

    config_df = pd.read_excel(config_filepath, sheet_name=context_name)
    matrix: pd.DataFrame = process_proteomics_data(matrix_filepath)

    groups = config_df["group"].unique().tolist()

    for group in groups:
        indices = np.where([g == group for g in config_df["group"]])
        sample_columns = [*np.take(config_df["sample_name"].to_numpy(), indices).ravel().tolist(), "gene_symbol"]
        matrix = matrix.loc[:, sample_columns]

        symbols_to_gene_ids = await load_gene_symbol_map(
            gene_symbols=matrix["gene_symbol"].tolist(),
            entrez_map=input_entrez_map,
        )
        matrix.dropna(subset=["gene_symbol"], inplace=True)
        if "uniprot" in matrix.columns:
            matrix.drop(columns=["uniprot"], inplace=True)

        matrix = matrix.groupby(["gene_symbol"]).agg("max")
        matrix["entrez_gene_id"] = symbols_to_gene_ids["gene_id"]
        matrix.dropna(subset=["entrez_gene_id"], inplace=True)
        matrix.set_index("entrez_gene_id", inplace=True)

        # bool_filepath = output_dir / f"bool_prot_Matrix_{context_name}_{group_name}.csv"
        abundance_to_bool_group(
            context_name=context_name,
            abundance_filepath=matrix_filepath,
            abundance_matrix=matrix,
            replicate_ratio=replicate_ratio,
            high_confidence_replicate_ratio=high_confidence_replicate_ratio,
            quantile=quantile,
            output_boolean_filepath=output_boolean_filepath,
            output_gaussian_png_filepath=output_gaussian_png_filepath,
            output_gaussian_html_filepath=output_gaussian_html_filepath,
            output_z_score_matrix_filepath=output_z_score_matrix_filepath,
        )
    to_bool_context(
        context_name=context_name,
        group_ratio=batch_ratio,
        hi_group_ratio=high_confidence_batch_ratio,
        group_names=groups,
    )
