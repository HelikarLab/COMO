from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from fast_bioservices.biodbnet import BioDBNet, Input, Output
from loguru import logger

from como import return_placeholder_data
from como.project import Config
from como.proteomics_preprocessing import protein_transform_main


# Load Proteomics
def process_proteomics_data(path: Path) -> pd.DataFrame:
    """Load proteomics data from a given context and filename."""
    # Preprocess data, drop na, duplicate ';' in symbol,
    matrix: pd.DataFrame = pd.read_csv(path)
    if "gene_symbol" not in matrix.columns:
        raise ValueError("No gene_symbol column found in proteomics data.")

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
        df.to_csv(filepath, index_label="gene_symbol")

    return df[~df.index.duplicated()]


def abundance_to_bool_group(
    context_name,
    group_name,
    abundance_matrix,
    replicate_ratio,
    high_confidence_replicate_ratio,
    quantile,
    output_boolean_filepath: Path,
):
    """Descrioption...."""
    config = Config()
    output_dir = config.result_dir / context_name / "proteomics"
    output_dir.mkdir(parents=True, exist_ok=True)

    # write group abundances to individual files
    abundance_filepath = (
        config.result_dir / context_name / "proteomics" / "".join(["protein_abundance_", group_name, ".csv"])
    )
    abundance_matrix.to_csv(abundance_filepath, index_label="entrez_gene_id")
    protein_transform_main(abundance_matrix, output_dir, group_name)

    # Logical Calculation
    abundance_matrix_nozero = abundance_matrix.replace(0, np.nan)
    thresholds = abundance_matrix_nozero.quantile(quantile, axis=0)
    testbool = pd.DataFrame(0, columns=list(abundance_matrix), index=abundance_matrix.index)

    for col in list(abundance_matrix):
        testbool.loc[abundance_matrix[col] > thresholds[col], [col]] = 1

    abundance_matrix["pos"] = (abundance_matrix > 0).sum(axis=1) / abundance_matrix.count(axis=1)
    abundance_matrix["expressed"] = 0
    abundance_matrix.loc[(abundance_matrix["pos"] >= rep_ratio), ["expressed"]] = 1
    abundance_matrix["high"] = 0
    abundance_matrix.loc[(abundance_matrix["pos"] >= hi_rep_ratio), ["high"]] = 1

    bool_filepath = output_dir / f"bool_prot_Matrix_{context_name}_{group_name}.csv"
    abundance_matrix.to_csv(bool_filepath, index_label="entrez_gene_id")


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
        raise FileNotFoundError(f"Error: file not found {inquiry_full_path}")

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
    input_entrez_map: Path | None = None,
    replicate_ratio: float = 0.5,
    batch_ratio: float = 0.5,
    high_confidence_replicate_ratio: float = 0.7,
    high_confience_batch_ratio: float = 0.7,
    quantile: int = 25,
):
    """Generate proteomics data."""
    if not config_filepath.exists():
        raise FileNotFoundError(f"Config file not found at {config_filepath}")
    if config_filepath.suffix not in (".xlsx", ".xls"):
        raise FileNotFoundError(f"Config file must be an xlsx or xls file at {config_filepath}")

    if not matrix_filepath.exists():
        raise FileNotFoundError(f"Matrix file not found at {matrix_filepath}")
    if matrix_filepath.suffix not in {".csv"}:
        raise FileNotFoundError(f"Matrix file must be a csv file at {matrix_filepath}")

    if quantile < 0 or quantile > 100:
        raise ValueError("Quantile must be an integer from 0 to 100")
    quantile /= 100

    logger.info(f"Config file is at '{prot_config_filepath}'")

    xl = pd.ExcelFile(prot_config_filepath)
    sheet_names = xl.sheet_names

    for context_name in sheet_names:
        datafilename = "".join(["protein_abundance_", context_name, ".csv"])
        config_sheet = pd.read_excel(prot_config_filepath, sheet_name=context_name)
        groups = config_sheet["group"].unique().tolist()

        for group in groups:
            group_idx = np.where([g == group for g in config_sheet["group"].tolist()])
            cols = [*np.take(config_sheet["sample_name"].to_numpy(), group_idx).ravel().tolist(), "gene_symbol"]

            proteomics_data = load_proteomics_data(datafilename, context_name)
            proteomics_data = proteomics_data.loc[:, cols]

            symbols_to_ids = await load_gene_symbol_map(gene_symbols=proteomics_data["gene_symbol"].tolist())
            proteomics_data.dropna(subset=["gene_symbol"], inplace=True)
            if "uniprot" in proteomics_data.columns:
                proteomics_data.drop(columns=["uniprot"], inplace=True)

            proteomics_data = proteomics_data.groupby(["gene_symbol"]).agg("max")
            proteomics_data["entrez_gene_id"] = symbols_to_ids["gene_id"]
            proteomics_data.dropna(subset=["entrez_gene_id"], inplace=True)
            proteomics_data.set_index("entrez_gene_id", inplace=True)

            # save proteomics data by test
            abundance_to_bool_group(context_name, group, proteomics_data, rep_ratio, hi_rep_ratio, quantile)
        to_bool_context(context_name, group_ratio, hi_group_ratio, groups)


        )
    )
