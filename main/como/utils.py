from __future__ import annotations

import asyncio
import contextlib
import functools
import inspect
import io
import sys
from collections.abc import Iterator
from concurrent.futures import ThreadPoolExecutor
from io import TextIOWrapper
from pathlib import Path
from typing import Union

import aiofiles
import numpy.typing as npt
import pandas as pd
import scanpy as sc
from fast_bioservices import BioDBNet, Output, Taxon
from fast_bioservices.pipeline import (
    determine_gene_type,
    ensembl_to_gene_id_and_symbol,
    gene_id_to_ensembl_and_gene_symbol,
    gene_symbol_to_ensembl_and_gene_id,
)
from loguru import logger

from como.data_types import LOG_FORMAT, Algorithm, LogLevel

__all__ = ["split_gene_expression_data", "stringlist_to_list", "suppress_stdout"]


def stringlist_to_list(stringlist: str | list[str]) -> list[str]:
    """Convert a string from the command line into a Python list.

    In doing so, we must deprecate the use of the current method

    If '[' and ']' are present in the first and last items of the list,
        assume we are using the "old" method of providing context names

    :param stringlist: The "string list" gathered from the command line. Example input: "['mat', 'xml', 'json']"
    """
    if isinstance(stringlist, list):
        return stringlist

    if not (stringlist.startswith("[") and stringlist.endswith("]")):
        return stringlist.split(" ")

    # Remove any brackets from the first and last items; replace quotation marks and commas with nothing
    new_list: list[str] = stringlist.strip("[]").replace("'", "").replace(" ", "").split(",")

    # Show a warning if more than one item is present in the list (this means we are using the old method)
    logger.critical(
        "DeprecationWarning: Please use the new method of providing context names, "
        "i.e. --output-filetypes 'type1 type2 type3'."
    )
    logger.critical(
        "If you are using COMO, this can be done by setting the 'context_names' variable to a "
        "simple string separated by spaces. Here are a few examples!"
    )
    logger.critical("context_names = 'cellType1 cellType2 cellType3'")
    logger.critical("output_filetypes = 'output1 output2 output3'")
    logger.critical(
        "\nYour current method of passing context names will be removed in the future. "
        "Update your variables above accordingly!\n\n"
    )

    return new_list


def split_gene_expression_data(
    expression_data: pd.DataFrame,
    recon_algorithm: Algorithm | None = None,
    entrez_as_index: bool = True,
):
    """Split the gene expression data into single-gene and multiple-gene names.

    :param expression_data: The gene expression data to map
    :param recon_algorithm: The recon algorithm used to generate the gene expression data
    :param entrez_as_index: Should the 'entrez_gene_id' column be set as the index
    :return:
    """
    expression_data.columns = [c.lower() for c in expression_data.columns]
    if recon_algorithm in {Algorithm.IMAT, Algorithm.TINIT}:
        expression_data.rename(columns={"combine_z": "active"}, inplace=True)

    expression_data = expression_data[["entrez_gene_id", "active"]]
    single_gene_names = expression_data[~expression_data["entrez_gene_id"].astype(str).str.contains("//")]
    multiple_gene_names = expression_data[expression_data["entrez_gene_id"].astype(str).str.contains("//")]
    split_gene_names = multiple_gene_names.assign(
        entrez_gene_id=multiple_gene_names["entrez_gene_id"].astype(str).str.split("///")
    ).explode("entrez_gene_id")

    gene_expressions = pd.concat([single_gene_names, split_gene_names], axis=0, ignore_index=True)
    if entrez_as_index:
        gene_expressions.set_index("entrez_gene_id", inplace=True)
    return gene_expressions


@contextlib.contextmanager
def suppress_stdout() -> Iterator[None]:
    """Suppress stdout output from the current context.

    :return: The context manager
    """
    with io.StringIO() as buffer:
        try:
            sys.stdout = buffer
            yield
        finally:
            sys.stdout = sys.__stdout__


async def _format_determination(
    biodbnet: BioDBNet, *, requested_output: Output | list[Output], input_values: list[str], taxon: Taxon
) -> pd.DataFrame:
    """Determine the data type of the given input values (i.e., Entrez Gene ID, Gene Symbol, etc.).

    :param biodbnet: The BioDBNet to use for deter
    :param requested_output: The data type to generate (of type `Output`)
    :param input_values: The input values to determine
    :param taxon: The Taxon ID
    :return: A pandas DataFrame
    """
    requested_output = [requested_output] if isinstance(requested_output, Output) else requested_output
    cohersion = (await biodbnet.db_find(values=input_values, output_db=requested_output, taxon=taxon)).drop(
        columns=["Input Type"]
    )
    cohersion.columns = pd.Index(["input_value", *[o.value.replace(" ", "_").lower() for o in requested_output]])
    return cohersion


async def _read_file(
    path: Path | io.StringIO | None,
    h5ad_as_df: bool = True,
    **kwargs,
) -> pd.DataFrame | sc.AnnData | None:
    """Asynchronously read a filepath and return a pandas DataFrame.

    If the provided path is None, None will also be returned.
    None may be provided to this function so that `asyncio.gather` can safely be used on all sources
        (trna, mrna, scrna, proteomics) without needing to check if the user has provided those sources

    :param path: The path to read from
    :param kwargs: Additional arguments to pass to pandas.read_csv, pandas.read_excel,
        or scanpy.read_h5ad, depending on the filepath provided
    :return: None, or a pandas DataFrame or AnnData
    """
    if not path:
        return None

    if not path.exists():
        logger.critical(f"File {path} does not exist")
        raise FileNotFoundError(f"File {path} does not exist")

def is_notebook() -> bool:
    """Check if the current environment is a Jupyter Notebook.

    :returns: True if the current environment is a Jupyter Notebook, False otherwise.
    """
    try:
        from IPython import get_ipython

        return get_ipython() is not None
    except ModuleNotFoundError:
        return False


async def convert_gene_data(values: list[str], taxon_id: int | str | Taxon) -> pd.DataFrame:
    gene_type = await determine_gene_type(values)
    if all(v == "gene_symbol" for v in gene_type.values()):
        return await gene_symbol_to_ensembl_and_gene_id(values, taxon=taxon_id)
    elif all(v == "ensembl_gene_id" for v in gene_type.values()):
        return await ensembl_to_gene_id_and_symbol(ids=values, taxon=taxon_id)
    elif all(v == "entrez_gene_id" for v in gene_type.values()):
        return await gene_id_to_ensembl_and_gene_symbol(ids=values, taxon=taxon_id)
    else:
        raise ValueError("Gene data must be of the same type (i.e., all Ensembl, Entrez, or Gene Symbols)")


def _listify(value):
    """Convert items into a list."""
    return [value] if not isinstance(value, list) else value


def _num_rows(item: pd.DataFrame | npt.NDArray) -> int:
    return item.shape[0]


def _num_columns(item: pd.DataFrame | npt.NDArray) -> int:
    return item.shape[1]


def return_placeholder_data() -> pd.DataFrame:
    return pd.DataFrame(data=0, index=pd.Index(data=[0], name="entrez_gene_id"), columns=["expressed", "top"])


def _set_up_logging(
    level: LogLevel,
    location: str | TextIOWrapper,
    formatting: str = LOG_FORMAT,
):
    with contextlib.suppress(ValueError):
        logger.remove(0)
        logger.add(sink=location, level=level.value, format=formatting)


def _log_and_raise_error(
    message: str,
    *,
    error: type[BaseException],
    level: LogLevel,
) -> None:
    caller = logger.opt(depth=1)
    match level:
        case LogLevel.ERROR:
            caller.error(message)
        case LogLevel.CRITICAL:
            caller.critical(message)
        case _:
            raise ValueError(f"When raising an error, LogLevel.ERROR or LogLevel.CRITICAL must be used. Got: {level}")

    raise error(message)
