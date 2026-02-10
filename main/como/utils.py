from __future__ import annotations

import contextlib
import io
import sys
from collections.abc import Iterator, Sequence
from pathlib import Path
from typing import Any, Literal, NoReturn, TextIO, TypeVar, overload

import numpy.typing as npt
import pandas as pd
import scanpy as sc
from loguru import logger

from como.data_types import LOG_FORMAT, Algorithm, LogLevel
from como.pipelines.identifier import convert

T = TypeVar("T")
__all__ = [
    "get_missing_gene_data",
    "log_and_raise_error",
    "num_columns",
    "num_rows",
    "read_file",
    "set_up_logging",
    "split_gene_expression_data",
    "stringlist_to_list",
    "suppress_stdout",
]


def stringlist_to_list(stringlist: str | list[str]) -> list[str]:
    """Convert a string from the command line into a Python list.

    In doing so, we must deprecate the use of the current method

    If '[' and ']' are present in the first and last items of the list,
        assume we are using the "old" method of providing context names

    Args:
        stringlist: The "string list" gathered from the command line. Example input: "['mat', 'xml', 'json']"

    Returns:
        A list of strings. Example output: ['mat', 'xml', 'json']

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
    identifier_column: Literal["ensembl_gene_id", "entrez_gene_id"],
    recon_algorithm: Algorithm | None = None,
    *,
    ensembl_as_index: bool = True,
) -> pd.DataFrame:
    """Split the gene expression data into single-gene and multiple-gene names.

    Arg:
        expression_data: The gene expression data to map
        identifier_column: The column containing the gene identifiers, either 'ensembl_gene_id'
        recon_algorithm: The recon algorithm used to generate the gene expression data
        ensembl_as_index: Should the 'ensembl_gene_id' column be set as

    Returns:
        A pandas DataFrame with the split gene expression data
    """
    expression_data.columns = [c.lower() for c in expression_data.columns]
    if "combine_z" in expression_data.columns:
        expression_data.rename(columns={"combine_z": "active"}, inplace=True)

    expression_data = expression_data[[identifier_column, "active"]]
    single_gene_names = expression_data[~expression_data[identifier_column].astype(str).str.contains("//")]
    multiple_gene_names = expression_data[expression_data[identifier_column].astype(str).str.contains("//")]
    split_gene_names = multiple_gene_names.assign(
        ensembl_gene_id=multiple_gene_names[identifier_column].astype(str).str.split("///")
    ).explode(identifier_column)

    gene_expressions = pd.concat([single_gene_names, split_gene_names], axis=0, ignore_index=True)
    if ensembl_as_index:
        gene_expressions.set_index(identifier_column, inplace=True)
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


def get_missing_gene_data(values: Sequence[str] | pd.DataFrame | sc.AnnData, taxon_id: int | str) -> pd.DataFrame:
    """Get missing gene data from a given set of values.

    This function will attempt to find gene identifiers in the provided values and convert them to a DataFrame
    containing "entrez_gene_id", "ensembl_gene_id", and "gene_symbol"

    :param values: The values to extract gene identifiers from.
        This can be a list of strings, a pandas DataFrame, or a scanpy AnnData object.
    :param taxon_id: The taxonomic identifier to use for gene consversion
    :return: A DataFrame containing "entrez_gene_id", "ensembl_gene_id", and "gene_symbol" for the provided values
    """
    # second isinstance required for static type check to be happy
    # if isinstance(values, list) and not isinstance(values, pd.DataFrame):
    if isinstance(values, list):
        return convert(ids=values, taxon=taxon_id)
    elif isinstance(values, pd.DataFrame):
        # raise error if duplicate column names exist
        if any(values.columns.duplicated(keep=False)):
            duplicate_cols = values.columns[values.columns.duplicated(keep=False)].unique().tolist()
            log_and_raise_error(
                message=(
                    f"Duplicate column names exist! This will result in an error processing data. "
                    f"Duplicates: {','.join(duplicate_cols)}"
                ),
                error=ValueError,
                level=LogLevel.CRITICAL,
            )

        names: list[str] = values.columns.tolist()
        if values.index.name is not None:
            names.append(str(values.index.name))
        if "gene_symbol" in names:
            return get_missing_gene_data(
                values["gene_symbol"].tolist() if "gene_symbol" in values.columns else values.index.tolist(),
                taxon_id=taxon_id,
            )
        elif "entrez_gene_id" in names:
            return get_missing_gene_data(
                values["entrez_gene_id"].tolist() if "entrez_gene_id" in values.columns else values.index.tolist(),
                taxon_id=taxon_id,
            )
        elif "ensembl_gene_id" in names:
            return get_missing_gene_data(
                values["ensembl_gene_id"].tolist() if "ensembl_gene_id" in values.columns else values.index.tolist(),
                taxon_id=taxon_id,
            )
        else:
            log_and_raise_error(
                message="Unable to find 'gene_symbol', 'entrez_gene_id', or 'ensembl_gene_id' in the input matrix.",
                error=ValueError,
                level=LogLevel.CRITICAL,
            )
    else:
        log_and_raise_error(
            message=f"Values must be a list of strings or a pandas DataFrame, got: {type(values)}",
            error=TypeError,
            level=LogLevel.CRITICAL,
        )


@overload
def read_file(path: None, h5ad_as_df: bool = True, **kwargs: Any) -> None: ...


@overload
def read_file(path: pd.DataFrame, h5ad_as_df: bool = True, **kwargs: Any) -> pd.DataFrame: ...


@overload
def read_file(path: io.StringIO, h5ad_as_df: bool = True, **kwargs: Any) -> pd.DataFrame: ...


@overload
def read_file(path: sc.AnnData, h5ad_as_df: Literal[False], **kwargs: Any) -> sc.AnnData: ...


@overload
def read_file(path: sc.AnnData, h5ad_as_df: Literal[True] = True, **kwargs: Any) -> pd.DataFrame: ...


@overload
def read_file(path: Path, h5ad_as_df: Literal[False], **kwargs: Any) -> pd.DataFrame | sc.AnnData: ...


@overload
def read_file(path: Path, h5ad_as_df: Literal[True] = True, **kwargs: Any) -> pd.DataFrame: ...


def read_file(  # noqa: C901
    path: Path | io.StringIO | pd.DataFrame | sc.AnnData | None,
    h5ad_as_df: bool = True,
    **kwargs: Any,
) -> pd.DataFrame | sc.AnnData | None:
    """Read a filepath and return pandas.DataFrame or scanpy.AnnData.

    If the provided path is None, None will also be returned.
    None may be provided to this function so that `asyncio.gather` can safely be used on all sources
        (trna, mrna, scrna, proteomics) without needing to check if the user has provided those sources

    Args:
        path: The path to read from
        h5ad_as_df: If True and the file is an h5ad, return a pandas DataFrame of the .X matrix
        kwargs: Additional arguments to pass to pandas.read_csv, pandas.read_excel, or scanpy.read_h5ad

    Returns:
        None, or a pandas DataFrame or AnnData

    """
    if isinstance(path, pd.DataFrame):
        return path
    elif isinstance(path, sc.AnnData):
        return path.to_df().T if h5ad_as_df else path
    elif isinstance(path, io.StringIO):
        return pd.read_csv(path, **kwargs)
    elif not path:
        return None

    if isinstance(path, Path) and not path.exists():
        log_and_raise_error(f"File {path} does not exist", error=FileNotFoundError, level=LogLevel.CRITICAL)

    match path.suffix:
        case ".csv" | ".tsv" | ".txt" | ".tab" | ".sf":
            kwargs.setdefault("sep", "," if path.suffix == ".csv" else "\t")  # set sep if not defined
            return pd.read_csv(path, **kwargs)
        case ".xlsx" | ".xls":
            return pd.read_excel(path, **kwargs)
        case ".h5ad":
            adata: sc.AnnData = sc.read_h5ad(path, **kwargs)
            if h5ad_as_df:
                df = adata.to_df().T
                if not df.index.name:
                    df.index.name = "gene_symbol"
                df.reset_index(inplace=True)
                return df
            return adata
        case _:
            log_and_raise_error(
                f"Unknown file extension '{path.suffix}'. "
                "Valid options are '.tsv', '.csv', '.xlsx', '.xls', or '.h5ad'",
                error=ValueError,
                level=LogLevel.CRITICAL,
            )


@overload
def _listify(value: list[T]) -> list[T]: ...


@overload
def _listify(value: T) -> list[T]: ...


def _listify(value: T | list[T]) -> list[T]:
    """Convert items into a list.

    Args:
        value: The value to convert to a list (if it isn't already)

    Returns:
        A list of the provided value
    """
    return value if isinstance(value, list) else [value]


def num_rows(item: pd.DataFrame | npt.NDArray) -> int:
    """Return the number of rows in an object.

    :param item: The object to check the number of rows of
    :return: The number of rows in the provided object
    """
    return item.shape[0]


def num_columns(item: pd.DataFrame | npt.NDArray) -> int:
    """Return the number of columns in an object.

    :param item: The object to check the number of columns of
    :return: The number of columns in the provided object
    """
    return item.shape[1]


def return_placeholder_data() -> pd.DataFrame:
    return pd.DataFrame(data=0, index=pd.Index(data=[0], name="entrez_gene_id"), columns=["expressed", "top"])


def set_up_logging(
    level: LogLevel | str,
    location: str | TextIO,
    formatting: str = LOG_FORMAT,
):
    """Set up logging for the application.

    :param level: The default logging level to use (e.g., LogLevel.INFO, LogLevel.DEBUG, etc.)
    :param location: The location to log to (e.g., a file path or sys.stdout)
    :param formatting: The log message format to use (default is LOG_FORMAT)
    """
    if isinstance(level, str):
        level = LogLevel[level.upper()]
    with contextlib.suppress(ValueError):
        logger.remove(0)
        logger.add(sink=location, level=level.value, format=formatting)


def log_and_raise_error(
    message: str,
    *,
    error: type[BaseException],
    level: LogLevel,
) -> NoReturn:
    """Log an error message and raise an exception.

    :param message: The error message to log and include in the raised exception
    :param error: The type of exception to raise (e.g., ValueError, File NotFoundError, etc.)
    :param level: The LogLevel at which to log the error message (e.g., LogLevel.ERROR, LogLevel.CRITICAL)
    """
    caller = logger.opt(depth=1)
    if level == LogLevel.ERROR:
        caller.error(message)
        raise error(message)
    if level == LogLevel.CRITICAL:
        caller.critical(message)
        raise error(message)

    raise ValueError(f"When raising an error, LogLevel.ERROR or LogLevel.CRITICAL must be used. Got: {level}")
