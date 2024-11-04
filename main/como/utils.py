import contextlib
import io
import sre_constants
import sys
from pathlib import Path
from typing import Iterator

import aiofiles
import pandas as pd
import scanpy as sc
from fast_bioservices import BioDBNet, Output, Taxon

__all__ = ["Compartments", "stringlist_to_list", "split_gene_expression_data", "suppress_stdout"]


class Compartments:
    """
    This enum will be used to convert from compartment "long-hand" to "short-hand"

    Shorthand from: https://cobrapy.readthedocs.io/en/latest/_modules/cobra/medium/annotations.html

    "Extracellular" -> "e"
    "golgi" -> "g"
    """

    SHORTHAND = {
        "ce": ["cell envelope"],
        "c": [
            "cytoplasm",
            "cytosol",
            "default",
            "in",
            "intra cellular",
            "intracellular",
            "intracellular region",
            "intracellular space",
        ],
        "er": ["endoplasmic reticulum"],
        "erm": ["endoplasmic reticulum membrane"],
        "e": [
            "extracellular",
            "extraorganism",
            "out",
            "extracellular space",
            "extra organism",
            "extra cellular",
            "extra-organism",
            "external",
            "external medium",
        ],
        "f": ["flagellum", "bacterial-type flagellum"],
        "g": ["golgi", "golgi apparatus"],
        "gm": ["golgi membrane"],
        "h": ["chloroplast"],
        "l": ["lysosome"],
        "im": ["mitochondrial intermembrane space"],
        "mm": ["mitochondrial membrane"],
        "m": ["mitochondrion", "mitochondria"],
        "n": ["nucleus"],
        "p": ["periplasm", "periplasmic space"],
        "x": ["peroxisome", "glyoxysome"],
        "u": ["thylakoid"],
        "vm": ["vacuolar membrane"],
        "v": ["vacuole"],
        "w": ["cell wall"],
        "s": ["eyespot", "eyespot apparatus", "stigma"],
    }

    _REVERSE_LOOKUP = {value.lower(): key for key, values in SHORTHAND.items() for value in values}

    @classmethod
    def get(cls, longhand: str) -> str:
        return cls._REVERSE_LOOKUP.get(longhand.lower(), None)


def stringlist_to_list(stringlist: str | list[str]) -> list[str]:
    """
    We are attempting to move to a new method of gathering a list of items from the command line
    In doing so, we must deprecate the use of the current method

    If '[' and ']' are present in the first and last items of the list, assume we are using the "old" method of providing context names

    :param stringlist: The "string list" gathered from the command line. Example input: "['mat', 'xml', 'json']"
    """
    if isinstance(stringlist, str):
        if stringlist.startswith("[") and stringlist.endswith("]"):
            # Remove any brackets from the first and last items; replace quotation marks and commas with nothing
            new_list: list[str] = stringlist.strip("[]").replace("'", "").replace(" ", "").split(",")

            # Show a warning if more than one item is present in the list (this means we are using the old method)
            print("DeprecationWarning: Please use the new method of providing context names, i.e. --output-filetypes 'type1 type2 type3'.")
            print(
                "If you are using COMO, this can be done by setting the 'context_names' variable to a simple string separated by spaces. Here are a few examples!"
            )
            print("context_names = 'cellType1 cellType2 cellType3'")
            print("output_filetypes = 'output1 output2 output3'")
            print("\nYour current method of passing context names will be removed in the future. Please update your variables above accordingly!\n\n")

        else:
            new_list: list[str] = stringlist.split(" ")

        return new_list
    return stringlist


def split_gene_expression_data(expression_data: pd.DataFrame, recon_algorithm: str = "GIMME"):
    """
    Splits genes that have mapped to multiple Entrez IDs are formated as "gene12//gene2//gene3"
    """
    expression_data.columns = [c.lower() for c in expression_data.columns]
    if recon_algorithm in {"IMAT", "TINIT"}:
        expression_data.rename(columns={"combine_z": "active"}, inplace=True)

    expression_data = expression_data[["entrez_gene_id", "active"]]
    expression_data["entrez_gene_id"] = expression_data["entrez_gene_id"].astype(str)
    single_gene_names = expression_data[~expression_data["entrez_gene_id"].str.contains("//")]
    multiple_gene_names = expression_data[expression_data["entrez_gene_id"].str.contains("//")]
    split_gene_names = multiple_gene_names.assign(entrez_gene_id=multiple_gene_names["entrez_gene_id"].str.split("///")).explode("entrez_gene_id")

    gene_expressions = pd.concat([single_gene_names, split_gene_names], axis=0, ignore_index=True)
    gene_expressions.set_index("entrez_gene_id", inplace=True)
    return gene_expressions


@contextlib.contextmanager
def suppress_stdout() -> Iterator[None]:
    """
    Suppresses stdout output from the current context

    :return: The context manager
    """
    with io.StringIO() as buffer:
        try:
            sys.stdout = buffer
            yield
        finally:
            sys.stdout = sys.__stdout__


async def _format_cohersion(biodbnet: BioDBNet, *, requested_output: Output | list[Output], input_values: list[str], taxon: Taxon) -> pd.DataFrame:
    requested_output = [requested_output] if isinstance(requested_output, Output) else requested_output
    cohersion = (await biodbnet.db_find(values=input_values, output_db=requested_output, taxon=taxon)).drop(columns=["Input Type"])
    cohersion.columns = pd.Index(["input_value", *[o.value.replace(" ", "_").lower() for o in requested_output]])
    return cohersion


async def _async_read_csv(path: Path, **kwargs) -> pd.DataFrame:
    """
    Asynchronously reads a CSV file and returns a pandas DataFrame.

    :param path: The path to read from
    :param kwargs: Additional arguments to pass to pandas.read_csv
    :return: A pandas DataFrame
    """
    if not path.exists():
        raise FileNotFoundError(f"File {path} does not exist")
    if path.suffix not in {".csv", ".tsv"}:
        raise ValueError(f"File {path} is not a CSV file")

    kwargs.setdefault("sep", "," if path.suffix == ".csv" else "\t")
    async with aiofiles.open(path, "r") as f:
        content = await f.read()
        return pd.read_csv(io.StringIO(content), **kwargs)


async def _async_read_excel(path: Path, **kwargs) -> pd.DataFrame:
    """
    Asynchronously reads an Excel file and returns a pandas DataFrame.

    :param path: The path to read from
    :param kwargs: Additional arguments to pass to pandas.read_excel
    :return: A pandas DataFrame
    """
    if not path.exists():
        raise FileNotFoundError(f"File {path} does not exist")
    if path.suffix not in {".xls", ".xlsx"}:
        raise ValueError(f"File {path} is not an Excel file")

    async with aiofiles.open(path, "rb") as f:
        content = await f.read()
        return pd.read_excel(io.StringIO(content.decode()), **kwargs)


def is_notebook() -> bool:
    """
    Checks if the current environment is a Jupyter Notebook.

    :returns: True if the current environment is a Jupyter Notebook, False otherwise.
    """
    try:
        from IPython import get_ipython

        return get_ipython() is not None
    except ModuleNotFoundError:
        return False
