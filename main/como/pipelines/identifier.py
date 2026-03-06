from collections.abc import Iterable, Iterator, Sequence
from typing import Any, Literal, overload

import pandas as pd
from bioservices.mygeneinfo import MyGeneInfo
from tqdm import tqdm

__all__ = [
    "build_gene_info",
    "determine_gene_type",
    "get_remaining_identifiers",
]

T_IDS = int | str | Iterable[int] | Iterable[str] | Iterable[int | str]
T_MG_SCOPE = Literal["entrezgene", "ensembl.gene", "symbol"]
T_MG_TRANSLATE = Literal["entrez_gene_id", "ensembl_gene_id", "gene_symbol"]
T_MG_RETURN = list[dict[T_MG_TRANSLATE, str]]


def _get_conversion(info: MyGeneInfo, values: T_IDS, taxon: str | int) -> list[dict[str, Any]]:
    value_list = sorted(map(str, [values] if isinstance(values, (int, str)) else values))
    data_type = determine_gene_type(value_list)
    if not all(v == data_type[value_list[0]] for v in data_type.values()):
        raise ValueError("All items in ids must be of the same type (Entrez, Ensembl, or symbols).")

    chunk_size = 1000
    taxon_str = str(taxon)
    scope: T_MG_SCOPE = next(iter(data_type.values()))
    data = []
    chunks = range(0, len(value_list), chunk_size)

    for i in tqdm(chunks, desc=f"Getting info for '{scope}'"):
        result = info.get_queries(
            query=",".join(map(str, value_list[i : i + chunk_size])),
            dotfield=True,
            scopes=scope,
            fields="ensembl.gene,entrezgene,symbol,genomic_pos.start,genomic_pos.end,taxid,notfound",
            species=taxon_str,
        )
        if isinstance(result, int) and result == 414:
            raise ValueError(
                f"Query too long. Reduce the number of IDs in each query chunk (current chunk size: {chunk_size})."
            )
        if not isinstance(result, list):
            raise TypeError(f"Expected results to be a list, but got {type(result)}")
        if not isinstance(result[0], dict):
            raise TypeError(f"Expected each result to be a dict, but got {type(result[0])}")
        data.extend(result)
    return data


def get_remaining_identifiers(ids: T_IDS, taxon: int | str, cache: bool = True):
    """Convert between genomic identifiers.

    This function will convert between the following components:
        - Entrez Gene ID
        - Ensembl Gene ID
        - Gene Symbol

    :param ids: IDs to be converted
    :param taxon: Taxonomic identifier
    :param: scope: The type of identifier provided in `ids`
    :param cache: Should local caching be used for queries
    :return: DataFrame with columns "entrez_gene_id", "ensembl_gene_id", and "gene_symbol"
    """
    my_geneinfo = MyGeneInfo(cache=cache)
    gene_data = _get_conversion(info=my_geneinfo, values=ids, taxon=taxon)
    df = (
        pd.json_normalize(gene_data)
        .rename(
            columns={
                "ensembl.gene": "ensembl_gene_id",
                "entrezgene": "entrez_gene_id",
                "symbol": "gene_symbol",
                "taxid": "taxon_id",
            }
        )
        .drop(
            columns=["query", "_id", "_score", "genomic_pos.end", "genomic_pos.start"],
            errors="ignore",
        )
    )
    df = df[df["taxon_id"] == taxon]
    df["taxon_id"] = df["taxon_id"].astype(int, copy=True)

    # BUG: For an unknown reason, some Ensembl IDs are actually Entrez IDs
    #   To filter these, two approaches can be done:
    #   1) Remove rows where Ensembl IDs are integers
    #   2) Remove rows where Ensembl IDs equal Entrez IDs
    # We are selecting option 1 because it goes for the root cause: Ensembl IDs are not pure integers
    mask = df["ensembl_gene_id"].astype(str).str.fullmatch(r"\d+").fillna(False)
    df = df[
        (df["ensembl_gene_id"].astype("string").notna())  # remove NA values
        & (~df["ensembl_gene_id"].astype("string").str.fullmatch(r"\d+"))  # remove Entrez IDs
    ]
    return df


def _to_scalar(val) -> int:
    """Calculate the distance between end (e) and start (s)."""
    if isinstance(val, list):
        return int(sum(val) / len(val)) if val else 0  # `if val` checks that the list contains items
    if pd.isna(val):
        return 0
    return int(val)


def build_gene_info(ids: T_IDS, taxon: int | str, cache: bool = True):
    """Get genomic information from a given set of IDs.

    The input should be of the same type, otherwise this function will fail.
    Expected types are:
        - Ensembl Gene ID
        - Entrez Gene ID
        - Gene Symbol

    The returned data frame will have the following columns:
        - ensembl_gene_id
        - entrez_gene_id
        - gene_symbol
        - size (distance between genomic end and start)

    :param ids: IDs to be converted
    :param taxon: Taxonomic identifier
    :param cache: Should local caching be used for queries
    :return: pandas.DataFrame
    """
    my_geneinfo = MyGeneInfo(cache=cache)
    gene_data = _get_conversion(info=my_geneinfo, values=ids, taxon=taxon)
    df = pd.json_normalize(gene_data).rename(columns={"taxid": "taxon_id"})
    df = df[df["taxon_id"] == taxon]
    df["taxon_id"] = df["taxon_id"].astype(int, copy=True)

    df["size"] = df["genomic_pos.end"].fillna(0).map(_to_scalar) - df["genomic_pos.start"].fillna(0).map(_to_scalar)
    df = (
        df[~(df["size"] == 0)]
        .drop(
            columns=[
                "query",
                "_id",
                "_score",
                "genomic_pos.start",
                "genomic_pos.end",
                "notfound",
            ],
            inplace=False,
            errors="ignore",
        )
        .rename(
            columns={
                "ensembl.gene": "ensembl_gene_id",
                "entrezgene": "entrez_gene_id",
                "symbol": "gene_symbol",
            }
        )
        .explode(column=["ensembl_gene_id"])
        .sort_values(by="ensembl_gene_id", inplace=False)
    )

    return df


@overload
def determine_gene_type(items: str, /) -> T_MG_SCOPE: ...


@overload
def determine_gene_type(items: Sequence[str], /) -> dict[str, T_MG_SCOPE]: ...


def determine_gene_type(items: str | Sequence[str], /) -> str | dict[str, T_MG_SCOPE]:
    """Determine the genomic data type.

    :param items: A string or list of strings representing gene identifiers.
        The function will determine whether each identifier is an Entrez Gene ID,
        Ensembl Gene ID, or a gene symbol based on its format.

    :return: A dictionary mapping each input item to its determined type, which can be one of:
        - "entrez_gene_id": If the item consists solely of digits.
        - "ensembl_gene_id": If the item starts with "ENS" and is
            followed by a specific format (length greater than 11 and the last 11 characters are digits).
        - "gene_symbol": If the item does not match the above criteria, it is assumed to be a gene symbol.
    """
    values = (items,) if isinstance(items, str) else items
    result: dict[str, Literal["entrezgene", "ensembl.gene", "symbol"]] = {}

    for i in values:
        s = str(i).partition(".")[0]  # remove any transcripts that may exist

        if s.startswith("ENS") and len(s) > 11 and s[-11:].isdigit():
            result[s] = "ensembl.gene"
        elif s.isdigit():
            result[s] = "entrezgene"
        else:
            result[s] = "symbol"

    if isinstance(items, str):
        return result[items]
    return result


def contains_identical_gene_types(values: dict[str, T_MG_SCOPE] | Sequence[T_MG_SCOPE]) -> bool:
    data = values if not isinstance(values, dict) else list(values.values())
    return all(v == data[0] for v in data)


if __name__ == "__main__":
    r = build_gene_info(
        ids=[
            "ENSG00000284484",
            "ENSG00000299311",
            "ENSG00000202151",
            "ENSG00000226053",
            "ENSG00000131831",
        ],
        taxon=9606,
        cache=True,
    )
    print(r.columns)
