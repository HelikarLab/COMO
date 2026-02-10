from collections.abc import Sequence
from typing import Literal

import pandas as pd
from bioservices.mygeneinfo import MyGeneInfo

__all__ = [
    "convert",
    "determine_gene_type",
]

T_MG_SCOPE = Literal["entrezgene", "ensembl.gene", "symbol"]
T_MG_TRANSLATE = Literal["entrez_gene_id", "ensembl_gene_id", "gene_symbol"]
T_MG_RETURN = list[dict[T_MG_TRANSLATE, str]]


def _get_conversion(info: MyGeneInfo, values: list[str], scope: T_MG_SCOPE, fields: str, taxon: str) -> T_MG_RETURN:
    value_str = ",".join(map(str, values))
    results = info.get_queries(query=value_str, dotfield=True, scopes=scope, fields=fields, species=taxon)
    if not isinstance(results, list):
        raise TypeError(f"Expected results to be a list, but got {type(results)}")
    if not isinstance(results[0], dict):
        raise TypeError(f"Expected each result to be a dict, but got {type(results[0])}")

    data: T_MG_RETURN = []
    for result in results:
        ensembl = result.get("query" if scope == "ensembl.gene" else "ensembl.gene")
        entrez = result.get("query" if scope == "entrezgene" else "entrezgene")
        symbol = result.get("query" if scope == "symbol" else "symbol")
        data.append({"ensembl_gene_id": ensembl, "entrez_gene_id": entrez, "gene_symbol": symbol})
    return data


def convert(ids: int | str | Sequence[int] | Sequence[str] | Sequence[int | str], taxon: int | str, cache: bool = True):
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
    chunk_size = 1000
    id_list = list(map(str, [ids] if isinstance(ids, (int, str)) else ids))
    chunks = list(range(0, len(id_list), chunk_size))

    data_type = determine_gene_type(id_list)
    if not all(v == data_type[id_list[0]] for v in data_type.values()):
        raise ValueError("All items in ids must be of the same type (Entrez, Ensembl, or symbols).")

    scope = next(iter(data_type.values()))
    fields = ",".join({"ensembl.gene", "entrezgene", "symbol"} - {scope})
    taxon_str = str(taxon)
    return pd.DataFrame(
        [
            row
            for i in chunks
            for row in _get_conversion(
                info=my_geneinfo,
                values=id_list[i : i + chunk_size],
                scope=scope,
                fields=fields,
                taxon=taxon_str,
            )
        ]
    )


def determine_gene_type(items: str | list[str], /) -> dict[str, T_MG_SCOPE]:
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
    items = [items] if isinstance(items, str) else items

    determine: dict[str, Literal["entrezgene", "ensembl.gene", "symbol"]] = {}
    for i in items:
        i_str = str(i).split(".")[0] if isinstance(i, float) else str(i)
        if i_str.isdigit():
            determine[i_str] = "entrezgene"
        elif i_str.startswith("ENS") and (len(i_str) > 11 and all(i.isdigit() for i in i_str[-11:])):
            determine[i_str] = "ensembl.gene"
        else:
            determine[i_str] = "symbol"

    return determine
