# tests/test_gene_conversion.py
from __future__ import annotations

import importlib
import sys
import types
from dataclasses import dataclass

import pandas as pd
import pytest

MODULE_NAME = "como.pipelines.identifier"


@pytest.fixture
def mod(monkeypatch):
    """Import como in a way that does NOT require bioservices to be installed (we stub it before import).

    If bioservices *is* installed, this still keeps tests hermetic and avoids
    any accidental network calls.
    """
    # Stub bioservices.mygeneinfo.MyGeneInfo so importing the module works
    bioservices = types.ModuleType("bioservices")
    mygeneinfo = types.ModuleType("bioservices.mygeneinfo")

    class DummyMyGeneInfo:
        def __init__(self, cache: bool = True):
            self.cache = cache

        def get_queries(self, **kwargs):
            raise AssertionError(
                "DummyMyGeneInfo.get_queries was called unexpectedly. "
                "Patch mod.MyGeneInfo or use _get_conversion tests with a fake info object."
            )

    mygeneinfo.MyGeneInfo = DummyMyGeneInfo
    bioservices.mygeneinfo = mygeneinfo

    monkeypatch.setitem(sys.modules, "bioservices", bioservices)
    monkeypatch.setitem(sys.modules, "bioservices.mygeneinfo", mygeneinfo)

    # Ensure a clean import each test
    if MODULE_NAME in sys.modules:
        del sys.modules[MODULE_NAME]

    return importlib.import_module(MODULE_NAME)


def test_determine_gene_type_accepts_single_string(mod):
    assert mod.determine_gene_type("1017") == {"1017": "entrezgene"}


def test_determine_gene_type_classifies_common_formats_and_float(mod):
    # Includes a float to exercise the float-specific split(".")[0] branch.
    out = mod.determine_gene_type(["1017", "ENSG00000141510", "TP53", 123.0])

    assert out["1017"] == "entrezgene"
    assert out["ENSG00000141510"] == "ensembl.gene"
    assert out["TP53"] == "symbol"
    # float 123.0 -> "123"
    assert out["123"] == "entrezgene"


def test_determine_gene_type_ensembl_format_requires_ens_prefix_and_11_trailing_digits(mod):
    # Should be Ensembl-like: starts with ENS, length > 11, last 11 chars are digits
    assert mod.determine_gene_type(["ENSX00000123456"])["ENSX00000123456"] == "ensembl.gene"
    # Not enough trailing digits / too short => symbol
    assert mod.determine_gene_type(["ENS123"])["ENS123"] == "symbol"
    # ENS but last 11 not all digits => symbol
    assert mod.determine_gene_type(["ENSGABCDEFGHIJK"])["ENSGABCDEFGHIJK"] == "symbol"


# ---------------------------------------------------------------------------
# _get_conversion tests
# ---------------------------------------------------------------------------


@dataclass
class FakeInfo:
    return_value: object
    calls: list[dict] = None

    def __post_init__(self):  # noqa: D105
        self.calls = []

    def get_queries(self, **kwargs):
        self.calls.append(kwargs)
        return self.return_value


def test_get_conversion_raises_typeerror_when_results_not_list(mod):
    info = FakeInfo(return_value={"not": "a list"})
    with pytest.raises(TypeError, match="Expected results to be a list"):
        mod._get_conversion(info=info, values=["1017"], scope="entrezgene", fields="ensembl.gene,symbol", taxon="9606")


def test_get_conversion_raises_typeerror_when_first_item_not_dict(mod):
    info = FakeInfo(return_value=["not a dict"])
    with pytest.raises(TypeError, match="Expected each result to be a dict"):
        mod._get_conversion(info=info, values=["1017"], scope="entrezgene", fields="ensembl.gene,symbol", taxon="9606")


def test_get_conversion_calls_get_queries_with_expected_arguments(mod):
    info = FakeInfo(
        return_value=[
            {"query": "1017", "ensembl.gene": "ENSG00000123374", "symbol": "CDK2"},
        ]
    )

    out = mod._get_conversion(
        info=info, values=["1017"], scope="entrezgene", fields="ensembl.gene,symbol", taxon="9606"
    )

    assert out == [
        {"ensembl_gene_id": "ENSG00000123374", "entrez_gene_id": "1017", "gene_symbol": "CDK2"},
    ]

    assert len(info.calls) == 1
    call = info.calls[0]
    assert call["query"] == "1017"
    assert call["dotfield"] is True
    assert call["scopes"] == "entrezgene"
    assert call["fields"] == "ensembl.gene,symbol"
    assert call["species"] == "9606"


@pytest.mark.parametrize(
    ("scope", "results", "expected"),
    [
        (
            "entrezgene",
            [
                {"query": "1017", "ensembl.gene": "ENSG00000123374", "symbol": "CDK2"},
                {"query": "1018", "ensembl.gene": "ENSG00000123400", "symbol": "CDK3"},
            ],
            [
                {"ensembl_gene_id": "ENSG00000123374", "entrez_gene_id": "1017", "gene_symbol": "CDK2"},
                {"ensembl_gene_id": "ENSG00000123400", "entrez_gene_id": "1018", "gene_symbol": "CDK3"},
            ],
        ),
        (
            "ensembl.gene",
            [
                {"query": "ENSG00000141510", "entrezgene": "7157", "symbol": "TP53"},
            ],
            [
                {"ensembl_gene_id": "ENSG00000141510", "entrez_gene_id": "7157", "gene_symbol": "TP53"},
            ],
        ),
        (
            "symbol",
            [
                {"query": "TP53", "ensembl.gene": "ENSG00000141510", "entrezgene": "7157"},
            ],
            [
                {"ensembl_gene_id": "ENSG00000141510", "entrez_gene_id": "7157", "gene_symbol": "TP53"},
            ],
        ),
    ],
)
def test_get_conversion_scope_specific_mapping(mod, scope, results, expected):
    info = FakeInfo(return_value=results)
    out = mod._get_conversion(info=info, values=[r["query"] for r in results], scope=scope, fields="x", taxon="9606")
    assert out == expected


# ---------------------------------------------------------------------------
# convert tests
# ---------------------------------------------------------------------------


def test_convert_rejects_mixed_id_types(mod):
    # Mixed: Entrez-like ("1017") + symbol ("TP53") => ValueError
    with pytest.raises(ValueError, match="All items in ids must be of the same type"):
        mod.get_remaining_identifiers(ids=["1017", "TP53"], taxon=9606)


def test_convert_single_id_passes_expected_scope_fields_taxon_and_cache(monkeypatch, mod):
    calls = []

    class CapturingMyGeneInfo:
        def __init__(self, cache: bool = True):
            self.cache = cache

    def fake_get_conversion(*, info, values, scope, fields, taxon):
        calls.append({"info": info, "values": values, "scope": scope, "fields": fields, "taxon": taxon})
        # Return one row per value
        return [{"ensembl_gene_id": "ENSG_TEST", "entrez_gene_id": v, "gene_symbol": "SYM_TEST"} for v in values]

    monkeypatch.setattr(mod, "MyGeneInfo", CapturingMyGeneInfo)
    monkeypatch.setattr(mod, "_get_conversion", fake_get_conversion)

    df = mod.get_remaining_identifiers(ids=1017, taxon=9606, cache=False)

    assert isinstance(df, pd.DataFrame)
    assert set(df.columns) == {"ensembl_gene_id", "entrez_gene_id", "gene_symbol"}
    assert len(df) == 1
    assert df.loc[0, "entrez_gene_id"] == "1017"
    assert df.loc[0, "ensembl_gene_id"] == "ENSG_TEST"
    assert df.loc[0, "gene_symbol"] == "SYM_TEST"

    assert len(calls) == 1
    call = calls[0]
    assert isinstance(call["info"], CapturingMyGeneInfo)
    assert call["info"].cache is False
    assert call["values"] == ["1017"]
    assert call["scope"] == "entrezgene"
    assert call["taxon"] == "9606"

    # Fields are computed via set arithmetic; order is not guaranteed.
    assert set(call["fields"].split(",")) == {"ensembl.gene", "symbol"}


def test_convert_chunks_inputs_over_1000(monkeypatch, mod):
    calls = []

    class CapturingMyGeneInfo:
        def __init__(self, cache: bool = True):
            self.cache = cache

    def fake_get_conversion(*, info, values, scope, fields, taxon):
        calls.append({"values_len": len(values), "scope": scope, "fields": fields, "taxon": taxon})
        return [
            {"ensembl_gene_id": f"ENSG{v.zfill(11)}", "entrez_gene_id": v, "gene_symbol": f"SYM{v}"} for v in values
        ]

    monkeypatch.setattr(mod, "MyGeneInfo", CapturingMyGeneInfo)
    monkeypatch.setattr(mod, "_get_conversion", fake_get_conversion)

    ids = list(range(1, 1002))  # 1001 ids => should create 2 chunks: 1000 and 1
    df = mod.get_remaining_identifiers(ids=ids, taxon="9606")

    assert len(df) == 1001
    assert set(df.columns) == {"ensembl_gene_id", "entrez_gene_id", "gene_symbol"}

    # Ensure chunking happened as expected
    assert [c["values_len"] for c in calls] == [1000, 1]
    for c in calls:
        assert c["scope"] == "entrezgene"
        assert c["taxon"] == "9606"
        assert set(c["fields"].split(",")) == {"ensembl.gene", "symbol"}

    # Spot-check first and last rows keep ordering
    assert df.loc[0, "entrez_gene_id"] == "1"
    assert df.loc[0, "gene_symbol"] == "SYM1"
    assert df.loc[1000, "entrez_gene_id"] == "1001"
    assert df.loc[1000, "gene_symbol"] == "SYM1001"


def test_convert_symbol_scope_fields_exclude_symbol(monkeypatch, mod):
    calls = []

    class CapturingMyGeneInfo:
        def __init__(self, cache: bool = True):
            self.cache = cache

    def fake_get_conversion(*, info, values, scope, fields, taxon):
        calls.append({"scope": scope, "fields": fields, "taxon": taxon, "values": values})
        # Return minimal row(s)
        return [{"ensembl_gene_id": None, "entrez_gene_id": None, "gene_symbol": v} for v in values]

    monkeypatch.setattr(mod, "MyGeneInfo", CapturingMyGeneInfo)
    monkeypatch.setattr(mod, "_get_conversion", fake_get_conversion)

    df = mod.get_remaining_identifiers(ids=["TP53", "BRCA1"], taxon=9606)

    assert len(df) == 2
    assert set(df.columns) == {"ensembl_gene_id", "entrez_gene_id", "gene_symbol"}

    assert len(calls) == 1
    assert calls[0]["scope"] == "symbol"
    assert calls[0]["taxon"] == "9606"
    # When scope is "symbol", requested fields should be the other two
    assert set(calls[0]["fields"].split(",")) == {"ensembl.gene", "entrezgene"}
