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
    assert mod.determine_gene_type("1017") == "entrezgene"


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
        mod._get_conversion(info=info, values=["1017"], taxon="9606")


def test_get_conversion_raises_typeerror_when_first_item_not_dict(mod):
    info = FakeInfo(return_value=["not a dict"])
    with pytest.raises(TypeError, match="Expected each result to be a dict"):
        mod._get_conversion(info=info, values=["1017"], taxon="9606")


def test_get_conversion_calls_get_queries_with_expected_arguments(mod):
    info = FakeInfo(
        return_value=[
            {"query": "1017", "ensembl.gene": "ENSG00000123374", "symbol": "CDK2"},
        ]
    )
    
    out = mod._get_conversion(info=info, values=["1017"], taxon="9606")
    
    # _get_conversion returns raw API results
    assert out == [
        {"query": "1017", "ensembl.gene": "ENSG00000123374", "symbol": "CDK2"},
    ]
    
    assert len(info.calls) == 1
    call = info.calls[0]
    assert call["query"] == "1017"
    assert call["dotfield"] is True
    assert call["scopes"] == "entrezgene"
    assert set(call["fields"].split(",")) == {
        "ensembl.gene",
        "entrezgene",
        "symbol",
        "genomic_pos.start",
        "genomic_pos.end",
        "taxid",
        "notfound",
    }
    assert call["species"] == "9606"


@pytest.mark.parametrize(
    ("results",),
    [
        # Test data format: list of raw API results to pass through
        (
                [
                    {"query": "1017", "ensembl.gene": "ENSG00000123374", "symbol": "CDK2"},
                    {"query": "1018", "ensembl.gene": "ENSG00000123400", "symbol": "CDK3"},
                ],
        ),
        (
                [
                    {"query": "ENSG00000141510", "entrezgene": "7157", "symbol": "TP53"},
                ],
        ),
        (
                [
                    {"query": "TP53", "ensembl.gene": "ENSG00000141510", "entrezgene": "7157"},
                ],
        ),
    ],
)
def test_get_conversion_returns_raw_results(mod, results):
    """Test that _get_conversion returns raw API results without transformation."""
    info = FakeInfo(return_value=results)
    out = mod._get_conversion(info=info, values=[r["query"] for r in results], taxon="9606")
    assert out == results


def test_convert_rejects_mixed_id_types(mod):
    # Mixed: Entrez-like ("1017") + symbol ("TP53") => ValueError
    with pytest.raises(ValueError, match="All items in ids must be of the same type"):
        mod.get_remaining_identifiers(ids=["1017", "TP53"], taxon=9606)


def test_convert_single_id_passes_expected_scope_fields_taxon_and_cache(monkeypatch, mod):
    calls = []
    
    class CapturingMyGeneInfo:
        def __init__(self, cache: bool = True):
            self.cache = cache
    
    def fake_get_conversion(*, info, values, taxon):
        calls.append({"info": info, "values": values, "taxon": taxon})
        # Return raw API format with taxid for filtering
        value_list = [values] if isinstance(values, (int, str)) else values
        return [
            {"query": str(v), "ensembl.gene": "ENSG_TEST", "entrezgene": str(v), "symbol": f"SYM{v}", "taxid": taxon}
            for v in value_list
        ]
    
    monkeypatch.setattr(mod, "MyGeneInfo", CapturingMyGeneInfo)
    monkeypatch.setattr(mod, "_get_conversion", fake_get_conversion)
    
    df = mod.get_remaining_identifiers(ids=1017, taxon=9606, cache=False)
    
    assert isinstance(df, pd.DataFrame)
    assert set(df.columns) == {"ensembl_gene_id", "entrez_gene_id", "gene_symbol", "taxon_id"}
    assert len(df) == 1
    assert df.loc[0, "entrez_gene_id"] == "1017"
    assert df.loc[0, "ensembl_gene_id"] == "ENSG_TEST"
    assert df.loc[0, "gene_symbol"] == "SYM1017"
    assert df.loc[0, "taxon_id"] == 9606
    
    assert len(calls) == 1
    call = calls[0]
    assert isinstance(call["info"], CapturingMyGeneInfo)
    assert call["info"].cache is False
    assert call["values"] == 1017
    assert call["taxon"] == 9606


def test_convert_large_input_list(monkeypatch, mod):
    """Test that get_remaining_identifiers handles large input lists correctly."""
    calls = []
    
    class CapturingMyGeneInfo:
        def __init__(self, cache: bool = True):
            self.cache = cache
    
    def fake_get_conversion(*, info, values, taxon):
        calls.append({"values_len": len(values), "taxon": taxon})
        value_list = [values] if isinstance(values, (int, str)) else values
        return [
            {
                "query": str(v),
                "ensembl.gene": f"ENSG{str(v).zfill(11)}",
                "entrezgene": str(v),
                "symbol": f"SYM{v}",
                "taxid": taxon,
            }
            for v in value_list
        ]
    
    monkeypatch.setattr(mod, "MyGeneInfo", CapturingMyGeneInfo)
    monkeypatch.setattr(mod, "_get_conversion", fake_get_conversion)
    
    ids = list(range(1, 1002))
    df = mod.get_remaining_identifiers(ids=ids, taxon="9606")
    
    assert len(df) == 1001
    assert set(df.columns) == {"ensembl_gene_id", "entrez_gene_id", "gene_symbol", "taxon_id"}
    assert len(calls) == 1
    assert calls[0]["values_len"] == 1001
    assert calls[0]["taxon"] == "9606"
    
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
    
    def fake_get_conversion(*, info, values, taxon):
        calls.append({"taxon": taxon, "values": values})
        return [
            {
                "query": v,
                "ensembl.gene": f"ENSG{v}_TEST",
                "entrezgene": "123",
                "symbol": v,
                "taxid": taxon,
            }
            for v in values
        ]
    
    monkeypatch.setattr(mod, "MyGeneInfo", CapturingMyGeneInfo)
    monkeypatch.setattr(mod, "_get_conversion", fake_get_conversion)
    
    df = mod.get_remaining_identifiers(ids=["TP53", "BRCA1"], taxon=9606)
    
    assert len(df) == 2
    assert set(df.columns) == {"ensembl_gene_id", "entrez_gene_id", "gene_symbol", "taxon_id"}
    assert df.loc[0, "gene_symbol"] == "TP53"
    assert df.loc[1, "gene_symbol"] == "BRCA1"
    assert df.loc[0, "taxon_id"] == 9606
    
    assert len(calls) == 1
    assert calls[0]["taxon"] == 9606
