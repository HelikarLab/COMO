import os
import sys

# Add parent directory to path, allows us to import the "project.py" file from the parent directory
# From: https://stackoverflow.com/a/30536516/13885200
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from async_bioservices.input_database import InputDatabase
from async_bioservices.output_database import OutputDatabase
from async_bioservices.taxon_ids import TaxonIDs


def test_input_database():
    """
    We've chosen a few random input databases to test.
    These values should never change.
    """
    assert InputDatabase.GENE_ID.value == "Gene ID"
    assert InputDatabase.GENE_SYMBOL.value == "Gene Symbol"
    assert InputDatabase.TAXON_ID.value == "Taxon ID"
    assert InputDatabase.UNIGENE_ID.value == "UniGene ID"
    assert InputDatabase.MIM_ID.value == "MIM ID"
    assert InputDatabase.HGNC_ID.value == "HGNC ID"
    assert InputDatabase.HOMOLOGENE_ID.value == "HomoloGene ID"


def test_output_database():
    """
    We've chosen a few random input databases to test.
    These values should never change.
    """
    assert OutputDatabase.CPDB_PROTEIN_INTERACTOR.value == "CPDB Protein Interactor"
    assert OutputDatabase.EC_NUMBER.value == "EC Number"
    assert OutputDatabase.GAD_DISEASE_INFO.value == "GAD Disease Info"
    assert OutputDatabase.GO_ID.value == "GO ID"
    assert OutputDatabase.HOMOLOG_HUMAN_ENS_GENE_ID.value == "Homolog - Human Ens Gene ID"
    assert OutputDatabase.KEGG_PATHWAY_ID.value == "KEGG Pathway ID"
    assert OutputDatabase.UNIPROT_ACCESSION.value == "UniProt Accession"


def test_taxon_ids():
    assert TaxonIDs.HUMAN.value == 9606
