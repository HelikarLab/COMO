import os
import sys

# Add parent directory to path, allows us to import the "project.py" file from the parent directory
# From: https://stackoverflow.com/a/30536516/13885200
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from async_bioservices import input_database, output_database, taxon_ids


class TestDatabases:
    """
    We've chosen a few random input databases to test.
    These values should never change.
    """
    def test_input_database(self):
        assert input_database.InputDatabase.GENE_ID.value == "Gene ID"
        assert input_database.InputDatabase.GENE_SYMBOL.value == "Gene Symbol"
        assert input_database.InputDatabase.TAXON_ID.value == "Taxon ID"
        assert input_database.InputDatabase.UNIGENE_ID.value == "UniGene ID"
        assert input_database.InputDatabase.MIM_ID.value == "MIM ID"
        assert input_database.InputDatabase.HGNC_ID.value == "HGNC ID"
        assert input_database.InputDatabase.HOMOLOGENE_ID.value == "HomoloGene ID"

    def test_output_database(self):
        assert output_database.OutputDatabase.CPDB_PROTEIN_INTERACTOR.value == "CPDB Protein Interactor"
        assert output_database.OutputDatabase.EC_NUMBER.value == "EC Number"
        assert output_database.OutputDatabase.GAD_DISEASE_INFO.value == "GAD Disease Info"
        assert output_database.OutputDatabase.GO_ID.value == "GO ID"
        assert output_database.OutputDatabase.HOMOLOG_HUMAN_ENS_GENE_ID.value == "Homolog - Human Ens Gene ID"
        assert output_database.OutputDatabase.KEGG_PATHWAY_ID.value == "KEGG Pathway ID"
        assert output_database.OutputDatabase.UNIPROT_ACCESSION.value == "UniProt Accession"

    def test_taxon_ids(self):
        assert taxon_ids.TaxonIDs.HUMAN.value == 9606


