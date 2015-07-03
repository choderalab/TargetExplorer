from targetexplorer.flaskapp import models
from targetexplorer.tests.utils import projecttest_context
from targetexplorer.ncbi_gene import GatherNCBIGene


def test_gather_ncbi_gene():
    with projecttest_context(set_up_project_stage='uniprot'):
        GatherNCBIGene(use_existing_gene2pubmed=True)
        first_ncbi_gene_row = models.NCBIGeneEntry.query.first()
        assert isinstance(first_ncbi_gene_row, models.NCBIGeneEntry)
        assert first_ncbi_gene_row.gene_id == 25
        first_publication_row = models.NCBIGenePublication.query.first()
        assert isinstance(first_publication_row, models.NCBIGenePublication)
        assert first_publication_row.pmid == 1281542
