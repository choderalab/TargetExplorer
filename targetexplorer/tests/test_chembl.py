__author__ = 'isikm'
from targetexplorer.flaskapp import models
from targetexplorer.tests.utils import projecttest_context
from targetexplorer.chembl import GatherChemblTarget
from nose.plugins.attrib import attr

@attr('unit')
def test_gather_chembl_target():
    with projecttest_context(set_up_project_stage='uniprot'):
        GatherChemblTarget()

        # first_mutation_row = models.CbioportalMutation.query.first()
        # assert first_mutation_row is not None
        # first_mutation_in_domain_row = models.CbioportalMutation.query.filter_by(
        #     in_uniprot_domain=True
        # ).first()
        # assert first_mutation_in_domain_row is not None
        # assert isinstance(first_mutation_in_domain_row.oncotator_aa_pos, int)
        # assert first_mutation_in_domain_row.uniprot_domain.description == 'SH2'

