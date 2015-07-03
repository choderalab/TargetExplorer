from targetexplorer.flaskapp import models
from targetexplorer.tests.utils import projecttest_context
from targetexplorer.cbioportal_tmp import GatherCbioportalData


def test_gather_cbioportal():
    with projecttest_context(set_up_project_stage='uniprot'):
        GatherCbioportalData()
        first_mutation_row = models.CbioportalMutation.query.first()
        first_mutation_in_domain_row = models.CbioportalMutation.query.filter_by(
            in_uniprot_domain=True
        ).first()
        assert isinstance(first_mutation_row, models.CbioportalMutation)
        assert isinstance(first_mutation_in_domain_row.oncotator_aa_pos, int)
        assert first_mutation_in_domain_row.uniprot_domain.targetid == 'ABL1_HUMAN_D0'
