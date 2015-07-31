from targetexplorer.flaskapp import models
from targetexplorer.tests.utils import projecttest_context
from targetexplorer.cbioportal import GatherCbioportalData, retrieve_extended_mutation_datatxt
from targetexplorer.oncotator import retrieve_oncotator_mutation_data_as_json
from nose.plugins.attrib import attr
from nose.plugins.skip import SkipTest


@attr('network')
def test_retrieve_oncotator_mutation_data_as_json():
    oncotator_result = retrieve_oncotator_mutation_data_as_json(
        7,
        55259515,
        55259515,
        'T',
        'G'
    )
    assert oncotator_result.get('transcript_id') == 'ENST00000275493.2'


@attr('network')
def test_retrieve_extended_mutation_datatxt():
    lines = retrieve_extended_mutation_datatxt(
        'gbm_tcga_all',
        'gbm_tcga_mutations',
        ['EGFR', 'PTEN']
    )
    assert len(lines) > 0


@attr('private_cbioportal')
def test_retrieve_extended_mutation_datatxt_private_portal():
    """
    So far seems to be impossible to access private data via web API.
    """
    raise SkipTest
    lines = retrieve_extended_mutation_datatxt(
        'gbm_tcga_all',
        'gbm_tcga_mutations',
        ['EGFR', 'PTEN'],
        portal_version='private'
    )
    print '\n'.join(lines)
    assert len(lines) > 0


@attr('unit')
def test_gather_cbioportal():
    with projecttest_context(set_up_project_stage='uniprot'):
        GatherCbioportalData(use_existing_cbioportal_data=True, use_existing_oncotator_data=True)
        first_mutation_row = models.CbioportalMutation.query.first()
        assert first_mutation_row is not None
        first_mutation_in_domain_row = models.CbioportalMutation.query.filter_by(
            in_uniprot_domain=True
        ).first()
        assert first_mutation_in_domain_row is not None
        assert isinstance(first_mutation_in_domain_row.oncotator_aa_pos, int)
        assert first_mutation_in_domain_row.uniprot_domain.description == 'SH2'


@attr('network')
def test_gather_cbioportal_using_network():
    with projecttest_context(set_up_project_stage='uniprot'):
        GatherCbioportalData()
        first_mutation_row = models.CbioportalMutation.query.first()
        first_mutation_in_domain_row = models.CbioportalMutation.query.filter_by(
            in_uniprot_domain=True
        ).first()
        assert isinstance(first_mutation_row, models.CbioportalMutation)
        assert isinstance(first_mutation_in_domain_row.oncotator_aa_pos, int)
        assert first_mutation_in_domain_row.uniprot_domain.description == 'SH2'
