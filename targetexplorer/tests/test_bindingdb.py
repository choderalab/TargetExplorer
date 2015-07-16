from targetexplorer.flaskapp import models
from targetexplorer.tests.utils import projecttest_context
from targetexplorer.bindingdb import GatherBindingDB
from nose.plugins.attrib import attr


@attr('unit')
def test_gather_bindingdb():
    with projecttest_context(set_up_project_stage='uniprot'):
        GatherBindingDB(use_existing_bindingdb_data=True)
        first_bioassay_row = models.BindingDBBioassay.query.first()
        assert isinstance(first_bioassay_row, models.BindingDBBioassay)
        assert first_bioassay_row.target_name == 'ABL1'


@attr('network')
@attr('slow')
def test_gather_bindingdb_using_network():
    with projecttest_context(set_up_project_stage='uniprot'):
        GatherBindingDB()
        first_bioassay_row = models.BindingDBBioassay.query.first()
        assert isinstance(first_bioassay_row, models.BindingDBBioassay)
        assert first_bioassay_row.target_name == 'ABL1'
