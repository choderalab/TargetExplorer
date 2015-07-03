import os
from targetexplorer.flaskapp import models
from targetexplorer.tests.utils import projecttest_context
from targetexplorer.bindingdb_tmp import GatherBindingDB


def test_gather_pdb():
    with projecttest_context(set_up_project_stage='uniprot'):
        import platform
        if platform.system() == 'Darwin':
            grep_path = '/usr/local/bin/grep'
            if not os.path.exists(grep_path):
                raise Exception('Please use Homebrew version of grep')
        else:
            grep_path = False

        GatherBindingDB(use_existing_bindingdb_data=True, grep_path=grep_path)
        first_bioassay_row = models.BindingDBBioassay.query.first()
        assert isinstance(first_bioassay_row, models.BindingDBBioassay)
        assert first_bioassay_row.target_name == 'ABL1'
