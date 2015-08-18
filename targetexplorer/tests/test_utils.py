import os
from nose.plugins.attrib import attr


@attr('unit')
def test_installation_testdir_filepath():
    from targetexplorer.utils import installation_testdir_filepath
    assert os.path.exists(installation_testdir_filepath)
