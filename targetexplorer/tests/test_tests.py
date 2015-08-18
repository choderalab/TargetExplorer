import os
from targetexplorer.tests.utils import projecttest_context
from nose.plugins.attrib import attr


@attr('unit')
def test_projecttest_context():
    with projecttest_context() as temp_dir:
        assert os.path.exists(temp_dir)


@attr('unit')
def test_testdir_file_writeable():
    from targetexplorer.utils import installation_testdir_filepath
    assert os.access(installation_testdir_filepath, os.W_OK)
