import os
import tempfile
import shutil
from nose.plugins import Plugin
import logging
from targetexplorer.utils import installation_testdir_filepath
from targetexplorer.initproject import InitProject


class SetUpTmpDbPlugin(Plugin):
    name = 'setup_tmp_db_plugin'
    enabled = True

    def configure(self, options, conf):
        if options.capture:
            logging.disable(logging.DEBUG)
            logging.disable(logging.INFO)

    def begin(self):
        self.temp_dir = tempfile.mkdtemp()
        with open(installation_testdir_filepath, 'w') as installation_testdir_file:
            installation_testdir_file.write(self.temp_dir)
        os.chdir(self.temp_dir)
        ip = InitProject(db_name='test', run_main=False)
        ip.setup()
        ip.mk_project_dirs()
        ip.mk_project_config_file()
        ip.write_wsgi_file()
        ip.create_db()

    def finalize(self, result):
        shutil.rmtree(self.temp_dir)
