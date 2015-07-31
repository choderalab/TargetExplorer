import os
import shutil
import datetime
import targetexplorer
from targetexplorer.utils import get_installed_resource_filepath
from targetexplorer.core import write_yaml_file, logger, project_config_filename, database_filename
from targetexplorer.core import external_data_dirpath, wsgi_filename, manual_overrides_filename


class InitProject(object):
    def __init__(self, db_name=None, project_path=None,
                 uniprot_query='EXAMPLE... mnemonic:ABL1_HUMAN',
                 uniprot_domain_regex='EXAMPLE... ^Protein kinase(?!; truncated)(?!; inactive)',
                 ncrawls_to_save=5,
                 run_main=True
                 ):
        self.db_name = db_name
        self.project_path = os.getcwd() if project_path is None else project_path
        self.uniprot_query = uniprot_query
        self.uniprot_domain_regex = uniprot_domain_regex
        self.ncrawls_to_save = ncrawls_to_save
        if run_main:
            self.setup()
            self.mk_project_dirs()
            self.mk_project_config_file()
            self.write_wsgi_file()
            self.write_manual_overrides_file()
            self.create_db()
            self.initialize_crawldata_and_datestamps()
            self.finish()

    def setup(self):
        self.targetexplorer_install_dir = os.path.abspath(
            os.path.dirname(targetexplorer.__file__)
        )
        self.sqlalchemy_database_uri = 'sqlite:///' + os.path.join(self.project_path, database_filename)

    def mk_project_dirs(self):
        if not os.path.exists(external_data_dirpath):
            os.mkdir(external_data_dirpath)

    def mk_project_config_file(self):
        if not os.path.exists(project_config_filename):
            config_data = {
                'db_name': self.db_name,
                'sqlalchemy_database_uri': self.sqlalchemy_database_uri,
                'dbapi_name': self.db_name + 'DBAPI',
                'ncrawls_to_save': self.ncrawls_to_save,
                'uniprot_query': self.uniprot_query,
                'uniprot_domain_regex': self.uniprot_domain_regex,
                'ignore_uniprot_pdbs': None,
            }
            write_yaml_file(config_data, project_config_filename)

    def write_wsgi_file(self):
        if not os.path.exists(wsgi_filename):
            template_wsgi_filepath = get_installed_resource_filepath(
                os.path.join('resources', 'template-wsgi.py')
            )
            shutil.copy(template_wsgi_filepath, wsgi_filename)

    def write_manual_overrides_file(self):
        if not os.path.exists(manual_overrides_filename):
            template_manual_overrides_filepath = get_installed_resource_filepath(
                os.path.join('resources', 'template-manual_overrides.yaml')
            )
            shutil.copy(template_manual_overrides_filepath, manual_overrides_filename)

    def create_db(self):
        from targetexplorer.flaskapp import db, app
        app.config.update(
            SQLALCHEMY_DATABASE_URI=self.sqlalchemy_database_uri
        )
        db.create_all()

    def initialize_crawldata_and_datestamps(self):
        """Add crawldata and empty datestamps row"""
        from targetexplorer.flaskapp import db, models
        current_crawl_number = 0
        safe_crawl_number = -1
        now = datetime.datetime.utcnow()
        crawldata_row = models.CrawlData(
            current_crawl_number=current_crawl_number,
            safe_crawl_number=safe_crawl_number,
            safe_crawl_datestamp=now,
        )
        db.session.add(crawldata_row)
        datestamps_row = models.DateStamps(crawl_number=current_crawl_number)
        db.session.add(datestamps_row)
        db.session.commit()

    def finish(self):
        logger.info('Done.')
