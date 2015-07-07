import os
import shutil
import datetime
import targetexplorer
from targetexplorer.utils import get_installed_resource_filepath
from targetexplorer.core import write_yaml_file, logger


class InitProject(object):
    def __init__(
            self, db_name=None, project_path=None,
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
            self.create_db()
            initialize_crawldata_and_datestamps()
            self.finish()

    def setup(self):
        self.targetexplorer_install_dir = os.path.abspath(os.path.join(os.path.dirname(targetexplorer.__file__)))
        self.project_config_filename = 'project_config.yaml'
        self.sqlalchemy_uri = 'sqlite:///' + os.path.join(self.project_path, 'database.db')
        self.wsgi_filepath = 'webapi-wsgi.py'

    def mk_project_dirs(self):
        if not os.path.exists('external-data'):
            os.mkdir('external-data')

    def mk_project_config_file(self):
        config_data = {
            'db_name': self.db_name,
            'sqlalchemy_uri': self.sqlalchemy_uri,
            'dbapi_name': self.db_name + 'DBAPI',
            'ncrawls_to_save': self.ncrawls_to_save,
            'uniprot_query': self.uniprot_query,
            'uniprot_domain_regex': self.uniprot_domain_regex,
            'ignore_uniprot_pdbs': None,
        }
        write_yaml_file(config_data, self.project_config_filename)

    def write_wsgi_file(self):
        if not os.path.exists(self.wsgi_filepath):
            template_wsgi_filepath = get_installed_resource_filepath(
                os.path.join('resources', 'template-wsgi.py')
            )

            shutil.copy(template_wsgi_filepath, self.wsgi_filepath)

    def create_db(self):
        from targetexplorer.flaskapp import db
        db.create_all()

    def finish(self):
        logger.info('Done.')


def initialize_crawldata_and_datestamps():
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
