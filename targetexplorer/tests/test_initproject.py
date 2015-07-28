import os
from targetexplorer.tests.utils import projecttest_context
from targetexplorer.core import project_config_filename, database_filename, external_data_dirpath
from targetexplorer.core import wsgi_filename, manual_overrides_filename
from targetexplorer.initproject import InitProject
from nose.plugins.attrib import attr
from targetexplorer.flaskapp import app, models


@attr('unit')
def test_init_project():
    with projecttest_context(set_up_project_stage='blank'):
        InitProject(db_name='test', project_path=os.getcwd())
        assert os.path.exists(project_config_filename)
        assert os.path.exists(database_filename)
        assert os.path.exists(external_data_dirpath)
        assert os.path.exists(wsgi_filename)
        assert os.path.exists(manual_overrides_filename)
        assert app.config.get('SQLALCHEMY_DATABASE_URI') is not None
        crawldata_row = models.CrawlData.query.first()
        assert crawldata_row is not None
        assert crawldata_row.current_crawl_number == 0
        datestamps_row = models.DateStamps.query.first()
        assert datestamps_row is not None
        assert datestamps_row.crawl_number == 0
