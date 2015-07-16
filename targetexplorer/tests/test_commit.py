from targetexplorer.flaskapp import models
from targetexplorer.tests.utils import projecttest_context, expected_failure
from targetexplorer.commit import Commit
from nose.plugins.attrib import attr
from targetexplorer.utils import DatabaseException


@attr('unit')
def test_commit():
    with projecttest_context(set_up_project_stage='cbioportal'):
        Commit()
        crawl_data_row = models.CrawlData.query.first()
        assert crawl_data_row.safe_crawl_number == 0


@attr('unit')
def test_premature_commit():
    with projecttest_context(set_up_project_stage='uniprot'):
        try:
            Commit()
        except DatabaseException:
            pass
