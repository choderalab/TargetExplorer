from targetexplorer.initproject import initialize_crawldata_and_datestamps
from targetexplorer.tests.utils import projecttest_context
from targetexplorer.flaskapp import models


def test_all_table_classes_have_crawl_number_attrib():
    from targetexplorer.flaskapp import models
    for table_class_name in models.table_class_names:
        if table_class_name == 'CrawlData':
            continue
        table_class = getattr(models, table_class_name)
        assert hasattr(table_class, 'crawl_number')


def test_crawl_data():
    crawl_data_row = models.CrawlData(
        current_crawl_number=0,
        safe_crawl_number=-1,
    )
    assert crawl_data_row.current_crawl_number == 0


def test_initialize_crawldata_and_datestamps():
    with projecttest_context(set_up_project_stage='blankdb'):
        initialize_crawldata_and_datestamps()
        crawl_data_from_db = models.CrawlData.query.first()
        assert crawl_data_from_db.current_crawl_number == 0
