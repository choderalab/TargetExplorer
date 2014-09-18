def test_all_table_classes_have_crawl_number_attrib():
    from targetexplorer.flaskapp import models
    for table_class_name in models.table_class_names:
        if table_class_name == 'CrawlData':
            continue
        table_class = getattr(models, table_class_name)
        assert hasattr(table_class, 'crawl_number')
