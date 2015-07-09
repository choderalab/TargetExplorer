import datetime
import targetexplorer
from targetexplorer.flaskapp import db, models
from targetexplorer.core import read_project_config, logger


class Commit(object):
    def __init__(self, run_main=True):
        self.project_config = read_project_config()
        self.crawldata_row = models.CrawlData.query.first()
        self.current_crawl_number = self.crawldata_row.current_crawl_number
        self.safe_crawl_datestamp = self.crawldata_row.safe_crawl_datestamp
        self.current_crawl_datestamps_row = models.DateStamps.query.filter_by(crawl_number=self.current_crawl_number).first()
        if run_main:
            self.check_all_gather_scripts_have_been_run()
            self.update_crawl_numbers()
            self.update_datestamps()
            self.delete_old_crawls()
            self.commit()

    def check_all_gather_scripts_have_been_run(self):
        """
        Test whether each of the gather scripts have been run,
        and whether they have been updated in the correct order
        """
        data_problem = False
        for data_type in ['uniprot', 'ncbi_gene', 'bindingdb', 'pdb', 'cbioportal']:
            datestamp_type = data_type + '_datestamp'
            current_crawl_datatype_datestamp = getattr(self.current_crawl_datestamps_row, datestamp_type)
            if current_crawl_datatype_datestamp == None:
                logger.info('data_type "%s" FAIL: no data found in db' % data_type)
                data_problem = True
            elif current_crawl_datatype_datestamp <= self.safe_crawl_datestamp:
                logger.info('data_type "%s" FAIL: current data (%s) is older than or as old as safe-crawl data (%s)' % (data_type, current_crawl_datatype_datestamp.strftime(targetexplorer.core.datestamp_format_string), self.safe_crawl_datestamp.strftime(targetexplorer.core.datestamp_format_string)))
                data_problem = True
            elif current_crawl_datatype_datestamp > self.safe_crawl_datestamp:
                logger.info('data_type "%s" PASS: current data (%s) is newer than safe-crawl data (%s)' % (data_type, current_crawl_datatype_datestamp.strftime(targetexplorer.core.datestamp_format_string), self.safe_crawl_datestamp.strftime(targetexplorer.core.datestamp_format_string)))

        if data_problem:
            raise Exception('Commit aborted.')
        else:
            logger.info('Proceeding to commit to master db...')

    def update_crawl_numbers(self):
        self.crawldata_row.safe_crawl_number = self.current_crawl_number
        self.crawldata_row.current_crawl_number = self.current_crawl_number + 1

    def update_datestamps(self):
        now = datetime.datetime.utcnow()
        self.current_crawl_datestamps_row.commit_datestamp = now
        new_datestamps_row = models.DateStamps(crawl_number=self.current_crawl_number+1)
        db.session.add(new_datestamps_row)

    def delete_old_crawls(self):
        crawl_numbers = [row.crawl_number for row in models.DateStamps.query.all()]
        if len(crawl_numbers) > self.project_config['ncrawls_to_save']:
            logger.info('More than %d crawls found.' % self.project_config['ncrawls_to_save'])
            crawl_numbers_sorted = sorted(crawl_numbers, reverse=True)
            crawls_to_delete = crawl_numbers_sorted[self.project_config['ncrawls_to_save']:]
            # iterate through crawls to delete
            for crawl_to_delete in crawls_to_delete:
                logger.info('Deleting crawl %d...' % crawl_to_delete)
                # iterate through tables
                for table_class_name in models.table_class_names:
                    if table_class_name == 'CrawlData':
                        continue
                    table = getattr(models, table_class_name)
                    rows_to_delete = table.query.filter_by(crawl_number=crawl_to_delete)
                    logger.info('  - %s - %d rows' % (table_class_name, rows_to_delete.count()))
                    rows_to_delete.delete()

    def commit(self):
        db.session.commit()
        logger.info('Database committed.')
        logger.info('New safe crawl number: {0}'.format(self.current_crawl_number))
        logger.info('New current crawl number: {0}'.format(self.current_crawl_number+1))
        logger.info('Done.')
