#!/usr/bin/env python
import re
import datetime
import argparse
import targetexplorer
from targetexplorer.flaskapp import db, models
import flask_sqlalchemy
import project_config

crawldata_row = models.CrawlData.query.first()
current_crawl_number = crawldata_row.current_crawl_number
safe_crawl_datestamp = crawldata_row.safe_crawl_datestamp
current_crawl_datestamps_row = models.DateStamps.query.filter_by(crawl_number=current_crawl_number).first()

data_problem = False

argparser = argparse.ArgumentParser(description='Commit database')
args = argparser.parse_args()

# ===================
# Test whether each of the scripts have been run since previous safe crawl.
# ===================
for data_type in ['uniprot', 'ncbi_gene', 'bindingdb', 'pdb', 'cbioportal']:
    datestamp_type = data_type + '_datestamp'
    current_crawl_datatype_datestamp = getattr(current_crawl_datestamps_row, datestamp_type)
    if current_crawl_datatype_datestamp == None:
        print 'data_type "%s" FAIL: no data found in db' % data_type
        data_problem = True
    elif current_crawl_datatype_datestamp <= safe_crawl_datestamp:
        print 'data_type "%s" FAIL: current data (%s) is older than or as old as safe-crawl data (%s)' % (data_type, current_crawl_datatype_datestamp.strftime(targetexplorer.core.datestamp_format_string), safe_crawl_datestamp.strftime(targetexplorer.core.datestamp_format_string))
        data_problem = True
    elif current_crawl_datatype_datestamp > safe_crawl_datestamp:
        print 'data_type "%s" PASS: current data (%s) is newer than safe-crawl data (%s)' % (data_type, current_crawl_datatype_datestamp.strftime(targetexplorer.core.datestamp_format_string), safe_crawl_datestamp.strftime(targetexplorer.core.datestamp_format_string))

if data_problem:
    raise Exception, 'Commit aborted.'
else:
    print 'Proceeding to commit to master db...'

# ===================
# Update crawl numbers
# ===================
crawldata_row.safe_crawl_number = current_crawl_number
crawldata_row.current_crawl_number = current_crawl_number + 1

# ===================
# Update datestamps data
# ===================
now = datetime.datetime.utcnow()
current_crawl_datestamps_row.commit_datestamp = now
new_datestamps_row = models.DateStamps(crawl_number=current_crawl_number+1)
db.session.add(new_datestamps_row)

# ===================
# Delete old crawls
# ===================
crawl_numbers = [row.crawl_number for row in models.DateStamps.query.all()]
if len(crawl_numbers) > project_config.ncrawls_to_save:
    print 'More than %d crawls found.' % project_config.ncrawls_to_save
    crawl_numbers_sorted = sorted(crawl_numbers, reverse=True)
    crawls_to_delete = crawl_numbers_sorted[project_config.ncrawls_to_save:]
    # iterate through crawls to delete
    for crawl_to_delete in crawls_to_delete:
        print 'Deleting crawl %d...' % crawl_to_delete
        table_classes = [value for value in models.__dict__.values() if isinstance(value, flask_sqlalchemy._BoundDeclarativeMeta)]
        for table in table_classes:
            if table == models.CrawlData:
                continue
            rows_to_delete = table.query.filter_by(crawl_number=crawl_to_delete)
            print '  - %s - %d rows' % (table.__name__, rows_to_delete.count())
            rows_to_delete.delete()

# ===================
# Write db to disk
# ===================
db.session.commit()

print 'Database committed.'
print 'New safe crawl number: %d' % (current_crawl_number)
print 'New current crawl number: %d' % (current_crawl_number+1)
print 'Done.'
