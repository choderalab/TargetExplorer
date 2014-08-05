#!/usr/bin/env python
import TargetExplorer
import re
from flaskapp import db, models

metadata_row = models.MetaData.query.first()
current_crawl_number = metadata_row.current_crawl_number
safe_crawl_datestamp = metadata_row.safe_crawl_datestamp

data_problem = False

# Test whether each of the scripts have been run, and whether they have been updated in the correct order
for data_type in ['uniprot']:
    datestamp_type = 'current_' + data_type + '_datestamp'
    current_crawl_datatype_datestamp = getattr(metadata_row, datestamp_type)
    if current_crawl_datatype_datestamp == None:
        print 'data_type "%s" FAIL: no data in stage db' % data_type
        data_problem = True
    elif current_crawl_datatype_datestamp <= safe_crawl_datestamp:
        print 'data_type "%s" FAIL: staged data (%s) is older than or as old as master (%s)' % (data_type, current_crawl_datatype_datestamp.strftime(TargetExplorer.core.datestamp_format_string), safe_crawl_datestamp.strftime(TargetExplorer.core.datestamp_format_string))
        data_problem = True
    elif current_crawl_datatype_datestamp > safe_crawl_datestamp:
        print 'data_type "%s" PASS: staged data (%s) is newer than master (%s)' % (data_type, current_crawl_datatype_datestamp.strftime(TargetExplorer.core.datestamp_format_string), safe_crawl_datestamp.strftime(TargetExplorer.core.datestamp_format_string))

if data_problem:
    raise Exception, 'Commit aborted.'
else:
    print 'Proceeding to commit to master db...'

# update crawl numbers
metadata_row.safe_crawl_number = current_crawl_number
metadata_row.current_crawl_number = current_crawl_number + 1

# write db to disk
db.session.commit()

print 'Updated safe crawl number to %d' % (current_crawl_number)
print 'Done.'
