#!/usr/bin/env python
import os
import argparse
import shutil
import datetime
import targetexplorer

argparser = argparse.ArgumentParser(description='Initialize TargetExplorer database')
argparser.add_argument('--db_name', type=str, required=True, help='Database name, without extension')
args = argparser.parse_args()

print 'Initializing database project directory...'

targetexplorer_lib_dir = os.path.abspath(os.path.join(os.path.dirname(targetexplorer.__file__)))

# make external-data dir
if not os.path.exists('external-data'):
    os.mkdir('external-data')

# make user-editable config file
if not os.path.exists('project_config.py'):
    with open('project_config.py', 'w') as new_config_file:
        new_config_file.write('''import os
# UniProt search options
uniprot_query_string = 'EXAMPLE... domain:"protein kinase" AND reviewed:yes'
uniprot_domain_regex = 'EXAMPLE... ^Protein kinase(?!; truncated)(?!; inactive)'

# General database options
ncrawls_to_save = 5

# Don't edit the code below here
db_name = ''' + '\'%s\'' % (args.db_name) + '''
project_basedir = os.path.abspath(os.path.dirname(__file__))
targetexplorer_install_dir = ''' + '\'%s\'' % targetexplorer_lib_dir)

# copy wsgi file
wsgi_filepath = args.db_name + '-wsgi.py'
if not os.path.exists(wsgi_filepath):
    wsgi_src_filepath = os.path.join(targetexplorer_lib_dir, 'resources', 'template-wsgi.py')
    shutil.copy(wsgi_src_filepath, wsgi_filepath)

from targetexplorer import flaskapp

# create database
flaskapp.db.create_all()

# add crawldata and empty datestamps row
current_crawl_number = 0
safe_crawl_number = -1
now = datetime.datetime.utcnow()
crawldata_row = flaskapp.models.CrawlData(current_crawl_number=current_crawl_number, safe_crawl_number=safe_crawl_number, safe_crawl_datestamp=now)
flaskapp.db.session.add(crawldata_row)
datestamps_row = flaskapp.models.DateStamps(crawl_number=current_crawl_number)
flaskapp.db.session.add(datestamps_row)
flaskapp.db.session.commit()

print 'Done.'
print 'Please now edit the UniProt search options in project_config.py before running the database generation scripts.'
