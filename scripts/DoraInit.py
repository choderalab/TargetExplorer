#!/usr/bin/env python
import os, argparse, shutil, datetime
argparser = argparse.ArgumentParser(description='Initialize TargetExplorer database')
argparser.add_argument('--db_name', type=str, required=True, help='Database name, without extension')
args = argparser.parse_args()

print 'Initializing database project directory...'

tedb_basedir = os.path.abspath( os.path.join( os.path.dirname(__file__), '..' ) )

# make external-data dir
if not os.path.exists('external-data'):
    os.mkdir('external-data')

# make user-editable config file
if not os.path.exists('project_config.py'):
    with open('project_config.py', 'w') as new_config_file:
        new_config_file.write('''import os
# Edit this part
uniprot_query_string = 'EXAMPLE... domain:"protein kinase" AND reviewed:yes'
uniprot_domain_regex = 'EXAMPLE... ^Protein kinase(?!; truncated)(?!; inactive)'

# Don't edit this part
db_name = ''' + '\'%s\'' % (args.db_name) + '''
project_basedir = os.path.abspath(os.path.dirname(__file__))
targetexplorer_install_dir = ''' + '\'%s\'' % tedb_basedir)

# copy wsgi file
wsgi_filepath = args.db_name + '-wsgi.py'
if not os.path.exists(wsgi_filepath):
    wsgi_src_filepath = os.path.join(tedb_basedir, 'resources', 'template-wsgi.py')
    shutil.copy(wsgi_src_filepath, wsgi_filepath)

import flaskapp, flaskapp_config
import TargetExplorer

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
print 'Please now edit the file project_config.py before running the database generation scripts.'

