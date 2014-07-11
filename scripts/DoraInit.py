#!/usr/bin/env python
import os, argparse
argparser = argparse.ArgumentParser(description='Initialize TargetExplorer database')
argparser.add_argument('--db_name', type=str, required=True, help='Database name, without extension')
args = argparser.parse_args()

print 'Initializing database project directory...'

tedb_basedir = os.path.abspath( os.path.join( os.path.dirname(__file__), '..' ) )

# make external-data dir
if not os.path.exists('external-data'):
    os.mkdir('external-data')

# make user-editable config file
if not os.path.exists('config.py'):
    with open('config.py', 'w') as new_config_file:
        new_config_file.write('''import os
# Edit this part
uniprot_query_string = 'EXAMPLE... domain:"protein kinase" AND reviewed:yes'
uniprot_domain_regex = 'EXAMPLE... ^Protein kinase(?!; truncated)(?!; inactive)'

# Don't edit this part
DB_NAME = ''' + '\'%s\'' % (args.db_name) + '''
dbapi_name = DB_NAME + 'DBAPI'
BASEDIR = os.path.abspath(os.path.dirname(__file__))
targetexplorer_install_dir = ''' + '\'%s\'' % tedb_basedir)

# copy wsgi file (should not need to be edited)
wsgi_filepath = args.db_name + '-wsgi.py'
if not os.path.exists(wsgi_filepath):
    wsgi_src_filepath = os.path.join(tedb_basedir, 'resources', 'template-wsgi.py')
    import shutil
    shutil.copy(wsgi_src_filepath, wsgi_filepath)

# from app_config import SQLALCHEMY_DATABASE_URI
# from app_config import SQLALCHEMY_MIGRATE_REPO
import app_master, app_stage

# create database
app_master.db.create_all()
app_stage.db.create_all()

# add empty version data
version_row = app_master.models.Version(version_id=0, uniprot_datestamp=None, pdb_datestamp=None)
app_master.db.session.add(version_row)
app_master.db.session.commit()

version_row = app_stage.models.Version(version_id=0, uniprot_datestamp=None, pdb_datestamp=None)
app_stage.db.session.add(version_row)
app_stage.db.session.commit()

print 'Done.'
print 'Please now edit the file config.py before running the database generation scripts.'

