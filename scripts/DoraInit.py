#!/usr/bin/env python
import os, argparse, shutil
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

# flaskapp.app.config['SQLALCHEMY_DATABASE_URI'] = flaskapp_config.sqlite_db_stage_path
TargetExplorer.core.select_stage_db()

# create database
flaskapp.db.create_all()

# add empty version data
version_row = flaskapp.models.Version(version_id=0, uniprot_datestamp=None, pdb_datestamp=None)
flaskapp.db.session.add(version_row)
flaskapp.db.session.commit()

# copy db_stage to db
shutil.copy(flaskapp_config.sqlite_db_stage_filename, flaskapp_config.sqlite_db_filename)

print 'Done.'
print 'Please now edit the file project_config.py before running the database generation scripts.'

