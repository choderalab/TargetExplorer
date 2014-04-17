#!/usr/bin/env python
import os.path, argparse
argparser = argparse.ArgumentParser(description='Initialize TargetExplorer database')
argparser.add_argument('--db_name', type=str, required=True)
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
uniprot_query_string = 'EXAMPLE... (taxonomy:9606 AND domain:"protein kinase") AND reviewed:yes'
uniprot_query_string_url = 'EXAMPLE... ?query=%28taxonomy%3A9606+AND+domain%3A%22protein+kinase%22%29+AND+reviewed%3Ayes&sort=score&format=xml'
uniprot_domain_regex = 'EXAMPLE... ^Protein kinase(?!; truncated)(?!; inactive)'

# Don't edit this part
DB_NAME = ''' + '\'%s\'' % (args.db_name) + '''
BASEDIR = os.path.abspath(os.path.dirname(__file__))
targetexplorer_install_dir = ''' + '\'%s\'' % tedb_basedir)

# copy wsgi file (should not need to be edited)
wsgi_filepath = args.db_name + '.wsgi'
if not os.path.exists(wsgi_filepath):
    wsgi_src_filepath = os.path.join(tedb_basedir, 'resources', 'db.wsgi')
    import shutil
    shutil.copy(wsgi_src_filepath, wsgi_filepath)

from app_config import SQLALCHEMY_DATABASE_URI
from app_config import SQLALCHEMY_MIGRATE_REPO
from app import db, models

# create database
db.create_all()

# add empty version data
version_row = models.Version()
db.session.add(version_row)
db.session.commit()

# TODO not certain yet whether using sqlalchemy-migrate; if not, will remove the following
# from migrate.versioning import api as migapi
# if not os.path.exists(SQLALCHEMY_MIGRATE_REPO):
#     migapi.create(SQLALCHEMY_MIGRATE_REPO, 'database repository')
#     migapi.version_control(SQLALCHEMY_DATABASE_URI, SQLALCHEMY_MIGRATE_REPO)
# else:
#     migapi.version_control(SQLALCHEMY_DATABASE_URI, SQLALCHEMY_MIGRATE_REPO, migapi.version(SQLALCHEMY_MIGRATE_REPO))

print 'Done.'
print 'Please now edit the file config.py before running the database generation scripts.'

