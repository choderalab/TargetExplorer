import os
basedir = os.path.abspath(os.path.dirname(__file__))

DB_NAME = 'kinome'

SQLALCHEMY_DATABASE_URI = 'sqlite:///' + os.path.join(basedir, DB_NAME + '.db')
SQLALCHEMY_MIGRATE_REPO = os.path.join(basedir, 'db_repository')

