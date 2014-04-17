import os
from config import BASEDIR
from config import DB_NAME

tedb_basedir = os.path.abspath(os.path.dirname(__file__))

SQLALCHEMY_DATABASE_URI = 'sqlite:///' + os.path.join(BASEDIR, DB_NAME + '.db')
SQLALCHEMY_MIGRATE_REPO = os.path.join(BASEDIR, 'db_repository')

