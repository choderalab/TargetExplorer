import os
from project_config import db_name, project_basedir

tedb_basedir = os.path.abspath(os.path.dirname(__file__))
sqlite_db_path = 'sqlite:///' + os.path.join(project_basedir, db_name) + '.db'
sqlite_db_stage_path = 'sqlite:///' + os.path.join(project_basedir, db_name) + '-stage.db'
dbapi_name = db_name + 'DBAPI'

SQLALCHEMY_DATABASE_URI = sqlite_db_path
