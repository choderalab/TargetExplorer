import os
from project_config import db_name, project_basedir

sqlite_db_filename = db_name + '.db'
sqlite_db_path = 'sqlite:///' + os.path.join(project_basedir, db_name) + '.db'
dbapi_name = db_name + 'DBAPI'

SQLALCHEMY_DATABASE_URI = sqlite_db_path