#!/usr/bin/env python
from migrate.versioning import api as migapi
from app_config import SQLALCHEMY_DATABASE_URI
from app_config import SQLALCHEMY_MIGRATE_REPO
from app import db, models
import os.path

db.create_all()

if not os.path.exists('external-data'):
    os.mkdir('external-data')

# add empty version data
version_row = models.Version()
db.session.add(version_row)
db.session.commit()

# if not os.path.exists(SQLALCHEMY_MIGRATE_REPO):
#     migapi.create(SQLALCHEMY_MIGRATE_REPO, 'database repository')
#     migapi.version_control(SQLALCHEMY_DATABASE_URI, SQLALCHEMY_MIGRATE_REPO)
# else:
#     migapi.version_control(SQLALCHEMY_DATABASE_URI, SQLALCHEMY_MIGRATE_REPO, migapi.version(SQLALCHEMY_MIGRATE_REPO))
