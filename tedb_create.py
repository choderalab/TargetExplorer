#!/usr/bin/env python
from migrate.versioning import api as migapi
from config import SQLALCHEMY_DATABASE_URI
from config import SQLALCHEMY_MIGRATE_REPO
from app import db
import os.path

db.create_all()

if not os.path.exists(SQLALCHEMY_MIGRATE_REPO):
    migapi.create(SQLALCHEMY_MIGRATE_REPO, 'database repository')
    migapi.version_control(SQLALCHEMY_DATABASE_URI, SQLALCHEMY_MIGRATE_REPO)
else:
    migapi.version_control(SQLALCHEMY_DATABASE_URI, SQLALCHEMY_MIGRATE_REPO, migapi.version(SQLALCHEMY_MIGRATE_REPO))
