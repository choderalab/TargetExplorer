from flask import Flask
from flask.ext.sqlalchemy import SQLAlchemy

app = Flask(__name__)
app.config.from_object('app_config')
db = SQLAlchemy(app)

from app_master import views, models
