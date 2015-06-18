import os
from flask import Flask
from flask.ext.sqlalchemy import SQLAlchemy

targetexplorer_flaskapp_dir = os.path.dirname(__file__)

app = Flask(__name__)
db = SQLAlchemy(app)
app.config.from_pyfile(
    os.path.join(
        targetexplorer_flaskapp_dir,
        'config.py'
    )
)

from targetexplorer.flaskapp import views, models
