from flask import render_template, abort
import config
import re
from app import app, db, models

@app.route('/')
@app.route('/index')
def index():
    return render_template('index.html', title=config.DB_NAME)

@app.route('/<string:ac>')
def get_dbentry(ac):
    # check ac is in proper UniProt format
    # TODO UniProt AC format will be extended some time after June 11 2014!
    if not (re.match('[A-N,R-Z][0-9][A-Z][A-Z,0-9][A-Z,0-9][0-9]', ac) or re.match('[O,P,Q][0-9][A-Z,0-9][A-Z,0-9][A-Z,0-9][0-9]', ac)):
        abort(404)

    uniprot = db.session.query(models.UniProt).filter_by(ac=ac).first()

    return render_template('dbentry.html', ac=ac, entry_name=uniprot.entry_name, family=uniprot.family)
