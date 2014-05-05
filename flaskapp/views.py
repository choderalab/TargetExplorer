from flask import abort, jsonify, request, make_response
import config
import re
from flaskapp import app, db, models

# ======
# error handlers
# ======

@app.errorhandler(AssertionError)
def assertion_error(error):
    print error.message
    return make_response(jsonify( { 'error': error.message } ), 404)

@app.errorhandler(404)
def not_found(error):
    return make_response(jsonify( { 'error': 'Not found' } ), 404)

# ======
# Get data for individual target
# ======

@app.route('/<string:ac>', methods = ['GET'])
def get_dbentry(ac):
    # check ac is in proper UniProt format
    # XXX TODO UniProt AC format will be extended some time after June 11 2014!
    try: assert (re.match('[A-N,R-Z][0-9][A-Z][A-Z,0-9][A-Z,0-9][0-9]', ac) or re.match('[O,P,Q][0-9][A-Z,0-9][A-Z,0-9][A-Z,0-9][0-9]', ac))
    except AssertionError as e: e.message = 'Incorrect UniProt AC format'; raise e

    uniprot = db.session.query(models.UniProt).filter_by(ac=ac).first()
    try: assert uniprot != None
    except AssertionError as e: e.message = 'Database entry not found'; raise e

    return jsonify( {'ac': uniprot.ac, 'entry_name': uniprot.entry_name, 'family': uniprot.family} )
