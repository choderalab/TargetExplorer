from flask import abort, jsonify, request, make_response, current_app
import config
import re
from datetime import timedelta
from functools import update_wrapper
from flaskapp import app, db, models

# ======
# HTTP access control decorator (CORS)
# ======

def crossdomain(origin=None, methods=None, headers=None,
                max_age=21600, attach_to_all=True,
                automatic_options=True):
    if methods is not None:
        methods = ', '.join(sorted(x.upper() for x in methods))
    if headers is not None and not isinstance(headers, basestring):
        headers = ', '.join(x.upper() for x in headers)
    if not isinstance(origin, basestring):
        origin = ', '.join(origin)
    if isinstance(max_age, timedelta):
        max_age = max_age.total_seconds()

    def get_methods():
        if methods is not None:
            return methods

        options_resp = current_app.make_default_options_response()
        return options_resp.headers['allow']

    def decorator(f):
        def wrapped_function(*args, **kwargs):
            if automatic_options and request.method == 'OPTIONS':
                resp = current_app.make_default_options_response()
            else:
                resp = make_response(f(*args, **kwargs))
            if not attach_to_all and request.method != 'OPTIONS':
                return resp

            h = resp.headers

            h['Access-Control-Allow-Origin'] = origin
            h['Access-Control-Allow-Methods'] = get_methods()
            h['Access-Control-Max-Age'] = str(max_age)
            if headers is not None:
                h['Access-Control-Allow-Headers'] = headers
            return resp

        f.provide_automatic_options = False
        f.required_methods = ['OPTIONS']
        return update_wrapper(wrapped_function, f)
    return decorator

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

@app.route('/<path:leadingpath>/<string:ac>', methods = ['GET'])
@crossdomain(origin='*', headers=["Origin", "X-Requested-With", "Content-Type", "Accept"])
def get_dbentry(leadingpath, ac):
    # note: leadingpath is ignored
    # check ac is in proper UniProt format
    # XXX TODO UniProt AC format will be extended some time after June 11 2014!
    try: assert (re.match('[A-N,R-Z][0-9][A-Z][A-Z,0-9][A-Z,0-9][0-9]', ac) or re.match('[O,P,Q][0-9][A-Z,0-9][A-Z,0-9][A-Z,0-9][0-9]', ac))
    except AssertionError as e: e.message = 'Incorrect UniProt AC format'; raise e

    uniprot = db.session.query(models.UniProt).filter_by(ac=ac).first()
    try: assert uniprot != None
    except AssertionError as e: e.message = 'Database entry not found'; raise e

    response = make_response( jsonify( {'ac': uniprot.ac, 'entry_name': uniprot.entry_name, 'family': uniprot.family} ) )
    return response
