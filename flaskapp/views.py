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
# URL query handler
# ======

# Examples:
# http://ec2-54-227-62-182.compute-1.amazonaws.com/kinomeDBAPI/P00519
# http://ec2-54-227-62-182.compute-1.amazonaws.com/kinomeDBAPI/?query=family:"TK"

@app.route('/<string:leadingpath>/<string:query_string>', methods = ['GET'])
@crossdomain(origin='*', headers=["Origin", "X-Requested-With", "Content-Type", "Accept"])
def query_handler(leadingpath, query_string):
    # note: leadingpath is ignored

    db_response = None

    # If query string is an individual UniProt AC
    # XXX TODO UniProt AC format will be extended some time after June 11 2014!
    if re.match('[A-N,R-Z][0-9][A-Z][A-Z,0-9][A-Z,0-9][0-9]', query_string) or re.match('[O,P,Q][0-9][A-Z,0-9][A-Z,0-9][A-Z,0-9][0-9]', query_string):
        ac = query_string
        db_response = get_dbentry(ac)

    try: assert db_response != None
    except AssertionError as e: e.message = 'Invalid query string'; raise e

    response = make_response( jsonify(db_response) )
    return response


# ======
# Get data from database for individual target
# ======

def get_dbentry(ac):
    uniprot = db.session.query(models.UniProt).filter_by(ac=ac).first()
    print ac
    try: assert uniprot != None
    except AssertionError as e: e.message = 'Database entry not found'; raise e

    dbentry = db.session.query(models.DBEntry).filter_by(id=uniprot.dbentry_id).first()

    target_obj = {
        'uniprot': {
            'ac': uniprot.ac,
            'entry_name': uniprot.entry_name,
            'family': uniprot.family,
        },
        'pdb': [],
        'hgnc': [],
        'ensembl_gene': [],
        'ncbi_gene': [],
    }

    # PDB
    for pdb in dbentry.pdbs:
        target_obj['pdb'].append({'pdbid': pdb.pdbid})

    # HGNC
    for entry in dbentry.hgnc_entries:
        target_obj['hgnc'].append({'gene_id': entry.gene_id, 'approved_symbol': entry.approved_symbol})

    # Ensembl Gene
    for entry in dbentry.ensembl_gene_entries:
        target_obj['ensembl_gene'].append({'gene_id': entry.gene_id})

    # NCBI Gene
    for entry in dbentry.ncbi_gene_entries:
        target_obj['ncbi_gene'].append({'gene_id': entry.gene_id})

    return target_obj

