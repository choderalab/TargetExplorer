from flask import abort, jsonify, request, make_response, current_app
import re
from datetime import timedelta
from functools import update_wrapper
from flaskapp import app, db, models
import flaskapp_config

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
# Get individual db entry
# ======

# Examples:
# http://.../[DB_NAME]DBAPI/entry?ac=P00519

@app.route('/%s/entry' % flaskapp_config.dbapi_name, methods = ['GET'])
@crossdomain(origin='*', headers=["Origin", "X-Requested-With", "Content-Type", "Accept"])
def get_dbentry():
    # note: leadingpath is ignored

    ac = request.args.get('ac')
    if ac == None:
        abort(404)
    # TODO UniProt AC format will be extended some time after June 11 2014!
    elif not (re.match('[A-N,R-Z][0-9][A-Z][A-Z,0-9][A-Z,0-9][0-9]', ac) or re.match('[O,P,Q][0-9][A-Z,0-9][A-Z,0-9][A-Z,0-9][0-9]', ac)):
        abort(404)

    uniprot = db.session.query(models.UniProt).filter_by(ac=ac).first()
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

    response = make_response( jsonify(target_obj) )
    return response


# ======
# Get multiple database entries given a query string
# ======

# Examples:
# http://.../[DB_NAME]DBAPI/search?query=family=TK AND db_target_rank<300

@app.route('/%s/search' % flaskapp_config.dbapi_name, methods = ['GET'])
@crossdomain(origin='*', headers=["Origin", "X-Requested-With", "Content-Type", "Accept"])
def query_db():
    frontend_query_string = request.args.get('query') # expecting SQLAlchemy syntax (wtih frontend-style field names)

    # Convert the frontend data fields to backend identifiers '[table].[column]'
    # And determine which tables will need to be queried
    sql_query_string = frontend_query_string
    query_tables = []
    for frontend_field_name, backend_data_list in models.frontend2backend_mappings.iteritems():
        sql_query_string = sql_query_string.replace(frontend_field_name, '.'.join(backend_data_list))
        if frontend_field_name in frontend_query_string:
            query_tables.append(backend_data_list[0])
    query_tables = set(query_tables)

    # Start with the DBEntry table, then carry out SQL joins with the other tables
    query = db.session.query(models.DBEntry)
    for query_table_name in query_tables:
        if query_table_name != 'DBEntry':
            query_table = models.__dict__[query_table_name]
            if not hasattr(query_table, 'dbentry_id'):
                raise Exception, 'ERROR: Cannot filter on table %s (no relationship to table DBEntry)' % query_table_name
            query = query.join(query_table, models.DBEntry.id==query_table.dbentry_id)

    # Use the query string to filter DBEntry rows
    results = query.filter(sql_query_string)

    # Build results object
    targets_obj = {'results': []}

    for db_entry in results:
        uniprot = db.session.query(models.UniProt).filter_by(dbentry_id=db_entry.id).first()
        target_obj = {
            'ac': uniprot.ac,
            'entry_name': uniprot.entry_name,
            'family': uniprot.family,
            'npdbs': db_entry.pdbs.count(),
        }
        targets_obj['results'].append(target_obj)

    # Return results in JSON format
    response = make_response( jsonify(targets_obj) )
    return response
