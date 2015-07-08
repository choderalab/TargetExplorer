from flask import abort, jsonify, request, make_response, current_app
import re
from datetime import timedelta
from functools import update_wrapper
from targetexplorer.flaskapp import app, db, models
config = app.config

# ======
# HTTP access control decorator - for cross-origin resource sharing (CORS)
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
# Get list of all db entries
# ======

# Examples:
# http://.../[DB_NAME]DBAPI/listall

@app.route('/%s/listall' % config.dbapi_name, methods = ['GET'])
@crossdomain(origin='*', headers=["Origin", "X-Requested-With", "Content-Type", "Accept"])
def listall():
    # = Get safe crawl number =
    crawldata = models.CrawlData.query.first()
    safe_crawl_number = crawldata.safe_crawl_number

    # = Query the UniProt table =
    uniprot_values = [values for values in models.UniProt.query.filter_by(crawl_number=safe_crawl_number).values(models.UniProt.ac, models.UniProt.entry_name)]

    # = Construct the data structure for holding the results, to be returned as JSON =
    results_obj = {
        'listall': [],
    }

    # = Add info from other tables =
    for entry_data in uniprot_values:
        entry_obj = {
            'ac': entry_data[0],
            'entry_name': entry_data[1],
        }
        results_obj['listall'].append(entry_obj)

    # = Return data in JSON format =
    response = make_response( jsonify(results_obj) )
    return response


# ======
# Get individual db entry
# ======

# Examples:
# http://.../[DB_NAME]DBAPI/entry?ac=P00519
# http://.../[DB_NAME]DBAPI/entry?entry_name=P00519

@app.route('/%s/entry' % config.dbapi_name, methods = ['GET'])
@crossdomain(origin='*', headers=["Origin", "X-Requested-With", "Content-Type", "Accept"])
def get_dbentry():
    # = Get safe crawl number =
    crawldata = models.CrawlData.query.first()
    safe_crawl_number = crawldata.safe_crawl_number

    for request_field in ['ac', 'entry_name']:
        if request_field in request.args:
            break

    request_value = request.args.get(request_field)
    if request_field == 'ac':
        # TODO UniProt AC format will be extended some time after June 11 2014!
        if not (re.match('[A-N,R-Z][0-9][A-Z][A-Z,0-9][A-Z,0-9][0-9]', request_value) or re.match('[O,P,Q][0-9][A-Z,0-9][A-Z,0-9][A-Z,0-9][0-9]', request_value)):
            abort(404)

    # = Search the UniProt table using the query AC =
    model_request_field = getattr(models.UniProt, request_field)
    uniprot = models.UniProt.query.filter(model_request_field==request_value, models.UniProt.crawl_number==safe_crawl_number).first()
    try: assert uniprot != None
    except AssertionError as e: e.message = 'Database entry not found'; raise e

    # = Retrieve the corresponding DBEntry rows =
    dbentry = db.session.query(models.DBEntry).filter_by(id=uniprot.dbentry_id, crawl_number=safe_crawl_number).first()

    # = Construct the data structure for holding the results, to be returned as JSON =
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
        'npubs': dbentry.npubs,
        'nbioassays': dbentry.nbioassays,
    }

    # = Add info from other tables =
    # PDB
    for pdb in dbentry.pdbs:
        target_obj['pdb'].append({'pdbid': pdb.pdbid})

    # HGNC
    for entry in dbentry.hgnc_entries:
        target_obj['hgnc'].append({'gene_id': entry.gene_id, 'approved_symbol': entry.approved_symbol})

    # Ensembl Gene
    for entry in dbentry.ensembl_genes:
        target_obj['ensembl_gene'].append({'gene_id': entry.gene_id})

    # NCBI Gene
    for entry in dbentry.ncbi_gene_entries:
        target_obj['ncbi_gene'].append({'gene_id': entry.gene_id})

    # = Return data in JSON format =
    response = make_response(jsonify(target_obj))
    return response


# ======
# Get multiple database entries given a query string
# ======

# Examples:
# http://.../[DB_NAME]DBAPI/search?query=family="TK" AND db_target_rank<300&return="domain_seqs"

# Example SQLAlchemy filter syntax:
# 'family is null AND species="Human"'

@app.route('/%s/search' % config.dbapi_name, methods = ['GET'])
@crossdomain(origin='*', headers=["Origin", "X-Requested-With", "Content-Type", "Accept"])
def query_db():
    frontend_query_string = request.args.get('query') # expecting SQLAlchemy syntax (wtih frontend-style field names)
    return_fields = request.args.get('return') if 'return' in request.args else []
    if type(return_fields) == str:
        return_fields = [return_fields]

    # = Get safe crawl number =
    crawldata = models.CrawlData.query.first()
    safe_crawl_number = crawldata.safe_crawl_number

    # = Use the query string to query the db =
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
    # (note that we only need to filter by crawl_number for the DBEntry rows)
    query = db.session.query(models.DBEntry).filter_by(crawl_number=safe_crawl_number)
    for query_table_name in query_tables:
        if query_table_name != 'DBEntry':
            query_table = models.__dict__[query_table_name]
            if not hasattr(query_table, 'dbentry_id'):
                raise Exception, 'ERROR: Cannot filter on table %s (no relationship to table DBEntry)' % query_table_name
            query = query.join(query_table, models.DBEntry.id==query_table.dbentry_id)

    # Use the query string to filter DBEntry rows
    results = query.filter(sql_query_string)

    # = Construct the data structure for holding the results, to be returned as JSON =
    targets_obj = {'results': []}

    for db_entry in results:
        uniprot = db.session.query(models.UniProt).filter_by(dbentry_id=db_entry.id).first()
        domain_targetids = [domain_row.targetid for domain_row in uniprot.domains]
        target_obj = {
            'ac': uniprot.ac,
            'entry_name': uniprot.entry_name,
            'family': uniprot.family,
            'npdbs': db_entry.npdbs,
            'npubs': db_entry.npubs,
            'nbioassays': db_entry.nbioassays,
            'domains': domain_targetids,
        }

        # Optional additional data
        if 'seqs' in return_fields:
            canon_isoform = uniprot.isoforms.filter_by(canonical=True).first()
            target_obj['sequence'] = canon_isoform.sequence
        if 'domain_seqs' in return_fields:
            domain_data = [{'targetid': domain_row.targetid, 'sequence': domain_row.sequence} for domain_row in uniprot.domains]
            target_obj['domains'] = domain_data
        if 'pdb_data' in return_fields:
            pdb_data = []
            for pdb_row in db_entry.pdbs:
                pdbchain_data = [{'chainid': pdb_chain_row.chain_id, 'domainid': pdb_chain_row.domain_id, 'seq_begin': pdb_chain_row.begin, 'seq_end': pdb_chain_row.end} for pdb_chain_row in pdb_row.chains]
                pdb_data.append({'pdbid': pdb_row.pdbid, 'pdbchains': pdbchain_data})

            target_obj['pdbs'] = pdb_data

        targets_obj['results'].append(target_obj)

    # = Return data in JSON format =
    response = make_response( jsonify(targets_obj) )
    return response

# ======
# Get database metadata
# ======

# Examples:
# http://.../[DB_NAME]DBAPI/get_metadata

@app.route('/%s/get_metadata' % config.dbapi_name, methods = ['GET'])
@crossdomain(origin='*', headers=["Origin", "X-Requested-With", "Content-Type", "Accept"])
def get_metadata():
    import project_config
    results_obj = {
        'uniprot_query_string': project_config.uniprot_query_string,
        'uniprot_domain_regex': project_config.uniprot_domain_regex,
    }

    # = Return data in JSON format =
    response = make_response(jsonify(results_obj))
    return response