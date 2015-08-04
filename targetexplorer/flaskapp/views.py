import re
from flask import abort, jsonify, request, make_response
from targetexplorer.flaskapp import app, db, models
from targetexplorer.core import read_project_config
from targetexplorer.flaskapp.webapi_utils import crossdomain

project_config = read_project_config()
dbapi_name = project_config['dbapi_name']


# ======
# Get database metadata
# ======

# Examples:
# http://.../[DB_NAME]DBAPI/get_metadata

@app.route('/%s/get_metadata' % dbapi_name, methods=['GET'])
@crossdomain(origin='*', headers=["Origin", "X-Requested-With", "Content-Type", "Accept"])
def get_metadata():
    project_config = read_project_config()
    results_obj = {
        'uniprot_query': project_config['uniprot_query'],
        'uniprot_domain_regex': project_config['uniprot_domain_regex'],
    }

    # = Return data in JSON format =
    response = make_response(jsonify(results_obj))
    return response


# ======
# Get list of all db entries
# ======

# Examples:
# http://.../[DB_NAME]DBAPI/listall

@app.route('/%s/listall' % dbapi_name, methods=['GET'])
@crossdomain(origin='*', headers=["Origin", "X-Requested-With", "Content-Type", "Accept"])
def listall():
    crawldata = models.CrawlData.query.first()
    safe_crawl_number = crawldata.safe_crawl_number

    uniprot_values = [
        values for values
        in models.UniProtEntry.query.filter_by(
            crawl_number=safe_crawl_number
        ).values(models.UniProtEntry.ac, models.UniProtEntry.entry_name)
    ]

    results_obj = {
        'listall': [],
    }

    for entry_data in uniprot_values:
        entry_obj = {
            'ac': entry_data[0],
            'entry_name': entry_data[1],
        }
        results_obj['listall'].append(entry_obj)

    response = make_response(jsonify(results_obj))
    return response


# ======
# Get individual db entry
# ======

# Examples:
# http://.../[DB_NAME]DBAPI/entry?ac=P00519
# http://.../[DB_NAME]DBAPI/entry?entry_name=ABL1_HUMAN

@app.route('/%s/entry' % dbapi_name, methods=['GET'])
@crossdomain(origin='*', headers=["Origin", "X-Requested-With", "Content-Type", "Accept"])
def get_db_entry():
    crawldata = models.CrawlData.query.first()
    safe_crawl_number = crawldata.safe_crawl_number

    for request_field in ['ac', 'entry_name']:
        if request_field in request.args:
            break

    request_value = request.args.get(request_field)
    if request_field == 'ac':
        # TODO UniProt AC format will be extended some time after June 11 2014!
        if not (re.match('[A-N,R-Z][0-9][A-Z][A-Z,0-9][A-Z,0-9][0-9]', request_value) or
                re.match('[O,P,Q][0-9][A-Z,0-9][A-Z,0-9][A-Z,0-9][0-9]', request_value)
                ):
            abort(404)

    # = Search the UniProt table =
    # e.g. uniprot_row = models.UniProtEntry.query.filter(models.UniProtEntry.ac == 'P00519', models.UniProtEntry.crawl_number == 0).first()
    model_request_field = getattr(models.UniProtEntry, request_field)
    uniprot_row = models.UniProtEntry.query.filter(
        model_request_field == request_value, models.UniProtEntry.crawl_number == safe_crawl_number
    ).first()

    if uniprot_row is None:
        raise UniProtEntryNotFoundError('{}={}'.format(request_field, request_value))

    # = Retrieve the corresponding DBEntry row =
    db_entry = uniprot_row.db_entry

    # = Construct the data structure for holding the results, to be returned as JSON =
    db_entry_obj = {
        'uniprot': {
            'ac': uniprot_row.ac,
            'entry_name': uniprot_row.entry_name,
            'family': uniprot_row.family,
        },
        'pdb': [],
        'hgnc': [],
        'ensembl_gene': [],
        'ncbi_gene': [],
        'npubs': db_entry.npubs,
        'nbioassays': db_entry.nbioassays,
    }

    # = Add info from other tables =
    # PDB
    for pdb in db_entry.pdbs:
        db_entry_obj['pdb'].append({'pdb_id': pdb.pdb_id})

    # HGNC
    for entry in db_entry.hgnc_entries:
        db_entry_obj['hgnc'].append(
            {'gene_id': entry.gene_id, 'approved_symbol': entry.approved_symbol}
        )

    # Ensembl Gene
    for entry in db_entry.ensembl_genes:
        db_entry_obj['ensembl_gene'].append({'gene_id': entry.gene_id})

    # NCBI Gene
    for entry in db_entry.ncbi_gene_entries:
        db_entry_obj['ncbi_gene'].append({'gene_id': entry.gene_id})

    # = Return data in JSON format =
    response = make_response(jsonify(db_entry_obj))
    return response


# ======
# Get multiple database entries given a query string
# ======

# Examples:
# http://.../[DB_NAME]DBAPI/search?query=family="TK" AND npubs>0&return="domain_seqs"
# http://.../[DB_NAME]DBAPI/search?query=family="TK" AND db_target_rank<300&return="domain_seqs"

# Example SQLAlchemy filter syntax:
# 'family is null AND species="Human"'

@app.route('/%s/search' % dbapi_name, methods=['GET'])
@crossdomain(origin='*', headers=["Origin", "X-Requested-With", "Content-Type", "Accept"])
def query_db():
    frontend_query_string = request.args.get('query') # expecting SQLAlchemy syntax (wtih frontend-style field names)
    return_fields = request.args.get('return') if 'return' in request.args else []
    if type(return_fields) == str:
        return_fields = [return_fields]

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
            # query_table = models.__dict__[query_table_name]
            query_table = getattr(models, query_table_name)
            if not hasattr(query_table, 'db_entry_id'):
                raise Exception(
                    'ERROR: Cannot filter on table {} (no relationship to table DBEntry)'.format(
                        query_table_name
                    )
                )
            query = query.join(query_table, models.DBEntry.id == query_table.db_entry_id)

    # Use the query string to filter DBEntry rows
    results = query.filter(sql_query_string)

    # = Construct the data structure for holding the results, to be returned as JSON =
    targets_obj = {'results': []}

    for db_entry in results:
        uniprot = db.session.query(models.UniProtEntry).filter_by(db_entry_id=db_entry.id).first()
        domain_descriptions = [domain_row.description for domain_row in uniprot.domains]
        target_domain_ids = [
            domain_row.target_id for domain_row in uniprot.domains if domain_row.is_target_domain
        ]
        target_obj = {
            'ac': uniprot.ac,
            'entry_name': uniprot.entry_name,
            'family': uniprot.family,
            'npdbs': db_entry.npdbs,
            'npubs': db_entry.npubs,
            'nbioassays': db_entry.nbioassays,
            'domains': domain_descriptions,
            'target_domains': target_domain_ids,
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
    response = make_response(jsonify(targets_obj))
    return response


# ======
# error handlers
# ======


class UniProtEntryNotFoundError(Exception):
    pass


@app.errorhandler(UniProtEntryNotFoundError)
def return_uniprot_entry_not_found_error(error):
    return make_response(
        jsonify({
            'error': 'UniProt entry not found',
            'message': error.message,
        })
    )


@app.errorhandler(404)
def not_found(error):
    return make_response(jsonify({'error': 'Not found'}), 404)
