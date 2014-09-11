#!/usr/bin/env python
#
# For each PDB in database-stage.xml, download the SIFTS residue-mapping .xml file
# Add experimental and resolved sequences to the DB (numbered according to the UniProt sequence)
# Also add alignments of these sequences against the UniProt sequence
#
# Daniel L. Parton <partond@mskcc.org> - 10 Apr 2013
#
# Perhaps get DSSP info from here: http://www.rcsb.org/pdb/rest/das/pdbchainfeatures/features?segment=5pti.A
#

#==============================================================================
# Imports
#==============================================================================

import sys
import os
import gzip
import re
import datetime
import copy
import urllib2
import traceback
import argparse
from lxml import etree
import targetexplorer
from targetexplorer.flaskapp import models, db
import Bio.PDB
import Bio.Data.SCOPData
from multiprocessing import Pool

# =================
# Parameters
# =================

external_data_dir = os.path.join('external-data')
local_pdb_dir = os.path.join(external_data_dir, 'PDB')
local_sifts_dir = os.path.join(external_data_dir, 'SIFTS')
if not os.path.exists(local_pdb_dir):
    os.mkdir(local_pdb_dir)
if not os.path.exists(local_sifts_dir):
    os.mkdir(local_sifts_dir)

argparser = argparse.ArgumentParser(description='Gather PDB data')
args = argparser.parse_args()

now = datetime.datetime.utcnow()

# get current crawl number
crawldata_row = models.CrawlData.query.first()
current_crawl_number = crawldata_row.current_crawl_number
print 'Current crawl number: %d' % current_crawl_number

parser = etree.XMLParser(remove_blank_text=True)

verbose = False




# #==============================================================================
# # Parameters
# #==============================================================================
#
# if '-stage' in sys.argv:
#     run_mode = 'stage'
# elif '-dev' in sys.argv:
#     run_mode = 'dev'
# else:
#     run_mode = 'nowrite'
#
# print 'Running in mode: %s' % run_mode
#
# database_dir = 'database'
# external_data_dir = 'external-data'
# local_pdb_dir = os.path.join(external_data_dir, 'PDB')
# local_sifts_dir = os.path.join(external_data_dir, 'SIFTS')
# if not os.path.exists(local_pdb_dir):
#     os.mkdir(local_pdb_dir)
# if not os.path.exists(local_sifts_dir):
#     os.mkdir(local_sifts_dir)
#
# DBstage_filepath = os.path.join(database_dir, 'database-stage.xml')
# if not os.path.exists(DBstage_filepath):
#     raise Exception, '%s not found.' % DBstage_filepath
#
# if run_mode != 'nowrite':
#     DB_out_filename = 'database-%(run_mode)s.xml' % vars()
#     DB_out_filepath = os.path.join(database_dir, DB_out_filename)
#
# verbose = False
#
# now = datetime.datetime.utcnow()
# datestamp = now.strftime(targetexplorer.DB.datestamp_format_string)
#
# parser = etree.XMLParser(remove_blank_text=True)

#==============================================================================
# Definitions
#==============================================================================

class Extract_PDB_Results(object):
    '''
    Container class returned by the method gather_pdb
    '''
    def __init__(self, UniProt_entry_name=None):
        self.entry_name = UniProt_entry_name
        self.structures = []
    def add_structure_results(self, structures):
        self.structures.append(structures)

class Structure_Data(object):
    '''
    Class for communicating PDB structure data
    '''
    def __init__(self, structureID=None, exception_message=None):
        self.structureID = structureID
        self.exception_message = exception_message
        self.chains = []
        self.expression_data = None
    def add_chain_results(self, chain_results):
        self.chains.append(chain_results)
    def add_expression_data(self, expression_data):
        self.expression_data = expression_data
    def add_pdbcompoundID(self, pdbcompoundID):
        self.pdbcompoundID = pdbcompoundID

class Chain_Data(object):
    '''
    Class for communicating PDB chain data
    '''
    def __init__(self, chainID=None, experimental_sequence=None, experimental_sequence_aln=None, experimental_sequence_aln_conflicts=None, observed_sequence_aln_exp=None, observed_sequence_aln=None, ss_aln=None, exception_message=None):
        self.chainID = chainID
        self.experimental_sequence=experimental_sequence
        self.experimental_sequence_aln=experimental_sequence_aln
        self.experimental_sequence_aln_conflicts=experimental_sequence_aln_conflicts
        self.observed_sequence_aln_exp=observed_sequence_aln_exp
        self.observed_sequence_aln=observed_sequence_aln
        self.ss_aln=ss_aln
        self.exception_message=exception_message

def extract_pdb_data(pdb_dict):
    '''Extract data for a single PDB structure
    '''
    # pdbid, ac, entry_name, seq, chain_ids
    pdb_row_id = pdb_dict['pdb_row_id']
    pdbid = pdb_dict['pdbid']
    ac = pdb_dict['ac']
    entry_name = pdb_dict['entry_name']
    seq = pdb_dict['seq']
    chain_data = pdb_dict['chain_data']

    # if entry_name != 'MLKL_HUMAN':
    #     return None
    #
    # if pdbid != '2ITN':
    #     return None

    # ========
    # Get PDB and SIFTS files
    # PDB files are used to extract expression system metadata
    # SIFTS files are used to extract sequence data
    # ========

    # TODO define this via project metadata .yaml file.
    structure_paths = ['/Users/partond/tmp/kinome-MSMSeeder/structures/pdb', '/Users/partond/tmp/kinome-MSMSeeder/structures/sifts']

    local_pdb_filepath = os.path.join('external-data', 'PDB', pdbid + '.pdb.gz')
    local_sifts_filepath = os.path.join('external-data', 'SIFTS', pdbid + '.xml.gz')

    # Check if PDB file/symlink already exists and is not empty
    search_for_pdb = True
    if os.path.exists(local_pdb_filepath):
        if os.path.getsize(local_pdb_filepath) > 0:
            search_for_pdb = False

    # If not, search any user-defined paths and create a symlink if found
    if search_for_pdb:
        for structure_dir in structure_paths:
            pdb_filepath = os.path.join(structure_dir, pdbid + '.pdb.gz')
            if os.path.exists(pdb_filepath):
                if os.path.getsize(pdb_filepath) > 0:
                    if os.path.exists(local_pdb_filepath):
                        os.remove(local_pdb_filepath)
                    os.symlink(pdb_filepath, local_pdb_filepath)
                    break

        # If still not found, download the PDB file
        if not os.path.exists(local_pdb_filepath):
            print 'Downloading PDB file and saving as:', local_pdb_filepath
            page = targetexplorer.PDB.retrieve_pdb(pdbid, compressed='yes')
            # download and write compressed file
            with open(local_pdb_filepath, 'wb') as local_pdb_file:
                local_pdb_file.write(page)

    # Check if SIFTS file already exists and is not empty
    search_for_sifts = True
    if os.path.exists(local_sifts_filepath):
        if os.path.getsize(local_sifts_filepath) > 0:
            search_for_sifts = False

    # If not, search any user-defined paths and create a symlink if found
    if search_for_sifts:
        for structure_dir in structure_paths:
            sifts_filepath = os.path.join(structure_dir, pdbid + '.xml.gz')
            if os.path.exists(sifts_filepath):
                if os.path.getsize(sifts_filepath) > 0:
                    if os.path.exists(local_sifts_filepath):
                        os.remove(local_sifts_filepath)
                    os.symlink(sifts_filepath, local_sifts_filepath)
                    break

        # If still not found, download the SIFTS XML file
        if not os.path.exists(local_sifts_filepath):
            print 'Downloading SIFTS file (compressed) and saving as:', local_sifts_filepath
            try:
                page = targetexplorer.PDB.retrieve_sifts(pdbid)
            except urllib2.URLError as urlerror:
                if urlerror.reason == 'ftp error: [Errno ftp error] 550 Failed to change directory.':
                    # Check the PDB file has definitely been downloaded. If so, then the problem is probably that the SIFTS people have not yet created the file for this PDB entry, or they have not added it to their server yet.
                    if os.path.exists(local_pdb_filepath):
                        # In this case, just add a message telling the script to delete this PDB structure from the DB. The continue clause skips to the end of the function.
                        print '%s SIFTS file could not be downloaded - this PDB entry will be deleted from the DB' % pdbid
                        # structure_results_obj = Structure_Data(structureID=pdbid)
                        # structure_results_obj.exception_message='DELETE_ME - SIFTS file could not be downloaded'
                        # results_obj.add_structure_results(structure_results_obj)
                        return [pdbid, 'DELETE']
                    else:
                        raise urlerror
                else:
                    raise urlerror

            with gzip.open(local_sifts_filepath, 'wb') as local_sifts_file:
                local_sifts_file.write(page)

    # ======
    # From PDB file, get EXPRESSION_SYSTEM and related fields, using Bio.PDB.PDBParser
    # ======

    db_chain_ids_lower = [chain_dict['chain_id'].lower() for chain_dict in chain_data]

    pdbparser = Bio.PDB.PDBParser(QUIET=True)
    with gzip.open(local_pdb_filepath) as local_pdb_file:
        pdbdata = pdbparser.get_structure(pdbid, local_pdb_file)
    pdbheader = pdbparser.get_header()
    # Bio PDB compound structure: {'compound': {'1': {'chain': 'a, b'}}}
    pdb_compounds = pdbheader['compound']
    for pdb_compound_id in pdb_compounds.keys():
        try:
            for pdb_chain_id in pdb_compounds[pdb_compound_id]['chain'].split(', '):
                if pdb_chain_id in db_chain_ids_lower:
                    matching_pdb_compound_id = pdb_compound_id
                    break
        except Exception as e:
            print 'ERROR for entry %s PDB %s. PDB header dict as parsed by BioPython follows:' % (entry_name, pdbid)
            print pdbheader
            print traceback.format_exc()
            raise e

    expression_data = {}
    # Bio PDB source structure: {'source': {'1': {'expression_system': 'escherichia coli'}}}
    pdbexpression_data = pdbheader['source'][matching_pdb_compound_id]
    for key in pdbexpression_data.keys():
        if key[0:10] == 'expression':
            # Make expression data upper-case again. I think it looks better for single-case text.
            expression_data[key.upper()] = pdbexpression_data[key].upper()
            # expression_data_obj = models.PDBExpressionData(expression_data_type=key.upper(), expression_data_value=pdbexpression_data[key].upper(), pdb=pdb_row)
            # db.session.add(expression_data_obj)

    # ======
    # Iterate through chains in PDBRow and extract sequence data from SIFTS file, and add to database
    # ======

    results = {'pdb_row_id': pdb_row_id, 'expression_data': expression_data, 'chain_dicts': {}}

    for chain_dict in chain_data:
        chain_row_id = chain_dict['chain_row_id']
        chain_id = chain_dict['chain_id']
        if verbose: print entry_name, ac, pdbid, chain_id
        pdb_chain_dict = targetexplorer.PDB.extract_sifts_seq(local_sifts_filepath, ac, entry_name, pdbid, chain_id, seq)
        results['chain_dicts'][chain_row_id] = pdb_chain_dict

    return results




# ====================
# Main
# ====================

if __name__ == '__main__':

    # first get a list of PDB rows from the db
    db_pdb_rows = models.PDB.query.filter_by(crawl_number=current_crawl_number).all()

    # prepare a dict for each PDB row to be passed to the function extract_pdb_data
    def dictify_pdb_row(pdb_row):
        dbentry = models.DBEntry.query.filter_by(id=pdb_row.dbentry_id).first()
        uniprot_row = dbentry.uniprot.first()
        canon_isoform_row = dbentry.uniprotisoforms.filter_by(canonical=True).first()
        chain_data = [{'chain_row_id': chain_row.id, 'chain_id': chain_row.chain_id} for chain_row in pdb_row.chains]
        pdb_dict = {
            'pdb_row_id': pdb_row.id,
            'pdbid': pdb_row.pdbid,
            'ac': uniprot_row.ac,
            'entry_name': uniprot_row.entry_name,
            'seq': canon_isoform_row.sequence,
            'chain_data': chain_data,
        }
        return pdb_dict

    db_pdb_dicts = map(dictify_pdb_row, db_pdb_rows)

    # Use multiprocessor pool to retrieve various data for each PDB
    pool = Pool()
    # results = pool.map(extract_pdb_data, [(db_pdb_row.id, current_crawl_number) for db_pdb_row in db_pdb_rows])
    results = pool.map(extract_pdb_data, db_pdb_dicts)

    for pdb_results in results:
        pdb_row_id = pdb_results['pdb_row_id']
        pdb_row = models.PDB.query.filter_by(id=pdb_row_id).first()
        for chain_row_id in pdb_results['chain_dicts']:
            chain_row = models.PDBChain.query.filter_by(id=chain_row_id).first()
            chain_dict = pdb_results['chain_dicts'][chain_row_id]
            for key in chain_dict:
                if key in ['chain_id', 'exception_message']:
                    continue
                setattr(chain_row, key, chain_dict[key])

            # Delete PDB structure and chain entries with @DELETE_ME attrib. These are cases where the sifts_uniprotAC does not match the uniprotAC in DB_root (derived from the UniProt entry by gather-uniprot.py), or where more than 90% of the experimental sequence is unobserved
            if chain_dict['exception_message'] == 'DELETE_ME':
                db.session.delete(chain_row)

        if pdb_row.chains.count() == 0:
            db.session.delete(pdb_row)

# ==============================================================
# Update db PDB datestamp
# ==============================================================

current_crawl_datestamp_row = models.DateStamps.query.filter_by(crawl_number=current_crawl_number).first()
current_crawl_datestamp_row.pdb_datestamp = now
db.session.commit()
print 'Done.'
