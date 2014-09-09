#!/usr/bin/env python
#
# Add bioassay data from BindingDB to database.
#
# Daniel L. Parton <partond@mskcc.org> - 25 Oct 2013
#

# =================
# Imports
# =================

import sys
import os
import subprocess
import datetime
import argparse
from lxml import etree
import targetexplorer
from targetexplorer.flaskapp import models, db
from multiprocessing import Pool

# =================
# Parameters
# =================

external_data_dir = os.path.join('external-data')
bindingdb_data_dir = os.path.join(external_data_dir, 'BindingDB')

if not os.path.exists(bindingdb_data_dir):
    os.mkdir(bindingdb_data_dir)

bindingdb_all_data_filepath = os.path.join(bindingdb_data_dir, 'BindingDB_All.tab')

argparser = argparse.ArgumentParser(description='Gather BindingDB data')
argparser.add_argument('--use_existing_bindingdb_data', help='Do not download a new BindingDB_All.tab file. Only works if an existing file is present.', action='store_true', default=False)

argparser.add_argument('--grep_path', help='Provide explicit path for grep. Otherwise, searches in the usual places.', type=str, default='grep')
args = argparser.parse_args()

now = datetime.datetime.utcnow()

# get current crawl number
crawldata_row = models.CrawlData.query.first()
current_crawl_number = crawldata_row.current_crawl_number
print 'Current crawl number: %d' % current_crawl_number

parser = etree.XMLParser(remove_blank_text=True)

# =================
# Download new BindingDB data
# =================

# Unless args.use_existing_bindingdb_data is set to True, retrieve a new file from BindingDB
if os.path.exists(bindingdb_all_data_filepath) and args.use_existing_bindingdb_data:
    print 'BindingDB data file found at:', bindingdb_all_data_filepath
else:
    print 'Retrieving new BindingDB data file from BindingDB server...'
    targetexplorer.BindingDB.retrieve_all_BindingDB_data(bindingdb_all_data_filepath, decompress=False)



# clab.DB.retrieve_external_data(DB_script_ID, forcedl=force_bindingdb_all_download)

# =================
# Use grep to extract information for each entry in the DB
# =================

# Overview:
#   * first use grep -E to extract all lines matching any kinase uniprot AC
#   * using multiprocessing for each kinase:
#       * grep through to create data file for each kinase
#       * parse each kinase data file and add to kinDB
#       * delete individual kinase data files

db_uniprot_acs = [value_tuple[0] for value_tuple in models.UniProt.query.filter_by(crawl_number=current_crawl_number).values(models.UniProt.ac)]

# db_uniprot_acs = ['P00519']




regex_search = '|'.join(db_uniprot_acs)
bindingdb_matches_filepath = os.path.join(bindingdb_data_dir, 'tmp-bindingdb-matches.tab')

# TODO use proper logic here
if True:
    print 'First pass through BindingDB data with grep... - extracting all lines which match UniProt ACs in the DB (performance can be variable - testing has indicated approx. 1 hr on a 2012 MBP, but only 2 mins or so on a Linux machine; probably due to different GNU grep versions)'
    with open(bindingdb_all_data_filepath, 'r') as bindingdb_all_data_file:
        #subprocess.call('%(args.grep_path)s -E "%(regex_search)s" %(bindingdb_all_data_filepath)s > %(bindingdb_matches_filepath)s' % vars(), shell=True)
        subprocess.call('%s -E "%s" %s > %s' % (args.grep_path, regex_search, bindingdb_all_data_filepath, bindingdb_matches_filepath), shell=True)

def extract_bindingdb(input_data):
    AC, grep_path, bindingdb_matches_filepath = input_data
    print AC
    DB_entry_bindingdb_data_filepath = os.path.join(bindingdb_data_dir, AC + '.tab')
    subprocess.call('%s -E "%s" %s > %s' % (grep_path, AC, bindingdb_matches_filepath, DB_entry_bindingdb_data_filepath), shell=True)
    bioassays_data = []
    with open(DB_entry_bindingdb_data_filepath, 'r') as bindingdb_file:
        for line in bindingdb_file:
            returnedACs = targetexplorer.BindingDB.get_ACs(line)
            for returnedAC in returnedACs:
                if returnedAC == AC:
                    bioassay_data = targetexplorer.BindingDB.get_bioassay_data(line) # returns a dict
                    bioassays_data.append( bioassay_data )
    os.remove(DB_entry_bindingdb_data_filepath)
    return (AC, bioassays_data)

if __name__ == '__main__':
    pool = Pool()
    input_data = [(AC, args.grep_path, bindingdb_matches_filepath) for AC in db_uniprot_acs]
    results = pool.map(extract_bindingdb, input_data)

for result in results:
    ac, bioassays_data = result
    db_uniprot_row = models.UniProt.query.filter_by(crawl_number=current_crawl_number, ac=ac).first()
    dbentry_id = db_uniprot_row.dbentry_id
    dbentry_row = models.DBEntry.query.filter_by(id=dbentry_id).first()
    dbentry_row.nbioassays = len(bioassays_data)

    for bioassay_data in bioassays_data:
        bindingdb_bioassay_obj = models.BindingDBBioassay(
            crawl_number=current_crawl_number,
            bindingdb_source=bioassay_data['BindingDB_source'],
            doi=bioassay_data['DOI'],
            pmid=bioassay_data['PMID'],
            temperature=bioassay_data['temperature'],
            ph=bioassay_data['pH'],
            target_name=bioassay_data['target_name'],
            ligand_bindingdb_id=bioassay_data['ligand_BindingDB_ID'],
            ligand_chembl_id=bioassay_data['ligand_ChEMBL_ID'],
            ligand_smiles_string=bioassay_data['ligand_SMILES_string'],
            ligand_zinc_id=bioassay_data['ligand_zinc_id'],

            ki=bioassay_data['Ki'],
            ic50=bioassay_data['IC50'],
            kd=bioassay_data['Kd'],
            ec50=bioassay_data['EC50'],
            kon=bioassay_data['kon'],
            koff=bioassay_data['koff'],

            dbentry=dbentry_row,
        )
        db.session.add(bindingdb_bioassay_obj)

# ==============================================================
# Update db BindingDB datestamp
# ==============================================================

current_crawl_datestamp_row = models.DateStamps.query.filter_by(crawl_number=current_crawl_number).first()
current_crawl_datestamp_row.bindingdb_datestamp = now
db.session.commit()
print 'Done.'
