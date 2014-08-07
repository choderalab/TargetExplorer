#!/usr/bin/env python
#
# Add data from NCBI Gene to database.
#
# New gene2pubmed.gz file downloaded from NCBI Gene only if existing one is > 7 days old, or if --forcedl flag is used.
#
# Daniel L. Parton <partond@mskcc.org> - 16 Oct 2013
#

# =================
# Imports
# =================

import os
import datetime
import argparse
import TargetExplorer
from flaskapp import models, db
import pandas as pd

# =================
# Parameters
# =================

external_data_dir = 'external-data'
ncbi_gene_data_dir = os.path.join(external_data_dir, 'NCBI_Gene')

if not os.path.exists(ncbi_gene_data_dir):
    os.mkdir(ncbi_gene_data_dir)

gene2pubmed_filepath = os.path.join(ncbi_gene_data_dir, 'gene2pubmed.gz')

argparser = argparse.ArgumentParser(description='Gather UniProt')
argparser.add_argument('--use_existing_gene2pubmed', help='Do not download a new gene2pubmed.gz file. Only works if an existing file is present.', action='store_true', default=False)
args = argparser.parse_args()

now = datetime.datetime.utcnow()

# get current crawl number
crawldata_row = models.CrawlData.query.first()
current_crawl_number = crawldata_row.current_crawl_number
print 'Current crawl number: %d' % current_crawl_number

# ==============================================================
# Download new gene2pubmed.gz file from NCBI Gene
# ==============================================================

# Unless args.use_existing_gene2pubmed is set to true, retrieve a new file from NCBI Gene
if os.path.exists(gene2pubmed_filepath) and args.use_existing_gene2pubmed:
    print 'gene2pubmed file found at:', gene2pubmed_filepath
else:
    print 'Retrieving new Gene2PubMed file from NCBI server...'
    TargetExplorer.NCBI_Gene.retrieve_gene2pubmed(gene2pubmed_filepath)

print ''

# ==============================================================
# Extract publication data from Gene2PubMed file
# ==============================================================

print 'Extracting publication data from gene2pubmed file...'

# first get list of Gene IDs from the db
db_gene_ids = [value_tuple[0] for value_tuple in models.NCBIGeneEntry.query.filter_by(crawl_number=current_crawl_number).values(models.NCBIGeneEntry.gene_id)]

# now load the gene2pubmed file using pandas
# skiprows=1 -- skips first row, which is a non-tab-separated column header. Assign column names manually
# usecols=[1, 2] -- skips first column (taxonomy ID)
# the returned object is a pandas Series (squeeze=True), indexed on gene_id
pmids_by_gene_id = pd.read_table(gene2pubmed_filepath, compression='gzip', skiprows=1, names=['gene_id', 'pmid'], usecols=[1, 2], index_col=0, squeeze=True)

# iterate through the db Gene IDs
i = 0
# dbentry_npubs: npubs for each DBEntry, hashed by DBEntry.id
dbentry_npubs = {}
for db_gene_id in db_gene_ids:
    if i % 100 == 0:
        print '%d/%d Gene IDs searched' % (i, len(db_gene_ids))

    # check if the db Gene ID is present in the gene2pubmed data
    if db_gene_id in pmids_by_gene_id:
        # if so, get the PMIDs corresponding to that Gene ID
        matching_pmids = pmids_by_gene_id[pmids_by_gene_id.index == db_gene_id]
        # and get the matching NCBIGeneEntry object from the db
        db_ncbi_gene_entry = models.NCBIGeneEntry.query.filter_by(crawl_number=current_crawl_number, gene_id=db_gene_id).first()

        # initialize the npubs counter for the corresponding DBEntry, if necessary
        if db_ncbi_gene_entry.dbentry_id not in dbentry_npubs:
            dbentry_npubs[db_ncbi_gene_entry.dbentry_id] = 0

        # matching_pmids may be a single PMID (type:int) or a series of PMIDs (type:pd.Series)
        if type(matching_pmids) == pd.Series:
            dbentry_npubs[db_ncbi_gene_entry.dbentry_id] += len(matching_pmids)
            for pmid in matching_pmids:
                # create NCBIGenePublication object and add to db
                ncbi_gene_publication_obj = models.NCBIGenePublication(crawl_number=current_crawl_number, pmid=pmid, ncbi_gene_entry=db_ncbi_gene_entry)
                db.session.add(ncbi_gene_publication_obj)
        else:
            dbentry_npubs[db_ncbi_gene_entry.dbentry_id] += 1
            # create NCBIGenePublication object and add to db
            ncbi_gene_publication_obj = models.NCBIGenePublication(crawl_number=current_crawl_number, pmid=matching_pmids, ncbi_gene_entry=db_ncbi_gene_entry)
            db.session.add(ncbi_gene_publication_obj)

    i += 1

for dbentry in models.DBEntry.query.filter_by(crawl_number=current_crawl_number).all():
    if dbentry.id in dbentry_npubs:
        dbentry.npubs = dbentry_npubs[dbentry.id]
    else:
        dbentry.npubs = 0

# ==============================================================
# Update db NCBI Gene datestamp
# ==============================================================

current_crawl_datestamp_row = models.DateStamps.query.filter_by(crawl_number=current_crawl_number).first()
current_crawl_datestamp_row.ncbi_gene_datestamp = now
db.session.commit()
print 'Done.'
