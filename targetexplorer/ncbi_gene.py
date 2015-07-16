import urllib2
import os
import datetime
from targetexplorer.flaskapp import models, db
import pandas as pd


class GatherNCBIGene(object):
    def __init__(self,
                 use_existing_gene2pubmed=False,
                 run_main=True
                 ):
        if run_main:
            external_data_dir = 'external-data'
            ncbi_gene_data_dir = os.path.join(external_data_dir, 'NCBI_Gene')

            if not os.path.exists(ncbi_gene_data_dir):
                os.mkdir(ncbi_gene_data_dir)

            gene2pubmed_filepath = os.path.join(ncbi_gene_data_dir, 'gene2pubmed.gz')

            now = datetime.datetime.utcnow()

            # get current crawl number
            crawldata_row = models.CrawlData.query.first()
            current_crawl_number = crawldata_row.current_crawl_number
            print 'Current crawl number: %d' % current_crawl_number

            # ==============================================================
            # Download new gene2pubmed.gz file from NCBI Gene
            # ==============================================================

            if os.path.exists(gene2pubmed_filepath) and use_existing_gene2pubmed:
                print 'gene2pubmed file found at:', gene2pubmed_filepath
            else:
                print 'Retrieving new Gene2PubMed file from NCBI server...'
                retrieve_gene2pubmed(gene2pubmed_filepath)

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

                if i % 100 == 0 or i == len(db_gene_ids):
                    print '%d/%d Gene IDs searched' % (i, len(db_gene_ids))

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



def retrieve_gene2pubmed(gene2pubmed_gzfilepath):
    '''
    Download compressed gene2pubmed from NCBI FTP site and write to file.
    '''
    url = 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz'
    response = urllib2.urlopen(url)
    page = response.read(100000000)
    with open(gene2pubmed_gzfilepath, 'w') as gene2pubmed_gzfile:
        gene2pubmed_gzfile.write(page)
    # # Update metadata
    # metadata_root = etree.parse(targetexplorer.DB.external_data_metadata_filepath, parser).getroot()
    # gene2pubmed_node = metadata_root.find('NCBI_Gene/gene2pubmed')
    # if gene2pubmed_node is None:
    #     NCBI_Gene_node = etree.SubElement(metadata_root, 'NCBI_Gene')
    #     gene2pubmed_node = etree.SubElement(NCBI_Gene_node, 'gene2pubmed')
    # gene2pubmed_node.set('filepath', gene2pubmed_gzfilepath)
    # now = datetime.datetime.utcnow()
    # datestamp = now.strftime(targetexplorer.DB.datestamp_format_string)
    # gene2pubmed_node.set('datestamp', datestamp)
    # with open(targetexplorer.DB.external_data_metadata_filepath, 'w') as xml_out_file:
    #     xml_out_file.write(etree.tostring(metadata_root, pretty_print=True))


def retrieve_xml(GeneID):
    '''
    Don't use this to retrieve publication data.
    '''
    url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id=%(GeneID)s&retmode=xml' % vars()
    response = urllib2.urlopen(url)
    page = response.read(100000000000)
    lines = page.splitlines()

    ofilepath = 'tmp-gene.xml'
    with open(ofilepath, 'w') as ofile:
        ofile.write(page)
