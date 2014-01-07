# Add data from NCBI Gene to database.
#
# New gene2pubmed.gz file downloaded from NCBI Gene only if existing one is > 7 days old, or if --forcedl flag is used.
#
# Daniel L. Parton <partond@mskcc.org> - 16 Oct 2013
#

# =================
# Imports
# =================

import sys, os, datetime
from lxml import etree
import choderalab as clab

# =================
# Parameters
# =================

if '-stage' in sys.argv:
    run_mode = 'stage'
elif '-dev' in sys.argv:
    run_mode = 'dev'
else:
    run_mode = 'nowrite'

print 'Running in mode: %s' % run_mode

database_dir = 'database'
external_data_dir = os.path.join(clab.kinome_rootdir, 'database', 'external-data')
ncbi_gene_data_dir = os.path.join(external_data_dir, 'NCBI_Gene')
gene2pubmed_filepath = os.path.join(ncbi_gene_data_dir, 'gene2pubmed.gz')

if not os.path.exists(ncbi_gene_data_dir):
    os.mkdir(ncbi_gene_data_dir)

DBstage_filepath = os.path.join(database_dir, 'database-stage.xml')
if not os.path.exists(DBstage_filepath):
    raise Exception, '%s not found.' % DBstage_filepath

if run_mode != 'nowrite':
    DB_out_filename = 'database-%(run_mode)s.xml' % vars()
    DB_out_filepath = os.path.join(database_dir, DB_out_filename)

if '--forcedl' in sys.argv:
    force_gene2pubmed_download = True
else:
    force_gene2pubmed_download = False

now = datetime.datetime.utcnow()
datestamp = now.strftime(clab.DB.datestamp_format_string)

parser = etree.XMLParser(remove_blank_text=True)

days_elapsed_for_force_download = 7

# ==============================================================
# Read in the existing db
# ==============================================================

print 'Reading', DBstage_filepath
DB_root = etree.parse(DBstage_filepath, parser).getroot()
nentries = len(DB_root)
print 'Number of entries:', nentries

# ==============================================================
# Download new gene2pubmed.gz file from NCBI Gene if necessary
# ==============================================================

# First check if local copy of gene2pubmed.gz already exists
if os.path.exists(gene2pubmed_filepath):
    print 'Gene2PubMed file found at:', gene2pubmed_filepath
# If not, download it
else:
    print 'Gene2PubMed file not found.'
    print 'Retrieving new Gene2PubMed file from NCBI server...'
    clab.NCBI_Gene.retrieve_gene2pubmed(gene2pubmed_filepath)

# Check when the Gene2PubMed file was retrieved and download a new one if time elapsed is more than days_elapsed_for_force_download
external_data_metadata_root = etree.parse(clab.DB.external_data_metadata_filepath, parser).getroot()
gene2pubmed_datestamp = external_data_metadata_root.find('NCBI_Gene/gene2pubmed').get('datestamp')
gene2pubmed_datestamp = datetime.datetime.strptime(gene2pubmed_datestamp, clab.DB.datestamp_format_string)
time_elapsed = now - gene2pubmed_datestamp
if (time_elapsed.days > days_elapsed_for_force_download) or (force_gene2pubmed_download == True):
    if time_elapsed.days > days_elapsed_for_force_download:
        print 'Gene2PubMed file is %d days old.' % time_elapsed.days
    if force_gene2pubmed_download:
        print 'Forcing retrieval of new Gene2PubMed file.'
        download_new_gene2pubmed = True
    else:
        while True:
            user_response = raw_input('Suggest retrieving new Gene2PubMed file from NCBI server. Proceed? [y] ')
            if user_response in ['y', '']:
                download_new_gene2pubmed = True
                break
            elif user_response == 'n':
                download_new_gene2pubmed = False
                break
            else:
                print 'User input not understood. Please try again.'

    if download_new_gene2pubmed:
        print 'Retrieving new Gene2PubMed file from NCBI server...'
        clab.NCBI_Gene.retrieve_gene2pubmed(gene2pubmed_filepath)

print ''

# ==============================================================
# Extract publication data from Gene2PubMed file
# ==============================================================

GeneIDs = [x.get('ID') for x in DB_root.findall('entry/NCBI_Gene/entry')]
#GeneIDs = ['5562']

#tax_ids = ['9606'] * len(GeneIDs)

print 'Extracting publication data from %s...' % gene2pubmed_filepath
# uses multiprocessing
PMIDs_by_GeneID = clab.NCBI_Gene.get_publications(gene2pubmed_filepath, GeneIDs)
print 'Done extracting publication data.'

#for xml_gene_node in publications_xml:
#    GeneID = xml_gene_node.get('GeneID')
#    DBentry_NCBI_Gene_entry_node = DB_root.find('entry/NCBI_Gene/entry[@ID="%s"]' % GeneID)
#    DBentry_NCBI_Gene_entry_node.append(xml_gene_node.find('publications'))

for GeneID in PMIDs_by_GeneID.keys():
    DBentry_NCBI_Gene_entry_node = DB_root.find('entry/NCBI_Gene/entry[@ID="%s"]' % GeneID)
    PMIDs = PMIDs_by_GeneID[GeneID]
    if len(PMIDs) > 0:
        publications_node = etree.SubElement(DBentry_NCBI_Gene_entry_node, 'publications')
    for PMID in PMIDs:
        etree.SubElement(publications_node, 'publication').set('PMID', PMID)

#==============================================================================
# If staging, update date_run
#==============================================================================

if run_mode == 'stage':
    DB_root.set('gather_ncbi_gene_last_run', datestamp)

# ==============================================================
# If staging, compare new and old DBs
# ==============================================================

data_modified = False

if run_mode == 'stage':
    # Parse the old DB
    DBold_root = etree.parse(DBstage_filepath, parser).getroot()

    # First a quick check to see whether the numbers of publication nodes match
    DB_publication_nodes = DB_root.findall('entry/NCBI_Gene/entry/publications/publication')
    DBold_publication_nodes = DBold_root.findall('entry/NCBI_Gene/entry/publications/publication')
    if len(DBold_publication_nodes) != len(DB_publication_nodes):
        print 'Comparison of latest NCBI Gene data with data in %s indicates changes. File will be overwritten with new data. Number of publications differs by: %d' % (DB_out_filepath, len(DB_publication_nodes) - len(DBold_publication_nodes))
        data_modified = True

    # If the numbers of publication nodes match, then proceed to a full comparison using GNU diff
    else:
        DB_comparison_string = etree.tostring(DB_root, pretty_print=True)
        DB_comparison_string = '\n'.join( DB_comparison_string.splitlines()[1:] ) # Ignore the first line, which contains datestamps (don't want to include these in the comparison)
        DBold_comparison_string = etree.tostring(DBold_root, pretty_print=True)
        DBold_comparison_string = '\n'.join( DBold_comparison_string.splitlines()[1:] ) # Ignore the first line, which contains datestamps (don't want to include these in the comparison)

        diff_output = clab.DB.diff_DB_comparison_strings(DBold_comparison_string, DB_comparison_string)

        if len(diff_output) > 0:
            print 'Comparison of latest NCBI Gene data with data in %s indicates changes. File will be overwritten. Lines in diff comparison: %s' % (DBstage_filepath, len(diff_output))
            data_modified = True
        else:
            print 'Comparison of latest NCBI Gene data with data in %s indicates no changes. File will be rewritten with updated gather_ncbi_gene_late_run attrib, but other data will not be modified.' % DB_out_filepath
            data_modified = False

#==============================================================================
# If staging and there have been modifications, update date_modified
#==============================================================================

if run_mode == 'stage' and data_modified:
    DB_root.set('gather_ncbi_gene_last_modif', datestamp)

# =================
# Output
# =================

if run_mode == 'stage' and data_modified:
    clab.DB.writeDB(DB_root, DB_out_filepath)

if run_mode == 'stage' and not data_modified:
    DBold_root = etree.parse(DB_out_filepath, parser).getroot()
    DBold_root.set('gather_ncbi_gene_last_run', datestamp)
    clab.DB.writeDB(DBold_root, DB_out_filepath)

elif run_mode == 'dev':
    clab.DB.writeDB(DB_root, DB_out_filepath)

