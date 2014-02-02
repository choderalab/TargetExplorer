# Add bioassay data from BindingDB to database.
#
# Daniel L. Parton <partond@mskcc.org> - 25 Oct 2013
#

# =================
# Imports
# =================

import sys, os, subprocess, datetime
from lxml import etree
import TargetExplorer as clab
from multiprocessing import Pool

# =================
# Parameters
# =================

DB_script_ID = 'BindingDB'

if '-stage' in sys.argv:
    run_mode = 'stage'
elif '-dev' in sys.argv:
    run_mode = 'dev'
else:
    run_mode = 'nowrite'

print 'Running in mode: %s' % run_mode

database_dir = 'database'
external_data_dir = os.path.join('external-data')
bindingdb_data_dir = os.path.join(external_data_dir, 'BindingDB')
if not os.path.exists(bindingdb_data_dir):
    os.mkdir(bindingdb_data_dir)
bindingdb_all_data_filepath = os.path.join(bindingdb_data_dir, 'BindingDB_All.tab')

DBstage_filepath = os.path.join(database_dir, 'database-stage.xml')
if not os.path.exists(DBstage_filepath):
    raise Exception, '%s not found.' % DBstage_filepath

if run_mode != 'nowrite':
    DB_out_filename = 'database-%(run_mode)s.xml' % vars()
    DB_out_filepath = os.path.join(database_dir, DB_out_filename)

if '--forcedl' in sys.argv:
    force_bindingdb_all_download = True
else:
    force_bindingdb_all_download = False

now = datetime.datetime.utcnow()
datestamp = now.strftime(clab.DB.datestamp_format_string)

parser = etree.XMLParser(remove_blank_text=True)

try:
    grep_path = sys.argv[sys.argv.index('-grep_path') + 1]
except ValueError:
    grep_path = 'grep'

debug = False

# =================
# Download new BindingDB data if necessary
# =================

clab.DB.retrieve_external_data(DB_script_ID, forcedl=force_bindingdb_all_download)

# =================
# Read in the existing DB
# =================

print 'Reading', DBstage_filepath
DB_root = etree.parse(DBstage_filepath, parser).getroot()
nentries = len(DB_root)
print 'Number of entries:', nentries

# =================
# Use grep to extract information for each entry in the DB
# =================

# Overview:
#   * first use grep -E to extract all lines matching any kinase uniprot AC
#   * using multiprocessing for each kinase:
#       * grep through to create data file for each kinase
#       * parse each kinase data file and add to kinDB
#       * delete individual kinase data files

DB_uniprotACs = [uniprot_node.get('AC') for uniprot_node in DB_root.findall('entry/UniProt')]

#DB_uniprotACs = ['P00519']

regex_search = '|'.join(DB_uniprotACs)
bindingdb_matches_filepath = os.path.join(bindingdb_data_dir, 'tmp-bindingdb-matches.tab')

#if not os.path.exists(bindingdb_matches_filepath): # XXX DEBUG
if True:
    print 'First pass through BindingDB data with grep... - extracting all lines which match UniProt ACs in the DB (performance can be variable - testing has indicated approx. 1 hr on a 2012 MBP, but only 2 mins or so on a Linux machine; probably due to different GNU grep versions)'
    with open(bindingdb_all_data_filepath, 'r') as bindingdb_all_data_file:
        subprocess.call('%(grep_path)s -E "%(regex_search)s" %(bindingdb_all_data_filepath)s > %(bindingdb_matches_filepath)s' % vars(), shell=True)

def extract_bindingdb(input_data):
    AC, grep_path, bindingdb_matches_filepath = input_data
    print AC
    DB_entry_bindingdb_data_filepath = os.path.join(bindingdb_data_dir, AC + '.tab')
    subprocess.call('%(grep_path)s "%(AC)s" %(bindingdb_matches_filepath)s > %(DB_entry_bindingdb_data_filepath)s' % vars(), shell=True)
    bioassays_data = []
    with open(DB_entry_bindingdb_data_filepath, 'r') as bindingdb_file:
        for line in bindingdb_file:
            returnedACs = clab.BindingDB.get_ACs(line)
            for returnedAC in returnedACs:
                if returnedAC == AC:
                    bioassay_data = clab.BindingDB.get_bioassay_data(line) # returns a dict
                    bioassays_data.append( bioassay_data )
    os.remove(DB_entry_bindingdb_data_filepath)
    return (AC, bioassays_data)

if __name__ == '__main__':
    pool = Pool()
    input_data = [(AC, grep_path, bindingdb_matches_filepath) for AC in DB_uniprotACs]
    results = pool.map(extract_bindingdb, input_data)

for result in results:
    AC, bioassays_data = result
    entry = DB_root.find('entry/UniProt[@AC="%s"]/..' % AC)
    # Create bioassays node
    bioassays_node = etree.SubElement(entry, 'bioassays')
    for bioassay_data in bioassays_data:
        bioassay_node = etree.SubElement(bioassays_node, 'bioassay')
        for data_key in bioassay_data.keys():
            bioassay_node.set(data_key, bioassay_data[data_key])

#==============================================================================
# If staging, update date_run
#==============================================================================

if run_mode == 'stage':
    DB_root.set('gather_bindingdb_last_run', datestamp)

# ==============================================================
# If staging, compare new and old DBs
# ==============================================================

data_modified = False

if run_mode == 'stage':
    # Parse the old DB
    DBold_root = etree.parse(DBstage_filepath, parser).getroot()

    # First a quick check to see whether the numbers of bioassay nodes match
    DB_bioassay_nodes = DB_root.findall('entry/bioassays/bioassay[@source="BindingDB"]')
    DBold_bioassay_nodes = DBold_root.findall('entry/bioassays/bioassay[@source="BindingDB"]')
    if len(DBold_bioassay_nodes) != len(DB_bioassay_nodes):
        print 'Comparison of latest BindingDB data with data in %s indicates changes. File will be overwritten with new data. Number of bioassays differs by: %d' % (DB_out_filepath, len(DB_bioassay_nodes) - len(DBold_bioassay_nodes))
        data_modified = True

    # If the numbers of bioassay nodes match, then proceed to a full comparison using GNU diff
    else:
        DB_comparison_string = etree.tostring(DB_root, pretty_print=True)
        DB_comparison_string = '\n'.join( DB_comparison_string.splitlines()[1:] ) # Ignore the first line, which contains datestamps (don't want to include these in the comparison)
        DBold_comparison_string = etree.tostring(DBold_root, pretty_print=True)
        DBold_comparison_string = '\n'.join( DBold_comparison_string.splitlines()[1:] ) # Ignore the first line, which contains datestamps (don't want to include these in the comparison)

        diff_output = clab.DB.diff_DB_comparison_strings(DBold_comparison_string, DB_comparison_string)

        if len(diff_output) > 0:
            print 'Comparison of latest BindingDB data with data in %s indicates changes. File will be overwritten. Lines in diff comparison: %s' % (DBstage_filepath, len(diff_output))
            data_modified = True
        else:
            print 'Comparison of latest BindingDB data with data in %s indicates no changes. File will be rewritten with updated gather_bindingdb_late_run attrib, but other data will not be modified.' % DB_out_filepath
            data_modified = False

#==============================================================================
# If staging and there have been modifications, update date_modified
#==============================================================================

if run_mode == 'stage' and data_modified:
    DB_root.set('gather_bindingdb_last_modif', datestamp)

# =================
# Output
# =================

if run_mode == 'stage' and data_modified:
    clab.DB.writeDB(DB_root, DB_out_filepath)

if run_mode == 'stage' and not data_modified:
    DBold_root = etree.parse(DB_out_filepath, parser).getroot()
    DBold_root.set('gather_bindingdb_last_run', datestamp)
    clab.DB.writeDB(DBold_root, DB_out_filepath)

elif run_mode == 'dev':
    clab.DB.writeDB(DB_root, DB_out_filepath)



