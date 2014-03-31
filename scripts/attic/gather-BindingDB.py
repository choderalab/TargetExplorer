# Add bioassay data from BindingDB to database.
#
# Daniel L. Parton <partond@mskcc.org> - 25 Oct 2013
#

# =================
# Imports
# =================

import sys, os, gzip, subprocess, re
from lxml import etree
import choderalab as clab
from multiprocessing import Pool

# =================
# Parameters
# =================

kinDB_path = os.path.join('database', 'kinDB-pdb-scop-mut-ncbi.xml')

parser = etree.XMLParser(remove_blank_text=True)
kinDB = etree.parse(kinDB_path, parser).getroot()

nkinases = len(kinDB)

okinDB_filepath = os.path.join('database', 'kinDB-pdb-scop-mut-ncbi-bindingdb.xml')
bindingdb_data_dir = os.path.join('external-data', 'BindingDB')
if not os.path.exists(bindingdb_data_dir):
    os.mkdir(bindingdb_data_dir)
bindingdb_all_data_filepath = os.path.join(bindingdb_data_dir, 'BindingDB_All.tab')

try:
    grep_path = sys.argv[sys.argv.index('-grep_path') + 1]
except ValueError:
    grep_path = 'grep'

debug = False

update_bindingdb_data = False

# =================
# Download all BindingDB data and decompress if necessary
# =================

if update_bindingdb_data:
    clab.BindingDB.retrieve_all_BindingDB_data(bindingdb_all_data_filepath + '.gz')

# decompress
if not os.path.exists(bindingdb_all_data_filepath) or update_bindingdb_data:
    subprocess.call('gunzip -c %(bindingdb_all_data_filepath)s.gz > %(bindingdb_all_data_filepath)s' % vars(), shell=True)

# =================
# Use grep to extract information about each kinase in the database
# =================

# Overview:
#   * first use grep -E to extract all lines matching any kinase uniprot AC
#   * using multiprocessing for each kinase:
#       * grep through to create data file for each kinase
#       * parse each kinase data file and add to kinDB
#       * delete individual kinase data files

kinase_uniprotACs = [uniprot_node.get('AC') for uniprot_node in kinDB.findall('kinase/uniprot')]

#kinase_uniprotACs = ['P00519']

regex_search = '|'.join(kinase_uniprotACs)
bindingdb_all_kinases_data_filepath = os.path.join(bindingdb_data_dir, 'tmp-bindingdb-all_kinases.tab')

if not os.path.exists(bindingdb_all_kinases_data_filepath) or update_bindingdb_data:
    with open(bindingdb_all_data_filepath, 'r') as bindingdb_all_data_file:
        subprocess.call('%(grep_path)s -E "%(regex_search)s" %(bindingdb_all_data_filepath)s > %(bindingdb_all_kinases_data_filepath)s' % vars(), shell=True)

def extract_bindingdb(input_data):
    AC, grep_path, bindingdb_all_kinases_data_filepath = input_data
    print AC
    kinase_bindingdb_data_filepath = os.path.join(bindingdb_data_dir, AC + '.tab')
    subprocess.call('%(grep_path)s "%(AC)s" %(bindingdb_all_kinases_data_filepath)s > %(kinase_bindingdb_data_filepath)s' % vars(), shell=True)
    bioassays_data = []
    with open(kinase_bindingdb_data_filepath, 'r') as bindingdb_file:
        for line in bindingdb_file:
            returnedACs = clab.BindingDB.get_ACs(line)
            for returnedAC in returnedACs:
                if returnedAC == AC:
                    bioassay_data = clab.BindingDB.get_bioassay_data(line) # returns a dict
                    bioassays_data.append( bioassay_data )
    os.remove(kinase_bindingdb_data_filepath)
    return (AC, bioassays_data)

if __name__ == '__main__':
    pool = Pool()
    input_data = [(AC, grep_path, bindingdb_all_kinases_data_filepath) for AC in kinase_uniprotACs]
    results = pool.map(extract_bindingdb, input_data)

for result in results:
    AC, bioassays_data = result
    kinase = kinDB.find('kinase/uniprot[@AC="%s"]/..' % AC)
    # Create bioassays node
    bioassays_node = etree.SubElement(kinase, 'bioassays')
    for bioassay_data in bioassays_data:
        bioassay_node = etree.SubElement(bioassays_node, 'bioassay')
        for data_key in bioassay_data.keys():
            bioassay_node.set(data_key, bioassay_data[data_key])

with open(okinDB_filepath, 'w') as okinDB_file:
    okinDB_file.write( etree.tostring(kinDB, pretty_print=True) )

