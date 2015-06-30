# Functions for retrieval of bioassay data from BindingDB
#
# Daniel L. Parton <partond@mskcc.org> - 25 Oct 2013

import urllib2, sys, os, subprocess
from lxml import etree

def retrieve_all_BindingDB_data(bindingdb_all_data_filepath, decompress=True):
    '''
    Retrieves all BindingDB data in compressed tab-separated format, then decompresses and writes to file.
    Uses gunzip for decompression.
    '''
    print 'Downloading BindingDB_All.tab.gz...'
    url = 'http://bindingdb.org/bind/downloads/BindingDB_All_2015m5.tsv.zip'
    response = urllib2.urlopen(url)
    page = response.read(1000000000)
    with open(bindingdb_all_data_filepath+'.gz', 'wb') as ofile:
        ofile.write(page)
    print 'Finished downloading BindingDB_All.tab.gz'

    # decompress
    if decompress:
        print 'Decompressing BindingDB_All.tab.gz...'
        subprocess.call('gunzip -f %s' % bindingdb_all_data_filepath+'.gz', shell=True)
        print 'Finished decompressing BindingDB_All.tab.gz'

def get_ACs(line):
    words = line.split('\t')
    UniProtACs = words[20].split(' ')
    return UniProtACs

def get_bioassay_data(line):
    words = line.split('\t')
    SMILES_string = words[0]
    BindingDB_monomerID = words[1]
    ChEMBL_ID = words[8]
    zinc_id = words[13]
    data_origin = words[15]
    target_biomolecule = words[16]
    #target_source_organism = words[17]
    #target_sequence = words[18]
    #pdbID = words[19]
    #UniProtAC = words[20]
    Ki = words[22]
    IC50 = words[23]
    Kd = words[24]
    EC50 = words[25]
    kon = words[26]
    koff = words[27]
    PMID = words[28]
    DOI = words[29]
    pH = words[32]
    temperature = words[33]

    bioassay_data = {}

    bioassay_data['ligand_SMILES_string'] = SMILES_string
    bioassay_data['ligand_BindingDB_ID'] = BindingDB_monomerID
    bioassay_data['ligand_ChEMBL_ID'] = ChEMBL_ID
    bioassay_data['ligand_zinc_id'] = zinc_id
    bioassay_data['BindingDB_source'] = data_origin
    bioassay_data['target_name'] = target_biomolecule
    #bioassay_data['target_sequence'] = target_sequence
    #bioassay_data['target_source_organism'] = target_source_organism

    bioassay_data['Ki'] = Ki # XXX remove
    bioassay_data['IC50'] = IC50
    bioassay_data['Kd'] = Kd
    bioassay_data['EC50'] = EC50
    bioassay_data['kon'] = kon
    bioassay_data['koff'] = koff
    bioassay_data['pH'] = pH
    bioassay_data['temperature'] = temperature
    # if Ki != 'n/a':
    #     bioassay_data['Ki'] = Ki
    # if IC50 != 'n/a':
    #     bioassay_data['IC50'] = IC50
    # if Kd != 'n/a':
    #     bioassay_data['Kd'] = Kd
    # if EC50 != 'n/a':
    #     bioassay_data['EC50'] = EC50
    # if kon != 'n/a':
    #     bioassay_data['kon'] = kon
    # if koff != 'n/a':
    #     bioassay_data['koff'] = koff
    # if pH != 'n/a':
    #     bioassay_data['pH'] = pH
    # if temperature != 'n/a':
    #     bioassay_data['temperature'] = temperature

    bioassay_data['PMID'] = PMID
    bioassay_data['DOI'] = DOI

    return bioassay_data

