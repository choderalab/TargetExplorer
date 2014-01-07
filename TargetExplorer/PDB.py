# Methods for retrieving data from the PDB
#
# Daniel L. Parton <partond@mskcc.org> - 16 Feb 2013

#==============================================================================
# IMPORTS
#==============================================================================

import urllib2

#==============================================================================
# METHODS
#==============================================================================

def retrieve_sifts(pdb_id):
    '''Retrieves a SIFTS .xml file, given a PDB ID. Works by modifying the PDBe download URL.
    File is retrieved compressed and returned uncompressed.
    Also removes annoying namespace stuff.
    '''
    import re, gzip, StringIO
    #sifts_download_base_url='http://www.ebi.ac.uk/pdbe-srv/view/files/sifts/'
    sifts_download_base_url='ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/'
    url = sifts_download_base_url + pdb_id.lower() + '.xml.gz'
    print url
    response = urllib2.urlopen(url)
    sifts_page = response.read(100000000) # Max 100MB
    # Decompress string
    sifts_page = gzip.GzipFile(fileobj=StringIO.StringIO(sifts_page)).read()

    # Remove all attribs from the entry tag, and the entire rdf tag and contents
    sifts_page_processed = ''
    skip_rdf_tag_flag = False
    for line in sifts_page.splitlines():
        if line[0:6] == '<entry':
            sifts_page_processed += '<entry>' + '\n'
        elif line[0:7] == '  <rdf:':
            skip_rdf_tag_flag = True
            pass
        elif line[0:8] == '  </rdf:':
            skip_rdf_tag_flag = False
            pass
        else:
            if skip_rdf_tag_flag:
                continue
            sifts_page_processed += line + '\n'
    return sifts_page_processed

def retrieve_pdb(pdb_id,compressed='no'):
    '''Retrieves a PDB file, given a PDB ID. Works by modifying the PDB download URL.
    '''
    pdb_download_base_url='http://www.rcsb.org/pdb/files/'
    url = pdb_download_base_url + pdb_id + '.pdb'
    if compressed == 'yes':
        url += '.gz'
    response = urllib2.urlopen(url)
    pdb_file = response.read(10000000) # Max 10MB
    return pdb_file

