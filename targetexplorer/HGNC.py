# Routines for retrieving HGNC nomenclature info
# Data can only be downloaded from HGNC as a flat text file
#
# Daniel L. Parton <partond@mskcc.org> - 11 Oct 2013

import urllib2, re, gzip, os, sys
from lxml import etree

def retrieve_data_file(HGNC_data_filepath):
    url = 'ftp://ftp.ebi.ac.uk/pub/databases/genenames/locus_groups/protein-coding_gene.txt.gz'
    response = urllib2.urlopen(url)
    page = response.read(1000000000)
    # Save compressed file
    with open(HGNC_data_filepath, 'w') as compressed_file:
        compressed_file.write(page)

class HGNC_data:
    def __init__(self, HGNC_data_filepath):
        # Open compressed file and decompress the contained data
        with gzip.open(HGNC_data_filepath, 'r') as decompressed_file:
            lines = decompressed_file.readlines()

        self.IDkeyed_data = {}

        for line in lines[1:]:
            words = line.split('\t')
            ID = words[0]
            self.IDkeyed_data[ID] = words



