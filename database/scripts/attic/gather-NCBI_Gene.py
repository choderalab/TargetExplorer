# Add data from NCBI Gene to database.
#
# Daniel L. Parton <partond@mskcc.org> - 16 Oct 2013
#

# =================
# IMPORTS
# =================

import sys, os
from lxml import etree
import choderalab as clab

# =================
# PARAMETERS
# =================

kinDB_path = os.path.join('database', 'kinDB-pdb-scop-mut.xml')

parser = etree.XMLParser(remove_blank_text=True)
kinDB = etree.parse(kinDB_path, parser).getroot()

nkinases = len(kinDB)

ncbi_gene_data_dir = os.path.join('external-data', 'NCBI_Gene')
if not os.path.exists(ncbi_gene_data_dir):
    os.mkdir(ncbi_gene_data_dir)

gene2pubmed_filepath = os.path.join(ncbi_gene_data_dir, 'gene2pubmed.gz')

okinDB_path = os.path.join('database', 'kinDB-pdb-scop-mut-ncbi.xml')

debug = False

download_new_gene2pubmed = False

# =================
# MAIN
# =================

if download_new_gene2pubmed:
    clab.NCBI_Gene.retrieve_gene2pubmed(gene2pubmed_filepath)

GeneIDs = [x.get('ID') for x in kinDB.findall('kinase/NCBI_Gene/entry')]
#GeneIDs = ['5562']

tax_ids = ['9606'] * len(GeneIDs)

print 'Extracting publication data from %s...' % gene2pubmed_filepath
publications_xml = clab.NCBI_Gene.get_publications(gene2pubmed_filepath, GeneIDs, tax_ids)
print 'Done extracting publication data.'

for xml_gene_node in publications_xml:
    GeneID = xml_gene_node.get('GeneID')
    kinase_NCBI_Gene_entry_node = kinDB.find('kinase/NCBI_Gene/entry[@ID="%(GeneID)s"]' % vars())
    kinase_NCBI_Gene_entry_node.append(xml_gene_node.find('publications'))
    

# =================
# OUTPUT
# =================

okinDB = open(okinDB_path, 'w')
okinDB.write( etree.tostring(kinDB, pretty_print=True) )
okinDB.close()

