# Add mutation data from cBioPortal to database.
#
# Daniel L. Parton <partond@mskcc.org> - 29 Aug 2013; 16 Oct 2013
#
# TODO: Implement a scheme for assigning IDs to mutants.
# TODO: Map mutation sequence positions to UniProt canonical sequence.

# =================
# IMPORTS
# =================

import sys, os
from lxml import etree
import choderalab as clab

# =================
# PARAMETERS
# =================

kinDB_path = os.path.join('database', 'kinDB-pdb-scop.xml')

parser = etree.XMLParser(remove_blank_text=True)
kinDB = etree.parse(kinDB_path, parser).getroot()

nkinases = len(kinDB)

okinDB_path = os.path.join('database', 'kinDB-pdb-scop-mut.xml')
cbioportal_data_dir = os.path.join('external-data', 'cBioPortal')
if not os.path.exists(cbioportal_data_dir):
    os.mkdir(cbioportal_data_dir)
output_cbioportal_xml_filepath = os.path.join(cbioportal_data_dir, 'cbioportal-mutations.xml')

debug = False

download_cbioportal_data = False

# =================
# RETRIEVE MUTANTS FROM CBIOPORTAL
# XXX IMPORTANT: sequence positions are not yet mapped to UniProt canonical sequence
# =================

# Get list of HGNC gene symbols from database
gene_list = [ x.get('Approved_Symbol') for x in kinDB.findall('kinase/HGNC/entry') ]

# Get list of all cancer studies available in cBioPortal
cancer_studies = clab.cBioPortal.get_cancer_studies()
#cancer_studies = ['kich_tcga']

if download_cbioportal_data == True:
    clab.cBioPortal.retrieve_mutants_xml(output_cbioportal_xml_filepath, cancer_studies, gene_list, debug=debug)

cbioportal_xml = etree.parse(output_cbioportal_xml_filepath, parser).getroot()

for gene_symbol in gene_list:
    kinDB_gene = kinDB.find('kinase/HGNC/entry[@Approved_Symbol="%(gene_symbol)s"]/../..' % vars())
    cbioportal_gene = cbioportal_xml.find('gene[@gene_symbol="%(gene_symbol)s"]' % vars())
    mutants_node = etree.SubElement(kinDB_gene, 'mutants')
    for mutant in cbioportal_gene:
        mutants_node.append(mutant)

# =================
# ADD THESE MUTANTS MANUALLY
# XXX TO BE DEPRECATED
# =================

#src = kinDB.find('kinase/uniprot[@AC="P12931"]/..')
#src_mutants = etree.SubElement(src, 'mutants')
#src_gatekeeper_mutant = etree.SubElement(src_mutants, 'mutant')
#src_gatekeeper_mutant.set('comment', 'Gatekeeper mutant')
#src_gatekeeper_mutations = etree.SubElement(src_gatekeeper_mutant, 'mutation')
#src_gatekeeper_mutations.text = 'T341I'
#
#mutant_id = 0
#WT_target_id = src.find('uniprot/pk_domain').get('kinDB_id')
#mutant_target_id = WT_target_id + '_M%d' % mutant_id
#src_gatekeeper_mutant.set('id', str(mutant_id))
#src_gatekeeper_mutant.set('target_id', mutant_target_id)
#src_gatekeeper_mutant.set('pk_domain_id', '0')
#
#abl1 = kinDB.find('kinase/uniprot[@AC="P00519"]/..')
#abl1_mutants = etree.SubElement(abl1, 'mutants')
#abl1_gatekeeper_mutant = etree.SubElement(abl1_mutants, 'mutant')
#abl1_gatekeeper_mutant.set('comment', 'Gatekeeper mutant')
#abl1_gatekeeper_mutations = etree.SubElement(abl1_gatekeeper_mutant, 'mutation')
#abl1_gatekeeper_mutations.text = 'T315I'
#
#mutant_id = 0
#WT_target_id = abl1.find('uniprot/pk_domain').get('kinDB_id')
#mutant_target_id = WT_target_id + '_M%d' % mutant_id
#abl1_gatekeeper_mutant.set('id', str(mutant_id))
#abl1_gatekeeper_mutant.set('target_id', mutant_target_id)
#abl1_gatekeeper_mutant.set('pk_domain_id', '0')

# =================
# OUTPUT
# =================

okinDB = open(okinDB_path, 'w')
okinDB.write( etree.tostring(kinDB, pretty_print=True) )
okinDB.close()

