# Add mutation data from cBioPortal to database.
#
# Daniel L. Parton <partond@mskcc.org> - 29 Aug 2013
#
# TODO: Write external-data metadata
# TODO: Implement a scheme for assigning IDs to mutants.
# TODO: Map mutation sequence positions to UniProt canonical sequence.

# =================
# Imports
# =================

import sys, os, datetime, copy
from lxml import etree
import TargetExplorer as clab

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
external_data_dir = 'external-data'
cbioportal_data_dir = os.path.join(external_data_dir, 'cBioPortal')
output_cbioportal_xml_filepath = os.path.join(cbioportal_data_dir, 'cbioportal-mutations.xml')
if not os.path.exists(cbioportal_data_dir):
    os.mkdir(cbioportal_data_dir)

DBstage_filepath = os.path.join(database_dir, 'database-stage.xml')
if not os.path.exists(DBstage_filepath):
    raise Exception, '%s not found.' % DBstage_filepath

if run_mode != 'nowrite':
    DB_out_filename = 'database-%(run_mode)s.xml' % vars()
    DB_out_filepath = os.path.join(database_dir, DB_out_filename)

verbose = False

download_cbioportal_data = False

now = datetime.datetime.utcnow()
datestamp = now.strftime(clab.DB.datestamp_format_string)

parser = etree.XMLParser(remove_blank_text=True)

# =================
# Read in existing database
# =================

DB_root = etree.parse(DBstage_filepath, parser).getroot()
nentries = len(DB_root)

# =================
# Retrieve mutants from cBioPortal
# XXX IMPORTANT: sequence positions are not yet mapped to UniProt canonical sequence
# =================

# Get list of HGNC gene symbols from database
DB_HGNC_symbols = [ x.get('Approved_Symbol') for x in DB_root.findall('entry/HGNC/entry') ]

# Get list of all cancer studies available in cBioPortal
cancer_studies = clab.cBioPortal.get_cancer_studies()
#cancer_studies = ['kich_tcga']

if download_cbioportal_data == True:
    clab.cBioPortal.retrieve_mutants_xml(output_cbioportal_xml_filepath, cancer_studies, DB_HGNC_symbols, verbose=verbose)

cbioportal_xml = etree.parse(output_cbioportal_xml_filepath, parser).getroot()

# Remove exisiting cBioPortal-derived data from DB_root

existing_cbioportal_mutants = DB_root.findall('entry/mutants/mutant[@source="cBioPortal"]')
for mutant_node in existing_cbioportal_mutants:
    mutants_node = mutant_node.getparent()
    mutants_node.remove(mutant_node)

# Now add back in cBioPortal-derived data from the local cbioportal xml

for HGNC_symbol in DB_HGNC_symbols:
    DB_entry = DB_root.find('entry/HGNC/entry[@Approved_Symbol="%(HGNC_symbol)s"]/../..' % vars())
    cbioportal_gene = cbioportal_xml.find('gene[@gene_symbol="%(HGNC_symbol)s"]' % vars())
    if cbioportal_gene != None:
        mutants_node = DB_entry.find('mutants')
        if mutants_node == None:
            mutants_node = etree.SubElement(DB_entry, 'mutants')
        for mutant in cbioportal_gene:
            mutants_node.append(mutant)

# =================
# Add these mutants manually
# XXX DEPRECATED
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

# =======================
# If staging, update date_run
# =======================

if run_mode == 'stage':
    DB_root.set('gather_cbioportal_last_run', datestamp)

#==============================================================================
# If staging, compare new and old DBs
#==============================================================================

data_modified = False

if run_mode == 'stage':
    # Parse the old DB
    DBold_root = etree.parse(DB_out_filepath, parser).getroot()

    # First a quick check to see whether the total numbers of cBioPortal-derived mutants match
    DB_seq_nodes = DB_root.findall('entry/mutants/mutant[@source="cBioPortal"]')
    DBold_seq_nodes = DBold_root.findall('entry/mutants/mutant[@source="cBioPortal"]')
    if len(DB_seq_nodes) != len(DBold_seq_nodes):
        print 'Comparison of latest cBioPortal-derived data with data in %s indicates changes. DB will be re-written with new data. Number of mutant nodes differs by : %d' % (DB_out_filepath, len(DB_seq_nodes) - len(DBold_seq_nodes))
        data_modified = True

    # If the numbers of PDB chain nodes match, then proceeed to a full comparison of PDB-derived data using diff
    else:
        DB_comparison_root = etree.Element('database')
        DBold_comparison_root = etree.Element('database')

        # XXX Below will need to be modified if other scripts ever add to the mutations node

        for entry in DB_root:
            mutations_node = entry.find('mutations')
            if mutations_node != None:
                DB_comparison_entry_node = etree.SubElement(DB_comparison_root, 'entry')
                DB_comparison_entry_node.append(copy.deepcopy(mutations_node))

        DB_comparison_string = etree.tostring(DB_comparison_root, pretty_print=True)
        DB_comparison_string = '\n'.join( DB_comparison_string.splitlines()[1:] ) # Ignore the first line, which contains datestamps (don't want to include these in the comparison)

        for entry in DBold_root:
            mutations_node = entry.find('mutations')
            if mutations_node != None:
                DBold_comparison_entry_node = etree.SubElement(DBold_comparison_root, 'entry')
                DBold_comparison_entry_node.append(copy.deepcopy(mutations_node))

        DBold_comparison_string = etree.tostring(DBold_comparison_root, pretty_print=True)
        DBold_comparison_string = '\n'.join( DBold_comparison_string.splitlines()[1:] ) # Ignore the first line, which contains datestamps (don't want to include these in the comparison)

        # Now compare the two comparison strings using GNU diff
        diff_output = clab.DB.diff_DB_comparison_strings(DBold_comparison_string, DB_comparison_string)

        if len(diff_output) > 0:
            print 'Comparison of latest cBioPortal-derived data with data in %s indicates changes. File will be rewritten with new data. Lines in diff comparison: %s' % (DB_out_filepath, len(diff_output))
            data_modified = True
            if verbose:
                print DB_comparison_string.split('\n')[0:20]
                print DBold_comparison_string.split('\n')[0:20]

        else:
            print 'Comparison of latest cBioPortal-derived data with data in %s indicates no changes. File will be rewritten with updated gather_cbioportal_last_run attrib, but other data will not be modified.' % DB_out_filepath
            data_modified = False


#==============================================================================
# If staging and there have been modifications, update date_modified
#==============================================================================

if run_mode == 'stage' and data_modified:
    DB_root.set('gather_cbioportal_last_modif', datestamp)

# =======================
# write the XML DB
# =======================
if run_mode == 'stage' and data_modified:
    clab.DB.writeDB(DB_root, DB_out_filepath)

elif run_mode == 'stage' and not data_modified:
    DBold_root = etree.parse(DB_out_filepath, parser).getroot()
    DBold_root.set('gather_pdb_last_run', datestamp)
    clab.DB.writeDB(DBold_root, DB_out_filepath)

elif run_mode == 'dev':
    clab.DB.writeDB(DB_root, DB_out_filepath)

print ''
print 'Done.'



