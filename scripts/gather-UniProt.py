#!/usr/bin/env python
#
# Extract pertinent data from UniProt XML document and store as an XML database.
#
# PDB cross-refs are added if the structure is deemed to contain the PK domain.
# This is determined from whether the PDB sequence span can include the span of
# the PK domain less 30 residues at each end. Only PDB ID, chain ID, and resi
# span are added. Use gather-pdb.py to add further info from sifts.
#
# New UniProt XML document downloaded only if existing one is > 7 days old, or
# if --forcedl flag is used.
#
# Daniel L. Parton <partond@mskcc.org> - 7 Mar 2013
#
# Dependencies: diff (GNU)
#

#==============================================================================
# IMPORTS
#==============================================================================

import sys, datetime, os, copy
import TargetExplorer
from lxml import etree

TargetExplorer.core.check_correct_working_dir()

#==============================================================================
# PARAMETERS
#==============================================================================

if '-stage' in sys.argv:
    run_mode = 'stage'
elif '-dev' in sys.argv:
    run_mode = 'dev'
else:
    run_mode = 'nowrite'

print 'Running in mode: %s' % run_mode

database_dir = 'database'
external_data_dir = os.path.join('external-data')
if not os.path.exists(external_data_dir):
    os.mkdir(external_data_dir)
external_data_metadata_filepath = os.path.join(external_data_dir, 'metadata.xml')
uniprot_data_dir = os.path.join(external_data_dir, 'UniProt')

if not os.path.exists(uniprot_data_dir):
    os.mkdir(uniprot_data_dir)

uniprot_xml_out_filepath = os.path.join(uniprot_data_dir, 'uniprot-search.xml')

DB_out_filename = 'database-%(run_mode)s.xml' % vars()
DB_out_filepath = os.path.join(database_dir, DB_out_filename)

uniprot_search_string = '(taxonomy:40674 AND domain:"protein kinase") AND reviewed:yes'
uniprot_search_string_query = '?query=%28taxonomy%3A40674+AND+domain%3A%22protein+kinase%22%29+AND+reviewed%3Ayes&sort=score&format=xml' # a bit crude, but urllib.urlencode might not be much better
# 1GQ5 is referenced by kinase P16234. The kinase is not in the actual structure.
ignore_uniprot_pdbs = ['1GQ5']

if '--forcedl' in sys.argv:
    force_uniprot_download = True
else:
    force_uniprot_download = False

now = datetime.datetime.utcnow()
datestamp = now.strftime(TargetExplorer.DB.datestamp_format_string)

parser = etree.XMLParser(remove_blank_text=True, huge_tree=True)

days_elapsed_for_force_download = 7

#==============================================================================
# RETRIEVE DATA FROM UNIPROT AND STORE TO LOCAL FILE
#==============================================================================

# First check if UniProt XML document already exists
if os.path.exists(uniprot_xml_out_filepath):
    print 'UniProt XML document found at:', uniprot_xml_out_filepath
# if not, search the UniProt database and save an XML document
else:
    print 'UniProt XML document not found.'
    print 'Retrieving new XML document from UniProt website.'
    new_xml_text = TargetExplorer.UniProt.retrieve_uniprot(uniprot_search_string_query)
    print 'Saving new XML document as:', uniprot_xml_out_filepath
    with open(uniprot_xml_out_filepath, 'w') as uniprot_xml_file:
        uniprot_xml_file.write(new_xml_text + '\n')
    TargetExplorer.UniProt.update_metadata_uniprot_search(datestamp, uniprot_xml_out_filepath)

# Read in the UniProt XML document
print 'Reading UniProt XML document:', uniprot_xml_out_filepath
uniprot_xml = etree.parse(uniprot_xml_out_filepath, parser).getroot()

# Check when the UniProt XML document was retrieved and call retrieve_uniprot() if time elapsed is more than days_elapsed_for_force_download
# TODO external_data_metadata_file may not exist yet - should handle this elegantly
metadata_root = etree.parse(external_data_metadata_filepath, parser).getroot()
retrieved = metadata_root.find('UniProt/uniprot_search').get('datestamp')
retrieved = datetime.datetime.strptime(retrieved, TargetExplorer.DB.datestamp_format_string) # turn into datetime object
now = datetime.datetime.now()
time_elapsed = now - retrieved
print 'UniProt XML document was retrieved: %s (%s days ago)' % (retrieved.strftime('%Y-%m-%d'), time_elapsed.days)
if (time_elapsed.days > days_elapsed_for_force_download) or (force_uniprot_download == True):
    if time_elapsed.days > days_elapsed_for_force_download:
        print 'UniProt XML document more than %d days old.' % days_elapsed_for_force_download
    if force_uniprot_download == True:
        print 'Forcing retrieval of new UniProt XML document.'
        download_new_uniprot_xml = True
    else:
        while True:
            user_response = raw_input('Suggest retrieving new XML document from UniProt. Proceed? [y] ')
            if user_response in ['y', '']:
                download_new_uniprot_xml = True
                break
            elif user_response == 'n':
                download_new_uniprot_xml = False
                break
            else:
                print 'User input not understood. Please try again.'

    if download_new_uniprot_xml:
        print 'Retrieving new XML document from UniProt...'
        old_xml = uniprot_xml
        new_xml_text = TargetExplorer.UniProt.retrieve_uniprot(uniprot_search_string_query)
        new_xml = etree.fromstring(new_xml_text, parser)

        # Print some basic statistics on the differences between the new and old XML documents
        print '\n=== COMPARISON OF NEW AND OLD XML DOCUMENTS ==='
        len_diff = TargetExplorer.UniProt.print_uniprot_xml_comparison(new_xml, old_xml)
        print ''

        # If no differences, continue (but update external-data/metadata.xml)
        if len_diff == 0:
            print 'No difference between new and old XML documents - continuing...'
            TargetExplorer.UniProt.update_metadata_uniprot_search(datestamp, uniprot_xml_out_filepath)
            pass

        # Otherwise prompt as to whether or not to write to file
        else:
            while True:
                user_response = raw_input('Write new XML document to %s? [y] ' % uniprot_xml_out_filepath)
                if user_response in ['y', '']:
                    print 'Saving new XML document as:', uniprot_xml_out_filepath
                    with open(uniprot_xml_out_filepath, 'w') as uniprot_xml_out_file:
                        uniprot_xml_out_file.write(new_xml_text)
                    TargetExplorer.UniProt.update_metadata_uniprot_search(datestamp, uniprot_xml_out_filepath)
                    uniprot_xml = new_xml
                    break
                elif user_response in ['n']:
                    print 'Continuing...'
                    break
                else:
                    print 'User input not understood. Please try again.'

print ''

uniprot_entries = uniprot_xml.findall('entry')
nuniprot_entries = len(uniprot_entries)
# Note that xpath querying is case-sensitive
print 'Number of entries in UniProt XML document:', nuniprot_entries
print 'Number of domains:' , len( uniprot_xml.xpath('./entry/feature[@type="domain"]') )
print 'Number of domains containing "kinase":' , len( uniprot_xml.xpath('./entry/feature[@type="domain"][contains(@description,"kinase")]') )
print 'Number of domains containing "Kinase":' , len( uniprot_xml.xpath('./entry/feature[@type="domain"][contains(@description,"Kinase")]') )
print 'Number of domains containing "Protein kinase":' , len( uniprot_xml.xpath('./entry/feature[@type="domain"][contains(@description,"Protein kinase")]') )
print 'Number of domains which are (exactly) "Protein kinase":' , len( uniprot_xml.xpath('./entry/feature[@type="domain"][@description="Protein kinase"]') )
print '= Domains which contain "Protein kinase" but do not equal "Protein kinase" are of the following types: =\nProtein kinase 1\nProtein kinase 2\nProtein kinase; truncated\nProtein kinase; inactive'
print '= Domains which contain "kinase" but do not equal "Protein kinase": ='
print 'Number of domains containing "Alpha-type protein kinase":' , len( uniprot_xml.xpath('./entry/feature[@type="domain"][contains(@description,"Alpha-type protein kinase")]') )
print 'Number of domains containing "AGC-kinase C-terminal":' , len( uniprot_xml.xpath('./entry/feature[@type="domain"][contains(@description,"AGC-kinase C-terminal")]') )
print 'Number of domains containing "Guanylate kinase-like":' , len( uniprot_xml.xpath('./entry/feature[@type="domain"][contains(@description,"Guanylate kinase-like")]') )
print 'Keeping only domains containing "Protein kinase"... (case sensitive)'

print ''

#==============================================================================
# CREATE NEW DATABASE XML TREE FROM DOWNLOADED UNIPROT DATA
#==============================================================================

# Root node to which all data will be added
DB_root = etree.Element('database')

# Iterate through each kinase from the UniProt XML document
for k in range(nuniprot_entries):

    # Create an entry element for each UniProt entry
    DBentry = etree.SubElement(DB_root, 'entry')
    # And a UniProt subelement
    DBentry_uniprot = etree.SubElement(DBentry, 'UniProt')

    # = IDs and names =
    AC = uniprot_entries[k].findtext('./accession')
    DBentry_uniprot.set('AC', AC)
    entry_name = uniprot_entries[k].findtext('./name')
    DBentry_uniprot.set('entry_name', entry_name)
    rec_name = uniprot_entries[k].findtext('./protein/recommendedName/fullName')
    etree.SubElement(DBentry_uniprot, 'protein_recommended_name').text = rec_name
    gene_name_nodes = uniprot_entries[k].findall('./gene/name')
    DBentry_uniprot_gene_names = etree.SubElement(DBentry_uniprot, 'gene_names')
    for gene_name_node in gene_name_nodes:
        DBgene_name_node = etree.SubElement(DBentry_uniprot_gene_names, 'gene_name')
        for key in gene_name_node.keys():
            DBgene_name_node.set(key, gene_name_node.get(key))
            DBgene_name_node.text = gene_name_node.text

    # = Date entry was last modified in UniProt =
    modified = uniprot_entries[k].attrib['modified']
    DBentry_uniprot.set('last_UniProt_update', modified)

    # XXX exception for SG196_HUMAN, which does not have protein kinase activity, and acts as a mannose kinase instead
    if entry_name == 'SG196_HUMAN':
        print 'Removing kinase as it does not have protein kinase activity (instead acts as a mannose kinase):', AC
        DBentry.getparent().remove(DBentry)
        continue

    # = Functions, disease associations =
    functions_node = etree.SubElement(DBentry_uniprot, 'functions')
    disease_associations_node = etree.SubElement(DBentry_uniprot, 'disease_associations')
    for x in uniprot_entries[k].findall('./comment[@type="function"]'):
        etree.SubElement(functions_node, 'function').text = TargetExplorer.core.twrap( x.findtext('./text') )
    for x in uniprot_entries[k].findall('./comment[@type="disease"]'):
        etree.SubElement(disease_associations_node, 'disease_association').text = TargetExplorer.core.twrap( x.findtext('./text') )

    # = Isoforms =
    # Canonical isoform is given the attrib type="displayed", meaning that the sequence is displayed in the HTML version of the entry
    # Example alt isoform:
    #     <isoform>
    #         <id>P00519-2</id>
    #         <name>IB</name>
    #         <sequence type="described" ref="VSP_004957"/>
    #         <note>Contains a N-myristoyl glycine at position 2.</note>
    #     </isoform>

    DBentry_uniprot_isoforms_node = etree.SubElement(DBentry_uniprot, 'isoforms')
    DBentry_uniprot_canonical_isoform_node = etree.SubElement(DBentry_uniprot_isoforms_node, 'canonical_isoform')
    for uniprot_isoform_node in uniprot_entries[k].findall('isoform'):
        isoform_AC = uniprot_isoform_node.findtext('id')
        notes = uniprot_isoform_node.findall('note')
        if uniprot_isoform_node.get('type') == 'displayed':
            DB_isoform_node = DBentry_uniprot_canonical_isoform_node
        else:
            DB_isoform_node = etree.SubElement(DBentry_uniprot_isoforms_node, 'alt_isoform')

        DB_isoform_node.set('AC', isoform_AC)
        for note in notes:
            DB_isoform_node.append(copy.deepcopy(note))

    # = Canonical sequence =
    # Returned UniProt XML contains sequence data only for the canonical isoform
    uniprot_canonical_sequence_node = uniprot_entries[k].find('./sequence[@length][@mass]')
    canonical_sequence = ''.join(uniprot_canonical_sequence_node.text.split())
    DBentry_uniprot_canonical_isoform_node.append( copy.deepcopy(uniprot_canonical_sequence_node) )

    # = UniProt "Protein kinase" domain annotations =
    # XXX TODO Generalize

    domains = uniprot_entries[k].xpath('./feature[@type="domain"][contains(@description,"Protein kinase")]')

    # XXX exceptions
    # These are the entries for which "Protein kinase" domains are known to be not found (case sensitive):
    # kinases_with_no_PK_domain = ['ALPK1_HUMAN', 'ALPK2_HUMAN', 'ALPK3_HUMAN', 'EF2K_HUMAN', 'TRPM6_HUMAN', 'TRPM7_HUMAN']
    # These are all alpha-kinases, which have no identity with typical protein kinases.
    # These kinases will therefore be deleted from the database.
    if len(domains) < 1:
        print 'Removing kinase as it does not possess a domain annotation containing "Protein kinase":', AC
        DBentry.getparent().remove(DBentry)
        continue
    # In cases where > 1 PK domain is found, add a warning to the DB entry. In some cases, a pseudokinase is present - these domains are not added.
    warnings_node = etree.Element('warnings')
    if len(domains) > 1:
        if uniprot_entries[k].findtext('name') == 'E2AK4_HUMAN':
            etree.SubElement(warnings_node,'warning').text = 'Kinase is annotated in UniProt wth both "Protein kinase 1" and "Protein kinase 2". "Protein kinase 1" is considered to be a pseudokinase domain. "Protein kinase 2" is considered active. Only the active PK domain is included in this DB.'
            domains.pop(0)
        elif uniprot_entries[k].findtext('name') in ['JAK1_HUMAN','JAK2_HUMAN','JAK3_HUMAN']:
            etree.SubElement(warnings_node,'warning').text = 'Kinase is annotated in UniProt wth both "Protein kinase 1" and "Protein kinase 2". Janus (Jak) tyrosine kinases (JAK1, JAK2 and JAK3) each contain a tyrosine kinase domain adjacent to a catalytically inactive pseudokinase domain. The pseudokinase domain interacts with and negatively regulates the active domain. The pseudokinase domain is the first one in the sequence. Only the active PK domain is included in this DB.'
            domains.pop(0)
        elif uniprot_entries[k].findtext('name') in ['KS6A1_HUMAN','KS6A2_HUMAN','KS6A3_HUMAN','KS6A4_HUMAN','KS6A5_HUMAN','KS6A6_HUMAN']:
            etree.SubElement(warnings_node,'warning').text = 'Kinase is annotated in UniProt wth both "Protein kinase 1" and "Protein kinase 2". Upon extracellular signal or mitogen stimulation, phosphorylated at Thr-573 in the C-terminal kinase domain (CTKD) by MAPK1/ERK2 and MAPK3/ERK1. The activated CTKD then autophosphorylates Ser-380, allowing binding of PDPK1, which in turn phosphorylates Ser-221 in the N-terminal kinase domain (NTKD) leading to the full activation of the protein and subsequent phosphorylation of the substrates by the NTKD. Both PK domains are included in this DB.'
        elif uniprot_entries[k].findtext('name') == 'KS6C1_HUMAN':
            etree.SubElement(warnings_node,'warning').text = 'Kinase is annotated in UniProt wth both "Protein kinase 1" and "Protein kinase 2". The first protein kinase domain appears to be a pseudokinase domain as it does not contain the classical characteristics, such as the ATP-binding motif, ATP-binding site and active site. Only "Protein kinase 2" is included in this DB.'
            domains.pop(0)
        elif uniprot_entries[k].findtext('name') == 'OBSCN_HUMAN':
            etree.SubElement(warnings_node,'warning').text = 'Kinase is annotated in UniProt wth both "Protein kinase 1" and "Protein kinase 2". Neither are described as pseudokinases, although are not specifically described as catalytically active either. Both PK domains are included in this DB.'
        elif uniprot_entries[k].findtext('name') == 'SPEG_HUMAN':
            etree.SubElement(warnings_node,'warning').text = 'Kinase is annotated in UniProt wth both "Protein kinase 1" and "Protein kinase 2". Neither are described as pseudokinases. Both PK domains are included in this DB.'
        elif uniprot_entries[k].findtext('name') == 'TAF1_HUMAN':
            etree.SubElement(warnings_node,'warning').text = 'Kinase is annotated in UniProt wth both "Protein kinase 1" and "Protein kinase 2". Neither are described as pseudokinases. Both PK domains are included in this DB.'
        elif uniprot_entries[k].findtext('name') == 'TYK2_HUMAN':
            etree.SubElement(warnings_node,'warning').text = 'Kinase is annotated in UniProt wth both "Protein kinase 1" and "Protein kinase 2". Neither are described as pseudokinases. Both PK domains are included in this DB.'
        else:
            raise Exception, 'More than 1 domain found containing "Protein kinase". Please check the following kinase and adjust the script: %s' % entry_name
    # And a couple of cases with one PK domain which are considered inactive. These kinase entries are removed completely.
    if domains[0].attrib['description'] == 'Protein kinase; truncated':
        # PLK5_HUMAN. Kinase considered inactive. Protein kinase domain is truncated. Remove it.
        print 'Removing kinase as PK domain is truncated and considered inactive:', AC
        DBentry.getparent().remove(DBentry)
        continue
    elif domains[0].attrib['description'] == 'Protein kinase; inactive':
        # PTK7_HUMAN. Kinase considered inactive. Remove it.
        print 'Removing kinase as PK domain is considered inactive:', AC
        DBentry.getparent().remove(DBentry)
        continue

    # Finally, add the domains to the new database
    DBentry_domains_node = etree.SubElement(DBentry_uniprot, 'domains')
    for x_iter,x in enumerate(domains):
        # First calculate the PK domain length and sequence
        pk_description = x.get('description')
        pk_begin = int( x.find('./location/begin').attrib['position'] )
        pk_end = int( x.find('./location/end').attrib['position'] )
        pk_length = pk_end - pk_begin + 1
        #domain = copy.deepcopy(x)
        domain = etree.Element('domain')
        domain.set('description', pk_description)
        domain.set('begin', str(pk_begin))
        domain.set('end', str(pk_end))
        domain.set('length', str(pk_length))
        domain.set('domainID', str(x_iter))
        domain.set('targetID', (entry_name + '_D' + str(x_iter)))

        #location = domain.find('./location')
        #etree.SubElement(location, 'length').text = str(pk_length)
        domain_sequence = TargetExplorer.core.seqwrap(canonical_sequence[pk_begin-1:pk_end])
        etree.SubElement(domain, 'sequence').text = '\n' + domain_sequence
        DBentry_domains_node.append(domain)

    # = References to other DBs =
    # NCBI Gene
    GeneIDs = [x.get('id') for x in uniprot_entries[k].findall('./dbReference[@type="GeneID"]')]
    # XXX: exceptions for kinases which have no GeneIDs annotated; LMTK3 RefSeq status is PROVISIONAL; RIPK4 presumably RefSeq sequence is not an exact match; SIK3 RefSeq status is VALIDATED
    # Will add these manually, since we are mainly using GeneID to collect publications currently
    if entry_name == 'LMTK3_HUMAN':
        GeneIDs = ['114783']
    if entry_name == 'RIPK4_HUMAN':
        GeneIDs = ['54101']
    if entry_name == 'SIK3_HUMAN':
        GeneIDs = ['23387']
    if len(GeneIDs) > 0:
        NCBI_Gene_node = etree.SubElement(DBentry, 'NCBI_Gene')
    for GeneID in GeneIDs:
        # XXX: exceptions for SGK3_HUMAN and TNI3K_HUMAN, which have two GeneIDs annotated; in each case, one is a readthrough fusion protein - ignore these GeneIDs
        if GeneID in ['100533105', '100526835']:
            continue
        NCBI_Gene_entry_node = etree.SubElement(NCBI_Gene_node, 'entry')
        NCBI_Gene_entry_node.set('ID', GeneID)

    # Ensembl
    EnsemblGeneIDs = uniprot_entries[k].findall('./dbReference[@type="Ensembl"]/property[@type="gene ID"]')
    EnsemblGeneIDs_set = set( [ id.attrib['value'] for id in EnsemblGeneIDs ] )
    DB_Ensembl_node = etree.SubElement(DBentry, 'Ensembl')
    for EnsemblGeneID in EnsemblGeneIDs_set:
        etree.SubElement(DB_Ensembl_node, 'GeneID').text = EnsemblGeneID

    # HGNC
    HGNC_dbRefs = uniprot_entries[k].findall('./dbReference[@type="HGNC"]')
    if len(HGNC_dbRefs) > 0:
        HGNC_element = etree.SubElement(DBentry, 'HGNC')
        for HGNC_dbRef in HGNC_dbRefs:
            ID = HGNC_dbRef.get('id')
            Approved_Symbol = HGNC_dbRef.find('property[@type="gene designation"]').get('value')
            HGNC_entry_element = etree.SubElement(HGNC_element, 'entry')
            HGNC_entry_element.set('ID', ID)
            HGNC_entry_element.set('Approved_Symbol', Approved_Symbol)

    # = Family information =
    similarity_comments = uniprot_entries[k].xpath('./comment[@type="similarity"]')
    family_found = False
    for s in similarity_comments:
        for f in TargetExplorer.UniProt.kinase_family_uniprot_similarity_text.keys():
            if f in s.findtext('text'):
                DBentry_uniprot.set('family', TargetExplorer.UniProt.kinase_family_uniprot_similarity_text[f])
                family_found = True
    if family_found == False:
        DBentry_uniprot.set('family', '')

    # = PDB entries (from UniProt XML) =
    pdbs = uniprot_entries[k].findall('./dbReference[@type="PDB"]')
    if len(pdbs) > 0:
        DBentry_pdb_node = etree.SubElement(DBentry, 'PDB')
    for p in pdbs:
        # Only keep XRC structures (no NMR or Model)
        if p.find('property[@type="method"]') == None:
            if p.attrib['id'] == '2LV6':
                continue  # 2LV6 has no method listed - it is actually an NMR structure, including only a very short fragment of the kinase, outside the PK domain
        elif p.find('property[@type="method"]').attrib['value'] == 'X-ray':
            pdb_structure_node = etree.Element('structure')
            pdbid = p.attrib['id']
            if pdbid in ignore_uniprot_pdbs:
                continue
            pdb_structure_node.set('ID', pdbid)
            resolution = p.find('property[@type="resolution"]').attrib['value']
            pdb_structure_node.set('resolution', resolution)
            chains_span_str = p.find('property[@type="chains"]').attrib['value']
            chains_span = TargetExplorer.UniProt.parse_uniprot_pdbref_chains(chains_span_str)
            chains_added = 0
            for c in chains_span.keys():
                chains = etree.Element('chain')
                chains.set('ID', c)
                pdb_begin = chains_span[c][0]
                pdb_end = chains_span[c][1]
                # Use the begin and end info to decide if this pdb chain includes the pk_domain. But we will get other sequence info from sifts XML files, using gather-pdb.py
                # Have to check against each PK domain
                for x_iter,x in enumerate(DBentry_domains_node):
                    pk_begin = int(x.get('begin'))
                    pk_end = int(x.get('end'))
                    if (pdb_begin < pk_begin+30) & (pdb_end > pk_end-30):
                        chains.set('domainID', str(x_iter))
                        chains.set('begin', str(pdb_begin))
                        chains.set('end', str(pdb_end))
                        pdb_structure_node.append(chains)
                        chains_added += 1
                    else:
                        continue

                if chains_added > 0:
                    DBentry_pdb_node.append(pdb_structure_node)

    # = Add the warnings node last (only if it contains any warnings) = #
    if len(warnings_node) > 0:
        DBentry.append(warnings_node)

print ''

#==============================================================================
# IF STAGING, UPDATE DATE_RUN
#==============================================================================

if run_mode == 'stage':
    DB_root.set('gather_uniprot_last_run', datestamp)

#==============================================================================
# IF STAGING, COMPARE NEW AND OLD DBS
#==============================================================================

data_modified = False

if run_mode == 'stage':

    if not os.path.exists(DB_out_filepath):
        data_modified = True

    else:
        # Parse the old DB
        DBold_root = etree.parse(DB_out_filepath, parser).getroot()

        # First a quick check to see whether the numbers of entries match
        if len(DBold_root) != len(DB_root):
            print 'Comparison of latest UniProt data with data in %s indicates changes. File will be overwritten with new data. Number of entries differs by: %d' % (DB_out_filepath, len(DB_root) - len(DBold_root))
            data_modified = True

        # If the numbers of entries match, then proceeed to a full comparison of UniProt-derived data using diff
        else:
            DB_comparison_string = etree.tostring(DB_root, pretty_print=True)
            DB_comparison_string = '\n'.join( DB_comparison_string.splitlines()[1:] ) # Ignore the first line, which contains datestamps (don't want to include these in the comparison)

            # A new tree will be built for DBold, containing only the data gathered by this script
            DBold_comparison_root = etree.Element('database')

            for DBold_entry in DBold_root:
                # make list of nodes containing only data derived from this script
                DBold_entry_uniprot_node = DBold_entry.find('UniProt')
                DBold_entry_HGNC_node = DBold_entry.find('HGNC')
                DBold_entry_Ensembl_node = DBold_entry.find('Ensembl')
                DBold_entry_NCBI_Gene_node = DBold_entry.find('NCBI_Gene')
                DBold_entry_PDB_node = DBold_entry.find('PDB')
                DBold_entry_warnings_node = DBold_entry.find('warnings')

                # For the NCBI Gene node, we need to remove everything apart from:
                # <NCBI_Gene>
                #   <entry ID= >
                if DBold_entry_NCBI_Gene_node != None:
                    entry_node = DBold_entry_NCBI_Gene_node.find('entry')
                    for child in entry_node.getchildren():
                        entry_node.remove(child)
                    for attrib_key in entry_node.keys():
                        if attrib_key not in ['ID']:
                            entry_node.attrib.pop(attrib_key)

                # For the PDB node, we need to remove everything apart from:
                # <PDB>
                #   <structure ID= resolution= >
                #     <chain ID= domainID= begin= end= >
                if DBold_entry_PDB_node != None:
                    structure_nodes = DBold_entry_PDB_node.findall('structure')
                    for structure_node in structure_nodes:

                        for attrib_key in structure_node.keys():
                            if attrib_key not in ['ID', 'resolution']:
                                structure_node.attrib.pop(attrib_key)

                        chain_nodes = structure_node.findall('chain')

                        for chain_node in chain_nodes:

                            for child in chain_node.getchildren():
                                chain_node.remove(child)

                            for attrib_key in chain_node.keys():
                                if attrib_key not in ['ID', 'domainID', 'begin', 'end']:
                                    chain_node.attrib.pop(attrib_key)

                # deepcopy these nodes to a new entry node
                # XXX important to add these in document order
                DBold_comparison_entry_node = etree.SubElement(DBold_comparison_root, 'entry')
                DBold_comparison_entry_node.append(copy.deepcopy(DBold_entry_uniprot_node))
                if DBold_entry_NCBI_Gene_node != None:
                    DBold_comparison_entry_node.append(DBold_entry_NCBI_Gene_node)
                if DBold_entry_Ensembl_node != None:
                    DBold_comparison_entry_node.append(copy.deepcopy(DBold_entry_Ensembl_node))
                if DBold_entry_HGNC_node != None:
                    DBold_comparison_entry_node.append(copy.deepcopy(DBold_entry_HGNC_node))
                if DBold_entry_PDB_node != None:
                    DBold_comparison_entry_node.append(DBold_entry_PDB_node)
                if DBold_entry_warnings_node != None:
                    DBold_comparison_entry_node.append(DBold_entry_warnings_node)

            # convert the new entry node to string
            DBold_comparison_string = etree.tostring(DBold_comparison_root, pretty_print=True)
            DBold_comparison_string = '\n'.join( DBold_comparison_string.splitlines()[1:] ) # Ignore the first line, which contains datestamps (don't want to include these in the comparison)

            # Now compare the two comparison strings using GNU diff
            diff_output = TargetExplorer.DB.diff_DB_comparison_strings(DBold_comparison_string, DB_comparison_string)

            if len(diff_output) > 0:
                print 'Comparison of latest UniProt data with data in %s indicates changes. File will be overwritten. Lines in diff comparison: %s' % (DB_out_filepath, len(diff_output))
                data_modified = True
            else:
                print 'Comparison of latest UniProt data with data in %s indicates no changes. File will be rewritten with updated gather_uniprot_late_run attrib, but other data will not be modified.' % DB_out_filepath
                data_modified = False


#==============================================================================
# IF STAGING AND THERE HAVE BEEN MODIFICATIONS, UPDATE DATE_MODIFIED
#==============================================================================

if run_mode == 'stage' and data_modified:
    DB_root.set('gather_uniprot_last_modif', datestamp)

#==============================================================================
# WRITE DATABASE TO FILE
#==============================================================================

if run_mode == 'stage' and data_modified:
    TargetExplorer.DB.writeDB(DB_root, DB_out_filepath)

elif run_mode == 'stage' and not data_modified:
    DBold_root = etree.parse(DB_out_filepath, parser).getroot()
    DBold_root.set('gather_uniprot_last_run', datestamp)
    TargetExplorer.DB.writeDB(DBold_root, DB_out_filepath)

elif run_mode == 'dev':
    TargetExplorer.DB.writeDB(DB_root, DB_out_filepath)

print ''
print 'Done.'

