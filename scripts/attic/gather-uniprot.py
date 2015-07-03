# Extract pertinent data from UniProt XML document and store as kinDB XML file.
#
# PDB cross-refs are added if the structure is deemed to contain the PK domain.
# This is determined from whether the PDB sequence span can include the span of
# the PK domain less 30 residues at each end. Only PDB ID, chain ID, and resi
# span are added. Use gather-protein_databank.py to add further info from sifts.
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

import sys, datetime, os
from choderalab.uniprot import retrieve_uniprot, parse_uniprot_pdbref_chains, kinase_family_uniprot_similarity_text
from choderalab.core import twrap, seqwrap
from lxml import etree
from copy import deepcopy
from subprocess import call

#==============================================================================
# PARAMETERS
#==============================================================================

database_dir = 'database'
external_data_dir = 'external-data'
uniprot_data_dir = os.path.join(external_data_dir, 'UniProt')

if not os.path.exists(uniprot_data_dir):
    os.mkdir(uniprot_data_dir)

uniprot_xml_out_filepath = os.path.join(uniprot_data_dir, 'uniprot-human_kinases.xml')
database_out_filepath = os.path.join(database_dir, 'kinDB.xml')

uniprot_search_string = '(organism:9606 AND domain:"protein kinase") AND reviewed:yes'
uniprot_search_string_query = '?query=%28organism%3A9606+AND+domain%3A%22protein+kinase%22%29+AND+reviewed%3Ayes&sort=score&format=xml' # a bit crude, but urllib.urlencode might not be much better
# 1GQ5 is referenced by kinase P16234. The kinase is not in the actual structure.
ignore_uniprot_pdbs = ['1GQ5']

if '--forcedl' in sys.argv:
    force_uniprot_download = True
else:
    force_uniprot_download = False

#==============================================================================
# MAIN
#==============================================================================

# First check if UniProt XML document already exists
if os.path.exists(uniprot_xml_out_filepath):
    print 'UniProt XML document found at:', uniprot_xml_out_filepath
# if not, search the UniProt database and save an XML document
else:
    print 'UniProt XML document not found.'
    print 'Retrieving new XML document from UniProt website.'
    new_xml = retrieve_uniprot(uniprot_search_string_query)
    print 'Saving new XML document as:', uniprot_xml_out_filepath
    with open(uniprot_xml_out_filepath, 'w') as uniprot_xml_file:
        uniprot_xml_file.write(new_xml + '\n')

# Read in the UniProt XML document
print 'Reading UniProt XML document:', uniprot_xml_out_filepath
parser = etree.XMLParser(remove_blank_text=True)
uniprot_xml = etree.parse(uniprot_xml_out_filepath, parser).getroot()

# Check when the UniProt XML document was retrieved and call retrieve_uniprot() if more than a week old
retrieved = uniprot_xml.find('.').attrib['retrieve_date']
retrieved = datetime.datetime.strptime(retrieved, '%Y-%m-%d') # turn into datetime object
now = datetime.datetime.now()
time_elapsed = now - retrieved
print 'UniProt XML document was retrieved: %s (%s days ago)' % (retrieved.strftime('%Y-%m-%d'), time_elapsed.days)
if (time_elapsed.days > 7) or (force_uniprot_download == True):
    if time_elapsed.days > 7:
        print 'UniProt XML document more than 7 days old.'
    elif force_uniprot_download == True:
        print 'Forcing download of new UniProt XML document.'
    print 'Retrieving new XML document from UniProt website.'
    new_xml = retrieve_uniprot(uniprot_search_string_query)
    # First save as a tmp file and carry out a diff with the current XML document
    # There is a python library libdiff, but it is extremely slow compared to the GNU tool
    # ?TODO instead of printing diff, save to a file, and print a comparison of stats (number of kinases; number of kinases with pdbs; number of kinases with pdbs and disease associations)
    print 'Conducting diff against current XML document.'
    new_uniprot_xml_path = 'tmp-new-uniprot.xml'
    with open(new_uniprot_xml_path, 'w') as new_uniprot_xml_file:
        new_uniprot_xml_file.write(new_xml + '\n')
    with open('tmp-diff', 'w') as diff_file:
        call(['diff', '--ignore-matching-lines=<uniprot retrieve_date', uniprot_xml_out_filepath, new_uniprot_xml_path], stdout=diff_file)
    with open('tmp-diff', 'r') as diff_file:
        diff = diff_file.readlines()
    if len(diff) > 0:
        print 'Differences found:\n==========\n'
        print ''.join(diff)
        print '\n==========\n'
    else:
        print '\n==========\nNo differences found. Continuing anyway.\n==========\n'
    print 'Saving new XML document as:', uniprot_xml_out_filepath
    os.rename('tmp-new-uniprot.xml', uniprot_xml_out_filepath)
    os.remove('tmp-diff')
    uniprot_xml = etree.parse(uniprot_xml_out_filepath, parser).getroot()

uniprot_kinases = uniprot_xml.xpath('./entry')
nkinases = len(uniprot_kinases)
# Note that xpath querying is case-sensitive
print 'Number of entries in UniProt XML document:', nkinases
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

# Data will be added to this XML root node
root = etree.Element('kinome')

# Iterate through each kinase from the UniProt XML document
for k in range(nkinases):
    kinase = etree.SubElement(root, 'kinase')
    kinase_uniprot = etree.SubElement(kinase, 'uniprot')

    # = Date entry was last modified in UniProt =
    modified = uniprot_kinases[k].attrib['modified']
    kinase_uniprot.set('last_uniprot_update', modified)

    # = IDs and names =
    AC = uniprot_kinases[k].findtext('./accession')
    kinase_uniprot.set('AC', AC)
    entry_name = uniprot_kinases[k].findtext('./name')
    kinase_uniprot.set('entry_name', entry_name)
    rec_name = uniprot_kinases[k].findtext('./protein/recommendedName/fullName')
    etree.SubElement(kinase_uniprot, 'protein_recommended_name').text = rec_name
    gene_name = uniprot_kinases[k].findtext('./gene/name[@type="primary"]')
    kinase_uniprot.set('primary_gene_name', gene_name)

    # XXX exception for SG196_HUMAN, which does not have protein kinase activity, and acts as a mannose kinase instead
    if entry_name == 'SG196_HUMAN':
        print 'Removing kinase as it does not have protein kinase acitivy (instead acts as a mannose kinase):', AC
        kinase.getparent().remove(kinase)
        continue

    # = Function, disease, complete sequence =
    for x in uniprot_kinases[k].xpath('./comment[@type="function"]'):
        etree.SubElement(kinase_uniprot, 'function').text = twrap( x.findtext('./text') )
    for x in uniprot_kinases[k].xpath('./comment[@type="disease"]'):
        etree.SubElement(kinase_uniprot, 'disease_association').text = twrap( x.findtext('./text') )
    for x in uniprot_kinases[k].xpath('./sequence[@length]'):
        kinase_uniprot.append(deepcopy(x)) # Can assume only one sequence per UniProt entry (if there are multiple isoforms, this sequence usually corresponds to the most common isoform)
        sequence = ''.join(x.text.split())

    # = UniProt "Protein kinase" domain annotations =
    PK_domains = uniprot_kinases[k].xpath('./feature[@type="domain"][contains(@description,"Protein kinase")]')
    # XXX exceptions
    # These are the entries for which "Protein kinase" domains are known to be not found (case sensitive):
    # kinases_with_no_PK_domain = ['ALPK1_HUMAN', 'ALPK2_HUMAN', 'ALPK3_HUMAN', 'EF2K_HUMAN', 'TRPM6_HUMAN', 'TRPM7_HUMAN']
    # These are all alpha-kinases, which have no identity with typical protein kinases.
    # These kinases will therefore be deleted from kinDB.
    if len(PK_domains) < 1:
        print 'Removing kinase as it does not possess a domain annotation containing "Protein kinase":', AC
        kinase.getparent().remove(kinase)
        continue
    # In cases where > 1 PK domain is found, add a warning to the DB entry. In some cases, a pseudokinase is present - these domains are not added.
    if len(PK_domains) > 1:
        if uniprot_kinases[k].findtext('name') == 'E2AK4_HUMAN':
            etree.SubElement(kinase,'warning').text = 'Kinase is annotated in UniProt wth both "Protein kinase 1" and "Protein kinase 2". "Protein kinase 1" is considered to be a pseudokinase domain. "Protein kinase 2" is considered active. Only the active PK domain is included in this DB.'
            PK_domains.pop(0)
        elif uniprot_kinases[k].findtext('name') in ['JAK1_HUMAN','JAK2_HUMAN','JAK3_HUMAN']:
            etree.SubElement(kinase,'warning').text = 'Kinase is annotated in UniProt wth both "Protein kinase 1" and "Protein kinase 2". Janus (Jak) tyrosine kinases (JAK1, JAK2 and JAK3) each contain a tyrosine kinase domain adjacent to a catalytically inactive pseudokinase domain. The pseudokinase domain interacts with and negatively regulates the active domain. The pseudokinase domain is the first one in the sequence. Only the active PK domain is included in this DB.'
            PK_domains.pop(0)
        elif uniprot_kinases[k].findtext('name') in ['KS6A1_HUMAN','KS6A2_HUMAN','KS6A3_HUMAN','KS6A4_HUMAN','KS6A5_HUMAN','KS6A6_HUMAN']:
            etree.SubElement(kinase,'warning').text = 'Kinase is annotated in UniProt wth both "Protein kinase 1" and "Protein kinase 2". Upon extracellular signal or mitogen stimulation, phosphorylated at Thr-573 in the C-terminal kinase domain (CTKD) by MAPK1/ERK2 and MAPK3/ERK1. The activated CTKD then autophosphorylates Ser-380, allowing binding of PDPK1, which in turn phosphorylates Ser-221 in the N-terminal kinase domain (NTKD) leading to the full activation of the protein and subsequent phosphorylation of the substrates by the NTKD. Both PK domains are included in this DB.'
        elif uniprot_kinases[k].findtext('name') == 'KS6C1_HUMAN':
            etree.SubElement(kinase,'warning').text = 'Kinase is annotated in UniProt wth both "Protein kinase 1" and "Protein kinase 2". The first protein kinase domain appears to be a pseudokinase domain as it does not contain the classical characteristics, such as the ATP-binding motif, ATP-binding site and active site. Only "Protein kinase 2" is included in this DB.'
            PK_domains.pop(0)
        elif uniprot_kinases[k].findtext('name') == 'OBSCN_HUMAN':
            etree.SubElement(kinase,'warning').text = 'Kinase is annotated in UniProt wth both "Protein kinase 1" and "Protein kinase 2". Neither are described as pseudokinases, although are not specifically described as catalytically active either. Both PK domains are included in this DB.'
        elif uniprot_kinases[k].findtext('name') == 'SPEG_HUMAN':
            etree.SubElement(kinase,'warning').text = 'Kinase is annotated in UniProt wth both "Protein kinase 1" and "Protein kinase 2". Neither are described as pseudokinases. Both PK domains are included in this DB.'
        elif uniprot_kinases[k].findtext('name') == 'TAF1_HUMAN':
            etree.SubElement(kinase,'warning').text = 'Kinase is annotated in UniProt wth both "Protein kinase 1" and "Protein kinase 2". Neither are described as pseudokinases. Both PK domains are included in this DB.'
        elif uniprot_kinases[k].findtext('name') == 'TYK2_HUMAN':
            etree.SubElement(kinase,'warning').text = 'Kinase is annotated in UniProt wth both "Protein kinase 1" and "Protein kinase 2". Neither are described as pseudokinases. Both PK domains are included in this DB.'
        else:
            raise Exception, 'More than 1 domain found containing "Protein kinase". Please check the following kinase and adjust the script: %s' % entry_name
    # And a couple of cases with one PK domain which are considered inactive. These kinase entries are removed completely.
    if PK_domains[0].attrib['description'] == 'Protein kinase; truncated':
        # PLK5_HUMAN. Kinase considered inactive. Protein kinase domain is truncated. Remove it.
        print 'Removing kinase as PK domain is truncated and considered inactive:', AC
        kinase.getparent().remove(kinase)
        continue
    elif PK_domains[0].attrib['description'] == 'Protein kinase; inactive':
        # PTK7_HUMAN. Kinase considered inactive. Remove it.
        print 'Removing kinase as PK domain is considered inactive:', AC
        kinase.getparent().remove(kinase)
        continue
    # Finally, add the PK domains to the new database
    for x_iter,x in enumerate(PK_domains):
        # First calculate the PK domain length and sequence
        pk_description = x.get('description')
        pk_begin = int( x.find('./location/begin').attrib['position'] )
        pk_end = int( x.find('./location/end').attrib['position'] )
        pk_length = pk_end - pk_begin + 1
        #PK_domain = deepcopy(x)
        PK_domain = etree.Element('pk_domain')
        PK_domain.set('description', pk_description)
        PK_domain.set('begin', str(pk_begin))
        PK_domain.set('end', str(pk_end))
        PK_domain.set('length', str(pk_length))
        PK_domain.set('id', str(x_iter))
        PK_domain.set('kinDB_id', (entry_name + '_' + AC + '_PK' + str(x_iter)))

        #location = PK_domain.find('./location')
        #etree.SubElement(location, 'length').text = str(pk_length)
        domain_sequence = seqwrap(sequence[pk_begin-1:pk_end])
        etree.SubElement(PK_domain, 'sequence').text = '\n' + domain_sequence
        kinase_uniprot.append(PK_domain)

    # = References to other DBs =
    # NCBI Gene
    GeneIDs = [x.get('id') for x in uniprot_kinases[k].findall('./dbReference[@type="GeneID"]')]
    # XXX: exceptions for kinases which have no GeneIDs annotated; LMTK3 RefSeq status is PROVISIONAL; RIPK4 presumably RefSeq sequence is not an exact match; SIK3 RefSeq status is VALIDATED
    # Will add these manually, since we are mainly using GeneID to collect publications currently
    if entry_name == 'LMTK3_HUMAN':
        GeneIDs = ['114783']
    if entry_name == 'RIPK4_HUMAN':
        GeneIDs = ['54101']
    if entry_name == 'SIK3_HUMAN':
        GeneIDs = ['23387']
    if len(GeneIDs) > 0:
        NCBI_Gene_node = etree.SubElement(kinase, 'NCBI_Gene')
    for GeneID in GeneIDs:
        # XXX: exceptions for SGK3_HUMAN and TNI3K_HUMAN, which have two GeneIDs annotated; in each case, one is a readthrough fusion protein - ignore these GeneIDs
        if GeneID in ['100533105', '100526835']:
            continue
        NCBI_Gene_entry_node = etree.SubElement(NCBI_Gene_node, 'entry')
        NCBI_Gene_entry_node.set('ID', GeneID)

    EnsemblGeneIDs = uniprot_kinases[k].findall('./dbReference[@type="Ensembl"]/property[@type="gene ID"]')
    EnsemblGeneIDs_set = set( [ id.attrib['value'] for id in EnsemblGeneIDs ] )
    for EnsemblGeneID in EnsemblGeneIDs_set:
        EnsemblGeneID_element = etree.Element('EnsemblGeneID')
        EnsemblGeneID_element.text = EnsemblGeneID
        kinase.append(EnsemblGeneID_element)

    HGNC_dbRefs = uniprot_kinases[k].findall('./dbReference[@type="HGNC"]')
    if len(HGNC_dbRefs) > 0:
        HGNC_element = etree.SubElement(kinase, 'HGNC')
        for HGNC_dbRef in HGNC_dbRefs:
            ID = HGNC_dbRef.get('id')
            Approved_Symbol = HGNC_dbRef.find('property[@type="gene designation"]').get('value')
            HGNC_entry_element = etree.SubElement(HGNC_element, 'entry')
            HGNC_entry_element.set('ID', ID)
            HGNC_entry_element.set('Approved_Symbol', Approved_Symbol)

    # = Family information =
    similarity_comments = uniprot_kinases[k].xpath('./comment[@type="similarity"]')
    family_found = False
    for s in similarity_comments:
        for f in kinase_family_uniprot_similarity_text.keys():
            if f in s.findtext('text'):
                kinase_uniprot.set('family', kinase_family_uniprot_similarity_text[f])
                family_found = True
    if family_found == False:
        kinase_uniprot.set('family', '')

    # = PDB entries (from UniProt XML) =
    pdbs = uniprot_kinases[k].findall('./dbReference[@type="PDB"]')
    kinase_PK_domains = kinase_uniprot.findall('pk_domain')
    for p in pdbs:
        # Only keep XRC structures (no NMR or Model)
        if p.find('property[@type="method"]') == None:
            if p.attrib['id'] == '2LV6':
                continue  # 2LV6 has no method listed - it is actually an NMR structure, including only a very short fragment of the kinase, outside the PK domain
        elif p.find('property[@type="method"]').attrib['value'] == 'X-ray':
            pk_pdb = etree.Element('pk_pdb')
            pdbid = p.attrib['id']
            if pdbid in ignore_uniprot_pdbs:
                continue
            pk_pdb.set('id', pdbid)
            resolution = p.find('property[@type="resolution"]').attrib['value']
            pk_pdb.set('resolution', resolution)
            chains_span_str = p.find('property[@type="chains"]').attrib['value']
            chains_span = parse_uniprot_pdbref_chains(chains_span_str)
            chains_added = 0
            for c in chains_span.keys():
                chains = etree.Element('chain')
                chains.set('id', c)
                pdb_begin = chains_span[c][0]
                pdb_end = chains_span[c][1]
                # Use the begin and end info to decide if this pdb chain includes the pk_domain. But we will get other sequence info from sifts XML files, using gather-protein_databank.py
                # Have to check against each PK domain
                for x_iter,x in enumerate(kinase_PK_domains):
                    pk_begin = int(x.get('begin'))
                    pk_end = int(x.get('end'))
                    if (pdb_begin < pk_begin+30) & (pdb_end > pk_end-30):
                        chains.set('pk_domain_id', str(x_iter))
                        chains.set('begin', str(pdb_begin))
                        chains.set('end', str(pdb_end))
                        pk_pdb.append(chains)
                        chains_added += 1
                    else:
                        continue

                if chains_added > 0:
                    kinase.append(pk_pdb)

# write the XML DB
database_out_filepath = open(database_out_filepath , 'w')
database_out_filepath.write( etree.tostring(root, pretty_print=True) )
database_out_filepath.close()

