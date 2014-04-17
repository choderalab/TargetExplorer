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

import sys, datetime, os, copy, yaml, argparse
import TargetExplorer
import config
from lxml import etree
from app import models, db

#==============================================================================
# PARAMETERS
#==============================================================================

# if '-stage' in sys.argv:
#     run_mode = 'stage'
# elif '-dev' in sys.argv:
#     run_mode = 'dev'
# else:
#     run_mode = 'nowrite'

# print 'Running in mode: %s' % run_mode

database_dir = 'database'
external_data_dir = 'external-data'
uniprot_data_dir = os.path.join(external_data_dir, 'UniProt')

if not os.path.exists(uniprot_data_dir):
    os.mkdir(uniprot_data_dir)

uniprot_xml_out_filepath = os.path.join(uniprot_data_dir, 'uniprot-search.xml')

# 1GQ5 is referenced by kinase P16234. The kinase is not in the actual structure.
ignore_uniprot_pdbs = ['1GQ5']

argparser = argparse.ArgumentParser(description='Gather UniProt')
argparser.add_argument('--forcedl')
args = argparser.parse_args()

now = datetime.datetime.utcnow()
# datestamp = now.strftime(TargetExplorer.DB.datestamp_format_string)

parser = etree.XMLParser(remove_blank_text=True, huge_tree=True)

# with open('config.yaml') as config_file:
#     config = yaml.load(config_file)

#==============================================================================
# RETRIEVE DATA FROM UNIPROT AND STORE TO LOCAL FILE
#==============================================================================

# If the UniProt external data does not already exist, download it
if not os.path.exists(uniprot_xml_out_filepath) or args.forcedl:
    print 'UniProt XML document not found.'
    print 'Retrieving new XML document from UniProt website.'
    new_xml_text = TargetExplorer.UniProt.retrieve_uniprot(config.uniprot_query_string_url)
    print 'Saving new XML document as:', uniprot_xml_out_filepath
    with open(uniprot_xml_out_filepath, 'w') as uniprot_xml_file:
        uniprot_xml_file.write(new_xml_text + '\n')
else:
    print 'UniProt XML document found at:', uniprot_xml_out_filepath

# Read in the UniProt XML document
print 'Reading UniProt XML document:', uniprot_xml_out_filepath
uniprot_xml = etree.parse(uniprot_xml_out_filepath, parser).getroot()



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






# ========
# Remove all existing data from db
# ========

print 'Deleting all existing content in db-stage'
print 'Deleting %d DBEntry rows...' % models.DBEntry.query.delete()
print 'Deleting %d UniProt rows...' % models.UniProt.query.delete()
print 'Deleting %d UniProtDomain rows...' % models.UniProtDomain.query.delete()
print 'Deleting %d PDB rows...' % models.PDB.query.delete()
print ''

# ========
# Iterate through each kinase from the UniProt XML document
# ========

for k in range(nuniprot_entries):

    # = IDs and names =
    ac = uniprot_entries[k].findtext('./accession')
    entry_name = uniprot_entries[k].findtext('./name')
    recommended_name = uniprot_entries[k].findtext('./protein/recommendedName/fullName')
    gene_name_nodes = uniprot_entries[k].findall('./gene/name')
    gene_name_data = []
    for gene_name_node in gene_name_nodes:
        gene_name = gene_name_node.text
        gene_name_type = gene_name_node.get('type')
        gene_name_obj = models.UniProtGeneName(gene_name=gene_name, gene_name_type=gene_name_type)
        gene_name_data.append(gene_name_obj)

    # = Date entry was last modified in UniProt =
    last_uniprot_update = uniprot_entries[k].attrib['modified']

    # XXX exception for SG196_HUMAN, which does not have protein kinase activity, and acts as a mannose kinase instead
    if entry_name == 'SG196_HUMAN':
        print 'Skipping kinase as it does not have protein kinase activity (instead acts as a mannose kinase):', ac
        continue

    # = Taxonomy =
    uniprot_organism_node = uniprot_entries[k].find('organism')
    ncbi_taxonid = uniprot_organism_node.find('dbReference[@type="NCBI Taxonomy"]').get('id')
    taxon_name_scientific = uniprot_organism_node.findtext('name[@type="scientific"]')
    taxon_name_common = uniprot_organism_node.findtext('name[@type="common"]')
    lineage = uniprot_organism_node.find('lineage')
    lineage_csv = ','.join([taxon.text for taxon in lineage.getchildren()])

    # = Functions, disease associations =
    functions = []
    disease_associations = []
    for x in uniprot_entries[k].findall('./comment[@type="function"]'):
        functions.append( models.UniProtFunction(function=x.findtext('./text')) )
    for x in uniprot_entries[k].findall('./comment[@type="disease"]'):
        disease_associations.append( models.UniProtDiseaseAssociation(disease_association=x.findtext('./text')) )

    # = Canonical isoform =

    isoforms = []

    # Returned UniProt XML contains sequence data only for the canonical isoform
    uniprot_canonical_sequence_node = uniprot_entries[k].find('./sequence[@length][@mass]')
    canonical_sequence = ''.join(uniprot_canonical_sequence_node.text.split())
    canseq_length = uniprot_canonical_sequence_node.get('length')
    canseq_mass = uniprot_canonical_sequence_node.get('mass')
    canseq_date_modified = uniprot_canonical_sequence_node.get('modified')
    canseq_version = uniprot_canonical_sequence_node.get('version')
    uniprotisoform = models.UniProtIsoform(ac=ac+'-1', canonical=True, length=canseq_length, mass=canseq_mass, date_modified=canseq_date_modified, version=canseq_version, sequence=canonical_sequence)
    isoforms.append((uniprotisoform, [])) # empty list for notes (which do not exist for the canonical sequence)

    # = Alternative isoforms =
    # Canonical isoform is given the attrib type="displayed", meaning that the sequence is displayed in the HTML version of the entry
    # Example alt isoform:
    #     <comment>
    #         <isoform>
    #             <id>P00519-2</id>
    #             <name>IB</name>
    #             <sequence type="described" ref="VSP_004957"/>
    #             <note>Contains a N-myristoyl glycine at position 2.</note>
    #         </isoform>
    #     </comment>

    for uniprot_isoform_node in uniprot_entries[k].findall('comment/isoform'):
        isoform_ac = uniprot_isoform_node.findtext('id')
        seq_node = uniprot_isoform_node.find('sequence')
        notes = [models.UniProtIsoformNote(note=node.text) for node in uniprot_isoform_node.findall('note')]
        if seq_node.get('type') != 'displayed':
            uniprotisoform = models.UniProtIsoform(ac=isoform_ac, canonical=False)

        isoforms.append((uniprotisoform, notes))

    # = UniProt "Protein kinase" domain annotations =
    # XXX TODO Generalize

    if config.uniprot_domain_regex != None:
        selected_domains = uniprot_entries[k].xpath('feature[@type="domain"][match_regex(@description, "%s")]' % config.uniprot_domain_regex, extensions = { (None, 'match_regex'): TargetExplorer.core.xpath_match_regex_case_sensitive })
    else:
        selected_domains = uniprot_entries[k].findall('feature[@type="domain"]')




    # XXX exceptions
    # These are the entries for which "Protein kinase" domains are known to be not found (case sensitive):
    # kinases_with_no_PK_domain = ['ALPK1_HUMAN', 'ALPK2_HUMAN', 'ALPK3_HUMAN', 'EF2K_HUMAN', 'TRPM6_HUMAN', 'TRPM7_HUMAN']
    # These are all alpha-kinases, which have no identity with typical protein kinases.
    # These kinases will therefore be deleted from the database.
    if len(selected_domains) < 1:
        print 'Skipping kinase as it does not possess a domain annotation containing "Protein kinase":', ac
        continue
    # In cases where > 1 PK domain is found, add a warning to the DB entry. In some cases, a pseudokinase is present - these domains are not added.
    warnings_node = etree.Element('warnings')
    if len(selected_domains) > 1:
        if uniprot_entries[k].findtext('name') == 'E2AK4_HUMAN':
            etree.SubElement(warnings_node,'warning').text = 'Kinase is annotated in UniProt wth both "Protein kinase 1" and "Protein kinase 2". "Protein kinase 1" is considered to be a pseudokinase domain. "Protein kinase 2" is considered active. Only the active PK domain is included in this DB.'
            selected_domains.pop(0)
        elif uniprot_entries[k].findtext('name') in ['JAK1_HUMAN','JAK2_HUMAN','JAK3_HUMAN']:
            etree.SubElement(warnings_node,'warning').text = 'Kinase is annotated in UniProt wth both "Protein kinase 1" and "Protein kinase 2". Janus (Jak) tyrosine kinases (JAK1, JAK2 and JAK3) each contain a tyrosine kinase domain adjacent to a catalytically inactive pseudokinase domain. The pseudokinase domain interacts with and negatively regulates the active domain. The pseudokinase domain is the first one in the sequence. Only the active PK domain is included in this DB.'
            selected_domains.pop(0)
        elif uniprot_entries[k].findtext('name') in ['KS6A1_HUMAN','KS6A2_HUMAN','KS6A3_HUMAN','KS6A4_HUMAN','KS6A5_HUMAN','KS6A6_HUMAN']:
            etree.SubElement(warnings_node,'warning').text = 'Kinase is annotated in UniProt wth both "Protein kinase 1" and "Protein kinase 2". Upon extracellular signal or mitogen stimulation, phosphorylated at Thr-573 in the C-terminal kinase domain (CTKD) by MAPK1/ERK2 and MAPK3/ERK1. The activated CTKD then autophosphorylates Ser-380, allowing binding of PDPK1, which in turn phosphorylates Ser-221 in the N-terminal kinase domain (NTKD) leading to the full activation of the protein and subsequent phosphorylation of the substrates by the NTKD. Both PK domains are included in this DB.'
        elif uniprot_entries[k].findtext('name') == 'KS6C1_HUMAN':
            etree.SubElement(warnings_node,'warning').text = 'Kinase is annotated in UniProt wth both "Protein kinase 1" and "Protein kinase 2". The first protein kinase domain appears to be a pseudokinase domain as it does not contain the classical characteristics, such as the ATP-binding motif, ATP-binding site and active site. Only "Protein kinase 2" is included in this DB.'
            selected_domains.pop(0)
        elif uniprot_entries[k].findtext('name') == 'OBSCN_HUMAN':
            etree.SubElement(warnings_node,'warning').text = 'Kinase is annotated in UniProt wth both "Protein kinase 1" and "Protein kinase 2". Neither are described as pseudokinases, although are not specifically described as catalytically active either. Both PK domains are included in this DB.'
        elif uniprot_entries[k].findtext('name') == 'SPEG_HUMAN':
            etree.SubElement(warnings_node,'warning').text = 'Kinase is annotated in UniProt wth both "Protein kinase 1" and "Protein kinase 2". Neither are described as pseudokinases. Both PK domains are included in this DB.'
        elif uniprot_entries[k].findtext('name') == 'TAF1_HUMAN':
            etree.SubElement(warnings_node,'warning').text = 'Kinase is annotated in UniProt wth both "Protein kinase 1" and "Protein kinase 2". Neither are described as pseudokinases. Both PK domains are included in this DB.'
        elif uniprot_entries[k].findtext('name') == 'TYK2_HUMAN':
            etree.SubElement(warnings_node,'warning').text = 'Kinase is annotated in UniProt wth both "Protein kinase 1" and "Protein kinase 2". Neither are described as pseudokinases. Both PK domains are included in this DB.'
        else:
            etree.SubElement(warnings_node,'warning').text = 'Kinase contains > 1 "Protein kinase*" domain. Not checked manually yet.'
            #raise Exception, 'More than 1 domain found containing "Protein kinase". Please check the following kinase and adjust the script: %s' % entry_name
    # And a couple of cases with one PK domain which are considered inactive. These kinase entries are removed completely.
    if selected_domains[0].attrib['description'] == 'Protein kinase; truncated':
        # PLK5_HUMAN. Kinase considered inactive. Protein kinase domain is truncated. Remove it.
        print 'Skipping kinase as PK domain is truncated and considered inactive:', ac
        continue
    elif selected_domains[0].attrib['description'] == 'Protein kinase; inactive':
        # PTK7_HUMAN. Kinase considered inactive. Remove it.
        print 'Skipping kinase as PK domain is considered inactive:', ac
        continue

    # Finally, add the domains to the new database
    domains_data = []
    for x_iter,x in enumerate(selected_domains):
        # First calculate the PK domain length and sequence
        description = x.get('description')
        begin = int( x.find('./location/begin').attrib['position'] )
        end = int( x.find('./location/end').attrib['position'] )
        length = end - begin + 1
        targetid = entry_name + '_D' + str(x_iter)
        domain_seq = canonical_sequence[begin-1:end]

        domain_obj = models.UniProtDomain(targetid=targetid, description=description, begin=begin, end=end, length=length, sequence=domain_seq)
        domains_data.append(domain_obj)

    # = References to other DBs =
    # NCBI Gene
    ncbi_gene_entries = []
    GeneIDs = [x.get('id') for x in uniprot_entries[k].findall('./dbReference[@type="GeneID"]')]
    # XXX: exceptions for kinases which have no GeneIDs annotated; LMTK3 RefSeq status is PROVISIONAL; RIPK4 presumably RefSeq sequence is not an exact match; SIK3 RefSeq status is VALIDATED
    # Will add these manually, since we are mainly using GeneID to collect publications currently
    if entry_name == 'LMTK3_HUMAN':
        GeneIDs = ['114783']
    if entry_name == 'RIPK4_HUMAN':
        GeneIDs = ['54101']
    if entry_name == 'SIK3_HUMAN':
        GeneIDs = ['23387']
    for GeneID in GeneIDs:
        # XXX: exceptions for SGK3_HUMAN and TNI3K_HUMAN, which have two GeneIDs annotated; in each case, one is a readthrough fusion protein - ignore these GeneIDs
        if GeneID in ['100533105', '100526835']:
            continue
        ncbi_gene_entries.append( models.NCBIGeneEntry(gene_id=GeneID) )

    # Ensembl
    ensembl_gene_entries = []
    ensembl_gene_ids = uniprot_entries[k].findall('./dbReference[@type="Ensembl"]/property[@type="gene ID"]')
    ensembl_gene_ids_set = set( [ id.attrib['value'] for id in ensembl_gene_ids ] )
    for ensembl_gene_id in ensembl_gene_ids_set:
        ensembl_gene_entries.append( models.EnsemblGeneEntry(gene_id=ensembl_gene_id) )

    # HGNC
    hgnc_entries = []
    hgnc_dbrefs = uniprot_entries[k].findall('./dbReference[@type="HGNC"]')
    for hgnc_dbref in hgnc_dbrefs:
        hgnc_gene_id = hgnc_dbref.get('id')
        approved_symbol = hgnc_dbref.find('property[@type="gene designation"]').get('value')
        hgnc_entries.append( models.HGNCEntry(gene_id=hgnc_gene_id, approved_symbol=approved_symbol) )

    # = Family information =
    similarity_comments = uniprot_entries[k].xpath('./comment[@type="similarity"]')
    family = False
    for s in similarity_comments:
        for f in TargetExplorer.UniProt.kinase_family_uniprot_similarity_text.keys():
            if f in s.findtext('text'):
                family = TargetExplorer.UniProt.kinase_family_uniprot_similarity_text[f]

    # = PDB entries (from UniProt XML) =
    pdbs = uniprot_entries[k].findall('./dbReference[@type="PDB"]')
    pdb_data = []
    for p in pdbs:
        # Only keep XRC structures (no NMR or Model) TODO should keep NMR models
        if p.find('property[@type="method"]') == None:
            if p.attrib['id'] == '2LV6':
                continue  # 2LV6 has no method listed - it is actually an NMR structure, including only a very short fragment of the kinase, outside the PK domain
        elif p.find('property[@type="method"]').attrib['value'] == 'X-ray':
            pdbid = p.attrib['id']
            if pdbid in ignore_uniprot_pdbs:
                continue
            resolution = p.find('property[@type="resolution"]').attrib['value']
            chains_span_str = p.find('property[@type="chains"]').attrib['value']
            chains_span = TargetExplorer.UniProt.parse_uniprot_pdbref_chains(chains_span_str)
            chains_added = 0
            for c in chains_span.keys():
                chainID = c
                pdb_begin = chains_span[c][0]
                pdb_end = chains_span[c][1]
                # Use the begin and end info to decide if this pdb chain includes the pk_domain. But we will get other sequence info from sifts XML files, using gather-pdb.py
                # Have to check against each PK domain
                for d,domain in enumerate(domains_data):
                    pk_begin = domain.begin
                    pk_end = domain.end
                    if (pdb_begin < pk_begin+30) & (pdb_end > pk_end-30):
                        domainID = str(d)
                        pdb_begin = str(pdb_begin)
                        pdb_end = str(pdb_end)
                        chains_added += 1
                    else:
                        continue

                if chains_added > 0:
                    pdb_obj = models.PDB(pdbid=pdbid)
                    pdb_data.append(pdb_obj)

    # # = Add the warnings node last (only if it contains any warnings) = #
    # TODO still need this?
    # if len(warnings_node) > 0:
    #     DBentry.append(warnings_node)



    # ========
    # Construct data objects and add to db
    # ========

    dbentry = models.DBEntry()
    db.session.add(dbentry)
    uniprot = models.UniProt(ac=ac, entry_name=entry_name, last_uniprot_update=last_uniprot_update, ncbi_taxonid=ncbi_taxonid, dbentry=dbentry, recommended_name=recommended_name, taxon_name_scientific=taxon_name_scientific, taxon_name_common=taxon_name_common, lineage=lineage_csv)
    if family:
        uniprot.family = family
    db.session.add(uniprot)
    for function_obj in functions:
        function_obj.dbentry = dbentry
        function_obj.uniprotentry = uniprot
        db.session.add(function_obj)
    for disease_association_obj in disease_associations:
        disease_association_obj.dbentry = dbentry
        disease_association_obj.uniprotentry = uniprot
        db.session.add(disease_association_obj)
    for isoform_data in isoforms:
        isoform_obj = isoform_data[0]
        notes = isoform_data[1]
        isoform_obj.dbentry = dbentry
        isoform_obj.uniprotentry = uniprot
        db.session.add(isoform_obj)
        for note_obj in notes:
            note_obj.uniprotisoform = isoform_obj
            db.session.add(note_obj)
    for domain_obj in domains_data:
        domain_obj.dbentry = dbentry
        domain_obj.uniprotentry = uniprot
        db.session.add(domain_obj)
    for pdb_obj in pdb_data:
        pdb_obj.dbentry = dbentry
        db.session.add(pdb_obj)
    for gene_name_obj in gene_name_data:
        gene_name_obj.dbentry = dbentry
        db.session.add(gene_name_obj)
    for NCBIGeneEntry in ncbi_gene_entries:
        NCBIGeneEntry.dbentry = dbentry
        db.session.add(NCBIGeneEntry)
    for EnsemblGeneEntry in ensembl_gene_entries:
        EnsemblGeneEntry.dbentry = dbentry
        db.session.add(EnsemblGeneEntry)
    for HGNCEntry in hgnc_entries:
        HGNCEntry.dbentry = dbentry
        db.session.add(HGNCEntry)

# update db UniProt datestamp
version_row = models.Version.query.all()[0]
version_row.uniprotdatestamp = now
db.session.commit()
print 'Done.'

