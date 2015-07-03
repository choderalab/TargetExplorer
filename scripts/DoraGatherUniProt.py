#!/usr/bin/env python
#==============================================================================
# Imports
#==============================================================================

import argparse
from targetexplorer.core import read_project_config
from targetexplorer.uniprot_tmp import GatherUniProt


def parse_arguments():
        argparser = argparse.ArgumentParser(description='Gather UniProt')
        argparser.add_argument('--use_existing_data', help='Do not download a new UniProt document. Only works if an existing document is present.', action='store_true', default=False)
        argparser.add_argument('--count_nonselected_domain_names', help='Count the number of occurrences of domain names which are not selected by the regex and write to the "selected_domain_names.txt" file (may take a little while)', action='store_true', default=False)
        return argparser.parse_args()

args = parse_arguments()

project_config = read_project_config()

GatherUniProt(
    use_existing_data=args.use_existing_data,
    count_nonselected_domain_names=args.count_nonselected_domain_names,
    uniprot_query=project_config['uniprot_query'],
    uniprot_domain_regex=project_config['uniprot_domain_regex'],
    ignore_uniprot_pdbs=project_config['ignore_uniprot_pdbs'],
)

print('Done.')






















#
# project_config = read_project_config()
# uniprot_query = project_config['uniprot_query']
# uniprot_domain_regex = project_config['uniprot_domain_regex']
#
# #==============================================================================
# # Parameters
# #==============================================================================
#
# external_data_dir = 'external-data'
# uniprot_data_dir = os.path.join(external_data_dir, 'UniProt')
#
# if not os.path.exists(uniprot_data_dir):
#     os.mkdir(uniprot_data_dir)
#
# uniprot_xml_out_filepath = os.path.join(uniprot_data_dir, 'uniprot-search.xml')
#
# domain_names_filename = 'selected_domain_names.txt'
#
# # 1GQ5 is referenced by kinase P16234. The kinase is not in the actual structure.
# ignore_uniprot_pdbs = ['1GQ5']
#
# argparser = argparse.ArgumentParser(description='Gather UniProt')
# argparser.add_argument('--use_existing_data', help='Do not download a new UniProt document. Only works if an existing document is present.', action='store_true', default=False)
# argparser.add_argument('--count_nonselected_domain_names', help='Count the number of occurrences of domain names which are not selected by the regex and write to the "selected_domain_names.txt" file (may take a little while)', action='store_true', default=False)
# args = argparser.parse_args()
#
# now = datetime.datetime.utcnow()
#
# parser = etree.XMLParser(remove_blank_text=True, huge_tree=True)
#
# # get current crawl number
# crawldata_row = models.CrawlData.query.first()
# current_crawl_number = crawldata_row.current_crawl_number
# print 'Current crawl number: %d' %  current_crawl_number
#
# #==============================================================================
# # Retrieve data from uniprot and store to local file
# #==============================================================================
#
# # Unless args.use_existing_data is set to true, retrieve a new set of data from UniProt
# if os.path.exists(uniprot_xml_out_filepath) and args.use_existing_data:
#     print 'UniProt XML document found at:', uniprot_xml_out_filepath
# else:
#     print 'Retrieving new XML document from UniProt website.'
#     xml_text = targetexplorer.UniProt.retrieve_uniprot(uniprot_query)
#     if len(xml_text) == 0:
#         raise Exception('UniProt search returned no entries.')
#     print 'Saving new XML document as:', uniprot_xml_out_filepath
#     with open(uniprot_xml_out_filepath, 'w') as uniprot_xml_file:
#         uniprot_xml_file.write(xml_text + '\n')
#
# # Read in the UniProt XML document
# print 'Reading UniProt XML document:', uniprot_xml_out_filepath
# uniprot_xml = etree.parse(uniprot_xml_out_filepath, parser).getroot()
#
# domain_names_str = 'Regex: %s\n' % uniprot_domain_regex
# uniprot_entries = uniprot_xml.findall('entry')
# nuniprot_entries = len(uniprot_entries)
# # Note that xpath querying is case-sensitive
# domain_names_str += 'Number of entries in UniProt XML document: %d\n' % nuniprot_entries
# all_domains = uniprot_xml.findall('./entry/feature[@type="domain"]')
# domain_names_str += 'Total number of domains: %d\n' % len(all_domains)
#
# selected_domains = uniprot_xml.xpath('entry/feature[@type="domain"][match_regex(@description, "%s")]' % uniprot_domain_regex, extensions = { (None, 'match_regex'): targetexplorer.core.xpath_match_regex_case_sensitive })
# domain_names_str += 'Number of domains matching regex: %d\n\n' % len(selected_domains)
#
# domain_names_str += '= Unique domain names which match regex =\n'
# selected_domain_names = list(set([ d.get('description') for d in selected_domains ]))
# selected_domain_name_counts = [ len( uniprot_xml.findall('entry/feature[@type="domain"][@description="%s"]' % name) ) for name in selected_domain_names ]
# for i in range(len(selected_domain_names)):
#     domain_names_str += '{:^{name_width}s} : {:>{pop_width}d}\n'.format(selected_domain_names[i],
#         selected_domain_name_counts[i],
#         name_width=max([len(n)+4 for n in selected_domain_names]),
#         pop_width=max([len(str(p))+1 for p in selected_domain_name_counts])
#     )
# domain_names_str += '\n'
# print domain_names_str,
#
# print '(Unique domain names which do not match regex will be output to %s)' % domain_names_filename
#
# domain_names_str += '= Unique domain names which do not match regex =\n'
# nonselected_domain_names = list(set([ d.get('description') for d in all_domains if d.get('description') not in selected_domain_names ]))
#
# if args.count_nonselected_domain_names:
#     nonselected_domain_name_counts = [ int( uniprot_xml.xpath('count(entry/feature[@type="domain"][@description="%s"])' % name) ) for name in nonselected_domain_names ]
#     for i in range(len(nonselected_domain_names)):
#         domain_names_str += '{:^{name_width}s} : {:>{pop_width}d}\n'.format(nonselected_domain_names[i],
#             nonselected_domain_name_counts[i],
#             name_width=max([len(n)+4 for n in nonselected_domain_names]),
#             pop_width=max([len(str(p))+1 for p in nonselected_domain_name_counts]),
#         )
# else:
#     for i in range(len(nonselected_domain_names)):
#         domain_names_str += '{:^{name_width}s}\n'.format(nonselected_domain_names[i],
#             name_width=max([len(n)+4 for n in nonselected_domain_names]),
#         )
#
# domain_names_str += '\n'
# with open(domain_names_filename, 'w') as domain_names_file:
#     domain_names_file.write(domain_names_str)
#
#
# # ========
# # Iterate through each kinase from the UniProt XML document
# # ========
#
# for k in range(nuniprot_entries):
#
#     # = IDs and names =
#     ac = uniprot_entries[k].findtext('./accession')
#     entry_name = uniprot_entries[k].findtext('./name')
#     recommended_name = uniprot_entries[k].findtext('./protein/recommendedName/fullName')
#     gene_name_nodes = uniprot_entries[k].findall('./gene/name')
#     gene_name_data = []
#     for gene_name_node in gene_name_nodes:
#         gene_name = gene_name_node.text
#         gene_name_type = gene_name_node.get('type')
#         gene_name_obj = models.UniProtGeneName(crawl_number=current_crawl_number, gene_name=gene_name, gene_name_type=gene_name_type)
#         gene_name_data.append(gene_name_obj)
#
#     # = Date entry was last modified in UniProt =
#     last_uniprot_update = uniprot_entries[k].get('modified')
#
#     # = Taxonomy =
#     uniprot_organism_node = uniprot_entries[k].find('organism')
#     ncbi_taxonid = uniprot_organism_node.find('dbReference[@type="NCBI Taxonomy"]').get('id')
#     taxon_name_scientific = uniprot_organism_node.findtext('name[@type="scientific"]')
#     taxon_name_common = uniprot_organism_node.findtext('name[@type="common"]')
#     lineage = uniprot_organism_node.find('lineage')
#     lineage_csv = ','.join([taxon.text for taxon in lineage.getchildren()])
#
#     # = Functions, disease associations, subcellular locations =
#     functions = []
#     disease_associations = []
#     subcellular_locations = []
#     for x in uniprot_entries[k].findall('./comment[@type="function"]'):
#         functions.append( models.UniProtFunction(crawl_number=current_crawl_number, function=x.findtext('./text')) )
#     for x in uniprot_entries[k].findall('./comment[@type="disease"]'):
#         disease_associations.append( models.UniProtDiseaseAssociation(crawl_number=current_crawl_number, disease_association=x.findtext('./text')) )
#     for x in uniprot_entries[k].findall('./comment[@type="subcellular location"]'):
#         subcellular_locations.append( models.UniProtSubcellularLocation(crawl_number=current_crawl_number, subcellular_location=x.findtext('./subcellularLocation/location')) )
#
#     # = Canonical isoform =
#
#     isoforms = []
#
#     # Returned UniProt XML contains sequence data only for the canonical isoform
#     uniprot_canonical_sequence_node = uniprot_entries[k].find('./sequence[@length][@mass]')
#     canonical_sequence = ''.join(uniprot_canonical_sequence_node.text.split())
#     canseq_length = uniprot_canonical_sequence_node.get('length')
#     canseq_mass = uniprot_canonical_sequence_node.get('mass')
#     canseq_date_modified = uniprot_canonical_sequence_node.get('modified')
#     canseq_version = uniprot_canonical_sequence_node.get('version')
#     uniprotisoform = models.UniProtIsoform(crawl_number=current_crawl_number, ac=ac+'-1', canonical=True, length=canseq_length, mass=canseq_mass, date_modified=canseq_date_modified, version=canseq_version, sequence=canonical_sequence)
#     isoforms.append((uniprotisoform, [])) # empty list for notes (which do not exist for the canonical sequence)
#
#     # = Alternative isoforms =
#     # Canonical isoform is given the attrib type="displayed", meaning that the sequence is displayed in the HTML version of the entry
#     # Example alt isoform:
#     #     <comment>
#     #         <isoform>
#     #             <id>P00519-2</id>
#     #             <name>IB</name>
#     #             <sequence type="described" ref="VSP_004957"/>
#     #             <note>Contains a N-myristoyl glycine at position 2.</note>
#     #         </isoform>
#     #     </comment>
#
#     for uniprot_isoform_node in uniprot_entries[k].findall('comment/isoform'):
#         isoform_ac = uniprot_isoform_node.findtext('id')
#         seq_node = uniprot_isoform_node.find('sequence')
#         notes = [models.UniProtIsoformNote(crawl_number=current_crawl_number, note=node.text) for node in uniprot_isoform_node.findall('note')]
#         if seq_node.get('type') != 'displayed':
#             uniprotisoform = models.UniProtIsoform(crawl_number=current_crawl_number, ac=isoform_ac, canonical=False)
#
#         isoforms.append((uniprotisoform, notes))
#
#     # = UniProt "Protein kinase" domain annotations =
#     # XXX TODO Generalize
#
#     if uniprot_domain_regex != None:
#         selected_domains = uniprot_entries[k].xpath('feature[@type="domain"][match_regex(@description, "%s")]' % uniprot_domain_regex, extensions = { (None, 'match_regex'): targetexplorer.core.xpath_match_regex_case_sensitive })
#     else:
#         selected_domains = uniprot_entries[k].findall('feature[@type="domain"]')
#
#
#
#
#     # = Exceptions =
#
#     # Skip if no matcing domains found
#     if len(selected_domains) < 1:
#         continue
#
#     # XXX exception for SG196_HUMAN, which does not have protein kinase activity, and acts as a mannose kinase instead
#     if entry_name == 'SG196_HUMAN':
#         print 'Skipping kinase as it does not have protein kinase activity (instead acts as a mannose kinase):', ac
#         continue
#
#     # In cases where > 1 PK domain is found, add a warning to the DB entry. In some cases, a pseudokinase is present - these domains are not added.
#     warnings_node = etree.Element('warnings')
#     if len(selected_domains) > 1:
#         if uniprot_entries[k].findtext('name') == 'E2AK4_HUMAN':
#             etree.SubElement(warnings_node,'warning').text = 'Kinase is annotated in UniProt wth both "Protein kinase 1" and "Protein kinase 2". "Protein kinase 1" is considered to be a pseudokinase domain. "Protein kinase 2" is considered active. Only the active PK domain is included in this DB.'
#             selected_domains.pop(0)
#         elif uniprot_entries[k].findtext('name') in ['JAK1_HUMAN','JAK2_HUMAN','JAK3_HUMAN']:
#             etree.SubElement(warnings_node,'warning').text = 'Kinase is annotated in UniProt wth both "Protein kinase 1" and "Protein kinase 2". Janus (Jak) tyrosine kinases (JAK1, JAK2 and JAK3) each contain a tyrosine kinase domain adjacent to a catalytically inactive pseudokinase domain. The pseudokinase domain interacts with and negatively regulates the active domain. The pseudokinase domain is the first one in the sequence. Only the active PK domain is included in this DB.'
#             selected_domains.pop(0)
#         elif uniprot_entries[k].findtext('name') in ['KS6A1_HUMAN','KS6A2_HUMAN','KS6A3_HUMAN','KS6A4_HUMAN','KS6A5_HUMAN','KS6A6_HUMAN']:
#             etree.SubElement(warnings_node,'warning').text = 'Kinase is annotated in UniProt wth both "Protein kinase 1" and "Protein kinase 2". Upon extracellular signal or mitogen stimulation, phosphorylated at Thr-573 in the C-terminal kinase domain (CTKD) by MAPK1/ERK2 and MAPK3/ERK1. The activated CTKD then autophosphorylates Ser-380, allowing binding of PDPK1, which in turn phosphorylates Ser-221 in the N-terminal kinase domain (NTKD) leading to the full activation of the protein and subsequent phosphorylation of the substrates by the NTKD. Both PK domains are included in this DB.'
#         elif uniprot_entries[k].findtext('name') == 'KS6C1_HUMAN':
#             etree.SubElement(warnings_node,'warning').text = 'Kinase is annotated in UniProt wth both "Protein kinase 1" and "Protein kinase 2". The first protein kinase domain appears to be a pseudokinase domain as it does not contain the classical characteristics, such as the ATP-binding motif, ATP-binding site and active site. Only "Protein kinase 2" is included in this DB.'
#             selected_domains.pop(0)
#         elif uniprot_entries[k].findtext('name') == 'OBSCN_HUMAN':
#             etree.SubElement(warnings_node,'warning').text = 'Kinase is annotated in UniProt wth both "Protein kinase 1" and "Protein kinase 2". Neither are described as pseudokinases, although are not specifically described as catalytically active either. Both PK domains are included in this DB.'
#         elif uniprot_entries[k].findtext('name') == 'SPEG_HUMAN':
#             etree.SubElement(warnings_node,'warning').text = 'Kinase is annotated in UniProt wth both "Protein kinase 1" and "Protein kinase 2". Neither are described as pseudokinases. Both PK domains are included in this DB.'
#         elif uniprot_entries[k].findtext('name') == 'TAF1_HUMAN':
#             etree.SubElement(warnings_node,'warning').text = 'Kinase is annotated in UniProt wth both "Protein kinase 1" and "Protein kinase 2". Neither are described as pseudokinases. Both PK domains are included in this DB.'
#         elif uniprot_entries[k].findtext('name') == 'TYK2_HUMAN':
#             etree.SubElement(warnings_node,'warning').text = 'Kinase is annotated in UniProt wth both "Protein kinase 1" and "Protein kinase 2". Neither are described as pseudokinases. Both PK domains are included in this DB.'
#         else:
#             etree.SubElement(warnings_node,'warning').text = 'Kinase contains > 1 "Protein kinase*" domain. Not checked manually yet.'
#             #raise Exception, 'More than 1 domain found containing "Protein kinase". Please check the following kinase and adjust the script: %s' % entry_name
#     # And a couple of cases with one PK domain which are considered inactive. These kinase entries are removed completely.
#     if selected_domains[0].get('description') == 'Protein kinase; truncated':
#         # PLK5_HUMAN. Kinase considered inactive. Protein kinase domain is truncated. Remove it.
#         print 'Skipping kinase as PK domain is truncated and considered inactive:', ac
#         continue
#     elif selected_domains[0].get('description') == 'Protein kinase; inactive':
#         # PTK7_HUMAN. Kinase considered inactive. Remove it.
#         print 'Skipping kinase as PK domain is considered inactive:', ac
#         continue
#
#     # Finally, add the domains to the new database
#     domains_data = []
#     for x_iter,x in enumerate(selected_domains):
#         # First calculate the PK domain length and sequence
#         description = x.get('description')
#         begin = int( x.find('./location/begin').get('position') )
#         end = int( x.find('./location/end').get('position') )
#         length = end - begin + 1
#         targetid = entry_name + '_D' + str(x_iter)
#         domain_seq = canonical_sequence[begin-1:end]
#
#         domain_obj = models.UniProtDomain(crawl_number=current_crawl_number, targetid=targetid, description=description, begin=begin, end=end, length=length, sequence=domain_seq)
#         domains_data.append(domain_obj)
#
#     # = References to other DBs =
#     # NCBI Gene
#     ncbi_gene_entries = []
#     GeneIDs = [x.get('id') for x in uniprot_entries[k].findall('./dbReference[@type="GeneID"]')]
#     # XXX: exceptions for kinases which have no GeneIDs annotated; LMTK3 RefSeq status is PROVISIONAL; RIPK4 presumably RefSeq sequence is not an exact match; SIK3 RefSeq status is VALIDATED
#     # Will add these manually, since we are mainly using GeneID to collect publications currently
#     if entry_name == 'LMTK3_HUMAN':
#         GeneIDs = ['114783']
#     if entry_name == 'RIPK4_HUMAN':
#         GeneIDs = ['54101']
#     if entry_name == 'SIK3_HUMAN':
#         GeneIDs = ['23387']
#     for GeneID in GeneIDs:
#         # XXX: exceptions for SGK3_HUMAN and TNI3K_HUMAN, which have two GeneIDs annotated; in each case, one is a readthrough fusion protein - ignore these GeneIDs
#         if GeneID in ['100533105', '100526835']:
#             continue
#         ncbi_gene_entries.append( models.NCBIGeneEntry(crawl_number=current_crawl_number, gene_id=GeneID) )
#
#     # Ensembl
#
#     # transcript_data = {
#     #     'ENSMUST00000003710':
#     #         {
#     #             'gene':
#     #                 'ENSG000...',
#     #             'protein':
#     #                 'ENSP000...',
#     #         }
#     # }
#
#     ensembl_transcript_nodes = uniprot_entries[k].findall('./dbReference[@type="Ensembl"]')
#
#     ensembl_data = {}
#     for transcript_node in ensembl_transcript_nodes:
#         ensembl_transcript_id = transcript_node.get('id')
#
#         ensembl_gene_nodes = transcript_node.findall('property[@type="gene ID"]')
#         if len(ensembl_gene_nodes) > 1:
#             print('WARNING: Ensembl transcript {0} linked with > 1 gene ID'.format(ensembl_transcript_id))
#         ensembl_gene_id = ensembl_gene_nodes[0].get('value')
#
#         ensembl_protein_nodes = transcript_node.findall('property[@type="protein sequence ID"]')
#         if len(ensembl_protein_nodes) > 1:
#             print('WARNING: Ensembl transcript {0} linked with > 1 protein ID'.format(ensembl_transcript_id))
#         ensembl_protein_id = ensembl_protein_nodes[0].get('value')
#
#         uniprot_isoform_molecule_node = transcript_node.find('molecule')
#         if uniprot_isoform_molecule_node is not None:
#             uniprot_isoform_ac = uniprot_isoform_molecule_node.get('id')
#         else:
#             uniprot_isoform_ac = None
#
#         ensembl_data[ensembl_transcript_id] = {
#             'gene': ensembl_gene_id,
#             'protein': ensembl_protein_id,
#             'uniprot_isoform_ac': uniprot_isoform_ac
#         }
#
#
#     # HGNC
#     hgnc_entries = []
#     hgnc_dbrefs = uniprot_entries[k].findall('./dbReference[@type="HGNC"]')
#     for hgnc_dbref in hgnc_dbrefs:
#         hgnc_gene_id = hgnc_dbref.get('id')
#         approved_symbol = hgnc_dbref.find('property[@type="gene designation"]').get('value')
#         hgnc_entries.append(models.HGNCEntry(crawl_number=current_crawl_number, gene_id=hgnc_gene_id, approved_symbol=approved_symbol))
#
#     # = Family information =
#     similarity_comments = uniprot_entries[k].xpath('./comment[@type="similarity"]')
#     family = False
#     for s in similarity_comments:
#         for f in targetexplorer.UniProt.kinase_family_uniprot_similarity_text.keys():
#             if f in s.findtext('text'):
#                 family = targetexplorer.UniProt.kinase_family_uniprot_similarity_text[f]
#
#     # = PDB entries (from UniProt XML) =
#     # keep X-ray and NMR structures (not "Model")
#     pdbs = uniprot_entries[k].xpath('./dbReference[@type="PDB"]/property[@type="method"][@value="X-ray" or @value="NMR"]/..')
#     pdb_data = []
#     for p in pdbs:
#         pdbid = p.get('id')
#         if pdbid in ignore_uniprot_pdbs:
#             continue
#         pdb_method = p.find('property[@type="method"]').get('value')
#         resolution_node = p.find('property[@type="resolution"]')
#         resolution = resolution_node.get('value') if resolution_node != None else None
#         chains_span_str = p.find('property[@type="chains"]').get('value')
#         chains_span = targetexplorer.UniProt.parse_uniprot_pdbref_chains(chains_span_str)
#         chain_objs = []
#         for c in chains_span.keys():
#             chain_id = c
#             pdb_begin = chains_span[c][0]
#             pdb_end = chains_span[c][1]
#             # Use the begin and end info to decide if this pdb chain includes the pk_domain. But we will get other sequence info from sifts XML files, using gather-pdb_tmp.py
#             # Have to check against each PK domain
#             for domain_id, domain in enumerate(domains_data):
#                 pk_begin = domain.begin
#                 pk_end = domain.end
#                 if (pdb_begin < pk_begin+30) & (pdb_end > pk_end-30):
#                     chain_obj = models.PDBChain(crawl_number=current_crawl_number, chain_id=chain_id, domain_id=domain_id, begin=pdb_begin, end=pdb_end)
#                     chain_objs.append(chain_obj)
#                 else:
#                     continue
#
#         if len(chain_objs) > 0:
#             pdb_obj = models.PDB(crawl_number=current_crawl_number, pdbid=pdbid, method=pdb_method, resolution=resolution)
#             pdb_data.append({'pdb_obj': pdb_obj, 'chain_objs': chain_objs})
#
#     # # = Add the warnings node last (only if it contains any warnings) = #
#     # TODO still need this?
#     # if len(warnings_node) > 0:
#     #     DBentry.append(warnings_node)
#
#
#
#     # ========
#     # Construct data objects and add to db
#     # ========
#
#     dbentry = models.DBEntry(
#         crawl_number=current_crawl_number,
#         npdbs=len(pdb_data),
#         ndomains=len(domains_data),
#         nisoforms=len(isoforms),
#         nfunctions=len(functions),
#         ndiseaseassociations=len(disease_associations),
#     )
#     db.session.add(dbentry)
#     uniprot = models.UniProt(
#         crawl_number=current_crawl_number,
#         ac=ac,
#         entry_name=entry_name,
#         last_uniprot_update=last_uniprot_update,
#         ncbi_taxonid=ncbi_taxonid,
#         dbentry=dbentry,
#         recommended_name=recommended_name,
#         taxon_name_scientific=taxon_name_scientific,
#         taxon_name_common=taxon_name_common,
#         lineage=lineage_csv,
#     )
#     if family:
#         uniprot.family = family
#     db.session.add(uniprot)
#     for function_obj in functions:
#         function_obj.dbentry = dbentry
#         function_obj.uniprot_entry = uniprot
#         db.session.add(function_obj)
#     for disease_association_obj in disease_associations:
#         disease_association_obj.dbentry = dbentry
#         disease_association_obj.uniprot_entry = uniprot
#         db.session.add(disease_association_obj)
#     for subcellular_location_obj in subcellular_locations:
#         subcellular_location_obj.dbentry = dbentry
#         subcellular_location_obj.uniprot_entry = uniprot
#         db.session.add(subcellular_location_obj)
#     for isoform_data in isoforms:
#         isoform_obj = isoform_data[0]
#         notes = isoform_data[1]
#         isoform_obj.dbentry = dbentry
#         isoform_obj.uniprot_entry = uniprot
#         db.session.add(isoform_obj)
#         for note_obj in notes:
#             note_obj.uniprotisoform = isoform_obj
#             db.session.add(note_obj)
#     for domain_obj in domains_data:
#         domain_obj.dbentry = dbentry
#         domain_obj.uniprot_entry = uniprot
#         db.session.add(domain_obj)
#     for pdb_data_dict in pdb_data:
#         pdb_obj = pdb_data_dict['pdb_obj']
#         chain_objs = pdb_data_dict['chain_objs']
#         pdb_obj.dbentry = dbentry
#         db.session.add(pdb_obj)
#         for chain_obj in chain_objs:
#             chain_obj.pdb = pdb_obj
#             db.session.add(chain_obj)
#     for gene_name_obj in gene_name_data:
#         gene_name_obj.dbentry = dbentry
#         db.session.add(gene_name_obj)
#     for NCBIGeneEntry in ncbi_gene_entries:
#         NCBIGeneEntry.dbentry = dbentry
#         db.session.add(NCBIGeneEntry)
#     for HGNCEntry in hgnc_entries:
#         HGNCEntry.dbentry = dbentry
#         db.session.add(HGNCEntry)
#     for ensembl_transcript_id in ensembl_data:
#         ensembl_gene_id = ensembl_data[ensembl_transcript_id]['gene']
#         ensembl_gene_row = models.EnsemblGene(
#             crawl_number=current_crawl_number,
#             gene_id=ensembl_gene_id,
#             dbentry=dbentry,
#         )
#         db.session.add(ensembl_gene_row)
#
#         ensembl_transcript_row = models.EnsemblTranscript(
#             crawl_number=current_crawl_number,
#             transcript_id=ensembl_transcript_id,
#             ensembl_gene=ensembl_gene_row,
#         )
#         ensembl_transcript_uniprot_isoform_ac = ensembl_data[ensembl_transcript_id]['uniprot_isoform_ac']
#         if ensembl_transcript_uniprot_isoform_ac is not None:
#             matching_uniprot_isoform_obj = [isoform[0] for isoform in isoforms if isoform[0].ac == ensembl_transcript_uniprot_isoform_ac]
#             if len(matching_uniprot_isoform_obj) != 0:
#                 ensembl_transcript_row.uniprotisoform = matching_uniprot_isoform_obj[0]
#         db.session.add(ensembl_transcript_row)
#
#         ensembl_protein_id = ensembl_data[ensembl_transcript_id]['protein']
#         ensembl_protein_row = models.EnsemblProtein(
#             crawl_number=current_crawl_number,
#             protein_id=ensembl_protein_id,
#             ensembl_gene=ensembl_gene_row,
#             ensembl_transcript=ensembl_transcript_row,
#         )
#         db.session.add(ensembl_protein_row)
#
# # update db UniProt datestamp
# current_crawl_datestamp_row = models.DateStamps.query.filter_by(crawl_number=current_crawl_number).first()
# current_crawl_datestamp_row.uniprot_datestamp = now
# db.session.commit()
# print 'Done.'
