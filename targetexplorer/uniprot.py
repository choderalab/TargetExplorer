import os
import urllib
import urllib2
import datetime
import re
from lxml import etree
from targetexplorer.core import external_data_dirpath, xml_parser, logger
from targetexplorer.core import xpath_match_regex_case_sensitive, read_manual_overrides
from targetexplorer.flaskapp import db, models

# This dict converts family information listed in UniProt in the similarity comments to codes similar to those used in the kinase.com poster
# Note that a UniProt "family" is equivalent to a Manning et al. "group". Also, there are a few additional families annotated in UniProt.
kinase_family_uniprot_similarity_text = {
    'AGC Ser/Thr protein kinase family': 'AGC',
    'CAMK Ser/Thr protein kinase family': 'CAMK',
    'CMGC Ser/Thr protein kinase family': 'CMGC',
    'CK1 Ser/Thr protein kinase family': 'CK1',
    'STE Ser/Thr protein kinase family': 'STE',
    'TKL Ser/Thr protein kinase family': 'TKL',
    'Tyr protein kinase family': 'TK',
    'NEK Ser/Thr protein kinase family': 'NEK',
    'RIO-type Ser/Thr kinase family': 'RIO-type'
}


class GatherUniProt(object):
    def __init__(self,
                 uniprot_query=None,
                 uniprot_domain_regex=None,
                 use_existing_data=False,
                 count_nonselected_domain_names=False,
                 run_main=True,
                 commit_to_db=True
                 ):
        """
        Parameters
        ----------
        uniprot_query: str
            e.g. 'mnemonic:ABL1_HUMAN'
        uniprot_domain_regex: str
            e.g. '^Protein kinase(?!; truncated)(?!; inactive)'
        use_existing_data: bool
        count_nonselected_domain_names: bool
        run_main: bool
        """
        self.commit_to_db = commit_to_db
        self.uniprot_query = uniprot_query
        self.uniprot_domain_regex = uniprot_domain_regex
        self.use_existing_data = use_existing_data
        self.count_nonselected_domain_names = count_nonselected_domain_names
        if run_main:
            self.setup()
            self.setup_manual_overrides()
            self.get_uniprot_data()
            self.extract_uniprot_entries_and_domains()
            if self.uniprot_domain_regex:
                self.analyze_domain_selections()
            for uniprot_entry in self.uniprot_entries:
                self.extract_detailed_uniprot_data(uniprot_entry)
            self.commit()

    def setup(self):
        self.uniprot_data_dir = os.path.join(external_data_dirpath, 'UniProt')
        if not os.path.exists(self.uniprot_data_dir):
            os.mkdir(self.uniprot_data_dir)
        self.uniprot_xml_out_filepath = os.path.join(self.uniprot_data_dir, 'uniprot-search.xml')
        self.domain_names_filename = 'selected_domain_names.txt'
        self.now = datetime.datetime.utcnow()
        # get current crawl number
        crawldata_row = models.CrawlData.query.first()
        self.current_crawl_number = crawldata_row.current_crawl_number
        logger.info('Current crawl number: {0}'.format(self.current_crawl_number))

    # TODO ?put this in __init__
    def setup_manual_overrides(self):
        self.manual_overrides = read_manual_overrides()
        self.skip_uniprot_entries = self.manual_overrides.get('skip_uniprot_entries')
        self.pseudodomain_manual_annotations = self.manual_overrides.get('pseudodomains')
        self.ncbi_gene_id_manual_annotations = self.manual_overrides.get('ncbi_gene_ids')
        self.skip_ncbi_gene_entries = self.manual_overrides.get('skip_ncbi_gene_entries')
        self.skip_pdbs = self.manual_overrides.get('skip_pdbs')

    def get_uniprot_data(self):
        if os.path.exists(self.uniprot_xml_out_filepath) and self.use_existing_data:
            logger.info('UniProt XML document found at: {0}'.format(self.uniprot_xml_out_filepath))
        else:
            logger.info('Retrieving new XML document from UniProt website.')
            xml_text = retrieve_uniprot(self.uniprot_query)
            if len(xml_text) == 0:
                raise Exception('UniProt search returned no entries.')
            logger.info('Saving new XML document as: {0}'.format(self.uniprot_xml_out_filepath))
            with open(self.uniprot_xml_out_filepath, 'w') as uniprot_xml_file:
                uniprot_xml_file.write(xml_text + '\n')

        logger.info('Reading UniProt XML document: {0}'.format(self.uniprot_xml_out_filepath))
        self.uniprot_xml = etree.parse(self.uniprot_xml_out_filepath, xml_parser).getroot()

    def extract_uniprot_entries_and_domains(self):
        self.uniprot_entries = self.uniprot_xml.findall('entry')
        self.nuniprot_entries = len(self.uniprot_entries)

        # TODO remove
        if self.uniprot_domain_regex is not None:
            self.selected_domains = self.uniprot_xml.xpath(
                'entry/feature[@type="domain"][match_regex(@description, "{0}")]'.format(
                    self.uniprot_domain_regex
                ),
                extensions={
                    (None, 'match_regex'): xpath_match_regex_case_sensitive
                }
            )
        else:
            self.selected_domains = []

    def analyze_domain_selections(self):
        """
        Prints useful info on the domains selected by uniprot_domain_regex
        """
        selected_domain_names = list(set(
            [d.get('description') for d in self.selected_domains]
        ))

        selected_domain_name_counts = [
            len(self.uniprot_xml.findall('entry/feature[@type="domain"][@description="%s"]' % name))
            for name in selected_domain_names
        ]

        domain_names_str = 'Regex: %s\n' % self.uniprot_domain_regex
        domain_names_str += 'Number of domains matching regex: %d\n\n' % len(self.selected_domains)
        domain_names_str += '= Unique domain names which match regex =\n'

        for i in range(len(selected_domain_names)):
            domain_names_str += '{:^{name_width}s} : {:>{pop_width}d}\n'.format(selected_domain_names[i],
                selected_domain_name_counts[i],
                name_width=max([len(n)+4 for n in selected_domain_names]),
                pop_width=max([len(str(p))+1 for p in selected_domain_name_counts])
            )
        domain_names_str += '\n'
        logger.info(domain_names_str)

        logger.info(
            '(Unique domain names which do not match regex will be output to {0})'.format(
                self.domain_names_filename
            )
        )

        all_domains = self.uniprot_xml.findall('./entry/feature[@type="domain"]')
        domain_names_str += '= Unique domain names which do not match regex =\n'
        nonselected_domain_names = list(set([ d.get('description') for d in all_domains if d.get('description') not in selected_domain_names ]))

        if self.count_nonselected_domain_names:
            nonselected_domain_name_counts = [
                int(
                    self.uniprot_xml.xpath(
                        'count(entry/feature[@type="domain"][@description="{0}"])'.format(name)
                   )
                ) for name in nonselected_domain_names
            ]
            for i in range(len(nonselected_domain_names)):
                domain_names_str += '{:^{name_width}s} : {:>{pop_width}d}\n'.format(
                    nonselected_domain_names[i],
                    nonselected_domain_name_counts[i],
                    name_width=max([len(n)+4 for n in nonselected_domain_names]),
                    pop_width=max([len(str(p))+1 for p in nonselected_domain_name_counts]),
                )
        else:
            for i in range(len(nonselected_domain_names)):
                domain_names_str += '{:^{name_width}s}\n'.format(
                    nonselected_domain_names[i],
                    name_width=max([len(n)+4 for n in nonselected_domain_names]),
                )
        domain_names_str += '\n'

        with open(self.domain_names_filename, 'w') as domain_names_file:
            domain_names_file.write(domain_names_str)

    def extract_detailed_uniprot_data(self, uniprot_entry_node):
        # = IDs and names =
        ac = uniprot_entry_node.findtext('./accession')
        entry_name = uniprot_entry_node.findtext('./name')
        if self.skip_uniprot_entries and entry_name in self.skip_uniprot_entries:
            skip_message = self.skip_uniprot_entries[entry_name]
            logger.info(
                'OVERRIDE: Skipping UniProt entry {0} - reason: {1}'.format(
                    entry_name, skip_message
                )
            )
            return
        recommended_name = uniprot_entry_node.findtext('./protein/recommendedName/fullName')
        gene_name_nodes = uniprot_entry_node.findall('./gene/name')
        gene_name_data = []
        for gene_name_node in gene_name_nodes:
            gene_name = gene_name_node.text
            gene_name_type = gene_name_node.get('type')
            gene_name_obj = models.UniProtGeneName(
                crawl_number=self.current_crawl_number,
                gene_name=gene_name,
                gene_name_type=gene_name_type
            )
            gene_name_data.append(gene_name_obj)

        # = Date entry was last modified in UniProt =
        last_uniprot_update = uniprot_entry_node.get('modified')

        # = Taxonomy =
        uniprot_organism_node = uniprot_entry_node.find('organism')
        ncbi_taxon_id = uniprot_organism_node.find('dbReference[@type="NCBI Taxonomy"]').get('id')
        taxon_name_scientific = uniprot_organism_node.findtext('name[@type="scientific"]')
        taxon_name_common = uniprot_organism_node.findtext('name[@type="common"]')
        lineage = uniprot_organism_node.find('lineage')
        lineage_csv = ','.join([taxon.text for taxon in lineage.getchildren()])

        # = Functions, disease associations, subcellular locations =
        functions = []
        disease_associations = []
        subcellular_locations = []
        for domain in uniprot_entry_node.findall('./comment[@type="function"]'):
            functions.append(
                models.UniProtFunction(
                    crawl_number=self.current_crawl_number,
                    function=domain.findtext('./text')
                )
            )
        for domain in uniprot_entry_node.findall('./comment[@type="disease"]'):
            disease_associations.append(
                models.UniProtDiseaseAssociation(
                    crawl_number=self.current_crawl_number,
                    disease_association=domain.findtext('./text')
                )
            )
        for domain in uniprot_entry_node.findall('./comment[@type="subcellular location"]'):
            subcellular_locations.append(
                models.UniProtSubcellularLocation(
                    crawl_number=self.current_crawl_number,
                    subcellular_location=domain.findtext('./subcellularLocation/location')
                )
            )

        # = Canonical isoform =

        isoforms = []

        # Returned UniProt XML contains sequence data only for the canonical isoform
        uniprot_canonical_sequence_node = uniprot_entry_node.find(
            './sequence[@length][@mass]'
        )
        canonical_sequence = ''.join(uniprot_canonical_sequence_node.text.split())
        canseq_length = uniprot_canonical_sequence_node.get('length')
        canseq_mass = uniprot_canonical_sequence_node.get('mass')
        canseq_date_modified = uniprot_canonical_sequence_node.get('modified')
        canseq_version = uniprot_canonical_sequence_node.get('version')
        uniprot_isoform = models.UniProtIsoform(
            crawl_number=self.current_crawl_number,
            ac=ac+'-1',
            is_canonical=True,
            length=canseq_length,
            mass=canseq_mass,
            date_modified=canseq_date_modified,
            version=canseq_version,
            sequence=canonical_sequence
        )
        # empty list for notes (which do not exist for the canonical sequence)
        isoforms.append((uniprot_isoform, []))

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

        for uniprot_isoform_node in uniprot_entry_node.findall('comment/isoform'):
            isoform_ac = uniprot_isoform_node.findtext('id')
            seq_node = uniprot_isoform_node.find('sequence')
            notes = [
                models.UniProtIsoformNote(
                    crawl_number=self.current_crawl_number, note=node.text
                ) for node in uniprot_isoform_node.findall('note')
                ]
            if seq_node.get('type') != 'displayed':
                uniprot_isoform = models.UniProtIsoform(
                    crawl_number=self.current_crawl_number,
                    ac=isoform_ac,
                    is_canonical=False
                )

            isoforms.append((uniprot_isoform, notes))

        # = UniProt "Protein kinase" domain annotations =
        # XXX TODO Generalize

        # if self.uniprot_domain_regex != None:
        #     selected_domains = uniprot_entry_node.xpath(
        #         'feature[@type="domain"][match_regex(@description, "{0}")]'.format(
        #             self.uniprot_domain_regex
        #         ),
        #         extensions={(None, 'match_regex'): xpath_match_regex_case_sensitive}
        #     )
        # else:
        domains = uniprot_entry_node.findall('feature[@type="domain"]')

        # Skip if no matching domains found
        if len(domains) < 1:
            return

        # Finally, add the domains to the new database
        domain_objs = []
        target_iter = 0
        for domain_id, domain in enumerate(domains):
            # First calculate the PK domain length and sequence
            domain_description = domain.get('description')
            if self.uniprot_domain_regex and re.match(self.uniprot_domain_regex, domain_description):
                is_target_domain = True
                target_id = entry_name + '_D' + str(target_iter)
                target_iter += 1
            else:
                is_target_domain = False
            begin = int(domain.find('./location/begin').get('position'))
            end = int(domain.find('./location/end').get('position'))
            length = end - begin + 1
            domain_seq = canonical_sequence[begin-1:end]

            if (self.pseudodomain_manual_annotations
                and entry_name in self.pseudodomain_manual_annotations
                and domain_description == self.pseudodomain_manual_annotations[entry_name].get('description')
                ):
                pseudodomain_notes = self.pseudodomain_manual_annotations[entry_name].get('message')
                logger.info(
                    'OVERRIDE: Labeling domain "{0}" as a pseudodomain - reason: {1}'.format(
                        target_id,
                        pseudodomain_notes
                    )
                )
                is_pseudodomain = True
            else:
                is_pseudodomain = False

            domain_obj = models.UniProtDomain(
                crawl_number=self.current_crawl_number,
                domain_id=domain_id,
                target_id=target_id if is_target_domain else None,
                is_target_domain=is_target_domain,
                description=domain_description,
                is_pseudodomain=is_pseudodomain,
                pseudodomain_notes=pseudodomain_notes if is_pseudodomain else None,
                begin=begin,
                end=end,
                length=length,
                sequence=domain_seq
            )
            domain_objs.append(domain_obj)

        # = References to other DBs =
        # NCBI Gene
        ncbi_gene_entries = []
        gene_ids = [
            int(domain.get('id')) for domain in uniprot_entry_node.findall('./dbReference[@type="GeneID"]')
        ]

        # manual annotations
        if self.ncbi_gene_id_manual_annotations and entry_name in self.ncbi_gene_id_manual_annotations:
            gene_ids = self.ncbi_gene_id_manual_annotations[entry_name].get('gene_ids')
            gene_ids_message = self.ncbi_gene_id_manual_annotations[entry_name].get('message')
            logger.info(
                'OVERRIDE: Manually annotating Gene IDs for entry {0} - reason: {1}'.format(
                    entry_name, gene_ids_message
                )
            )

        for gene_id in gene_ids:
            # manual override skips
            if self.skip_ncbi_gene_entries and gene_id in self.skip_ncbi_gene_entries:
                skip_gene_id_message = self.skip_ncbi_gene_entries[gene_id]
                logger.info(
                    'OVERRIDE: Skipping Gene ID {0} for entry {1} - reason: {2}'.format(
                        gene_id,
                        entry_name,
                        skip_gene_id_message
                    )
                )
                continue

            ncbi_gene_entries.append(
                models.NCBIGeneEntry(
                    crawl_number=self.current_crawl_number,
                    gene_id=gene_id
                )
            )

        # Ensembl

        # transcript_data = {
        #     'ENSMUST00000003710':
        #         {
        #             'gene':
        #                 'ENSG000...',
        #             'protein':
        #                 'ENSP000...',
        #         }
        # }

        ensembl_transcript_nodes = uniprot_entry_node.findall(
            './dbReference[@type="Ensembl"]'
        )

        ensembl_data = {}
        ensembl_transcript_matched_to_uniprot_isoform = False
        for transcript_node in ensembl_transcript_nodes:
            ensembl_transcript_id = transcript_node.get('id')

            ensembl_gene_nodes = transcript_node.findall('property[@type="gene ID"]')
            if len(ensembl_gene_nodes) > 1:
                logger.info(
                    'WARNING: Ensembl transcript {0} linked with > 1 gene ID'.format(
                        ensembl_transcript_id
                    )
                )
            ensembl_gene_id = ensembl_gene_nodes[0].get('value')

            ensembl_protein_nodes = transcript_node.findall('property[@type="protein sequence ID"]')
            if len(ensembl_protein_nodes) > 1:
                logger.info(
                    'WARNING: Ensembl transcript {0} linked with > 1 protein ID'.format(
                        ensembl_transcript_id
                    )
                )
            ensembl_protein_id = ensembl_protein_nodes[0].get('value')

            uniprot_isoform_molecule_node = transcript_node.find('molecule')
            if uniprot_isoform_molecule_node is not None:
                uniprot_isoform_ac = uniprot_isoform_molecule_node.get('id')
                if uniprot_isoform_ac == isoforms[0][0].ac:
                    ensembl_transcript_matched_to_uniprot_isoform = True
            elif uniprot_isoform_molecule_node is None and ensembl_transcript_matched_to_uniprot_isoform is False:
                uniprot_isoform_ac = isoforms[0][0].ac
                ensembl_transcript_matched_to_uniprot_isoform = True
            else:
                uniprot_isoform_ac = None

            ensembl_data[ensembl_transcript_id] = {
                'gene': ensembl_gene_id,
                'protein': ensembl_protein_id,
                'uniprot_isoform_ac': uniprot_isoform_ac
            }


        # HGNC
        hgnc_entries = []
        hgnc_dbrefs = uniprot_entry_node.findall('./dbReference[@type="HGNC"]')
        for hgnc_dbref in hgnc_dbrefs:
            hgnc_gene_id = hgnc_dbref.get('id')
            approved_symbol = hgnc_dbref.find('property[@type="gene designation"]').get('value')
            hgnc_entries.append(
                models.HGNCEntry(
                    crawl_number=self.current_crawl_number,
                    gene_id=hgnc_gene_id,
                    approved_symbol=approved_symbol
                )
            )

        # = Family information =
        similarity_comments = uniprot_entry_node.xpath('./comment[@type="similarity"]')
        family = False
        for s in similarity_comments:
            for f in kinase_family_uniprot_similarity_text.keys():
                if f in s.findtext('text'):
                    family = kinase_family_uniprot_similarity_text[f]

        # = PDB entries (from UniProt XML) =
        # keep X-ray and NMR structures (not "Model")
        pdbs = uniprot_entry_node.xpath(
            './dbReference[@type="PDB"]/property[@type="method"][@value="X-ray" or @value="NMR"]/..'
        )
        pdb_data = []
        for p in pdbs:
            pdb_id = p.get('id')
            if self.skip_pdbs and pdb_id in self.skip_pdbs:
                skip_pdb_message = self.skip_pdbs[pdb_id]
                logger.info(
                    'OVERRIDE: Skipping PDB {0} for entry {1} - reason: {2}'.format(
                        pdb_id, entry_name, skip_pdb_message
                    )
                )
                continue

            pdb_method = p.find('property[@type="method"]').get('value')
            resolution_node = p.find('property[@type="resolution"]')
            resolution = resolution_node.get('value') if resolution_node != None else None
            chains_span_str = p.find('property[@type="chains"]').get('value')
            chains_span = parse_uniprot_pdbref_chains(chains_span_str)
            chain_data_dicts = []
            for c in chains_span.keys():
                chain_id = c
                pdb_begin = chains_span[c][0]
                pdb_end = chains_span[c][1]
                # Use the begin and end info to decide if this pdb chain includes the pk_domain. But we will get other sequence info from sifts XML files, using gather-pdb.py
                # Have to check against each PK domain
                for domain in domain_objs:
                    pk_begin = domain.begin
                    pk_end = domain.end
                    if (pdb_begin < pk_begin+30) & (pdb_end > pk_end-30):
                        chain_data_dict = models.PDBChain(
                            crawl_number=self.current_crawl_number,
                            chain_id=chain_id,
                            begin=pdb_begin,
                            end=pdb_end
                        )
                        chain_data_dicts.append({
                            'chain_obj': chain_data_dict,
                            'domain_obj': domain
                        })
                    else:
                        continue

            if len(chain_data_dicts) > 0:
                pdb_obj = models.PDBEntry(
                    crawl_number=self.current_crawl_number,
                    pdb_id=pdb_id,
                    method=pdb_method,
                    resolution=resolution
                )
                pdb_data.append({'pdb_obj': pdb_obj, 'chain_data_dicts': chain_data_dicts})

        # ========
        # Construct data objects and add to db
        # ========

        db_entry = models.DBEntry(
            crawl_number=self.current_crawl_number,
            npdbs=len(pdb_data),
            ndomains=len(domain_objs),
            nisoforms=len(isoforms),
            nfunctions=len(functions),
            ndisease_associations=len(disease_associations),
        )
        db.session.add(db_entry)
        uniprot_entry = models.UniProtEntry(
            crawl_number=self.current_crawl_number,
            ac=ac,
            entry_name=entry_name,
            last_uniprot_update=last_uniprot_update,
            ncbi_taxon_id=ncbi_taxon_id,
            db_entry=db_entry,
            recommended_name=recommended_name,
            taxon_name_scientific=taxon_name_scientific,
            taxon_name_common=taxon_name_common,
            lineage=lineage_csv,
        )
        if family:
            uniprot_entry.family = family
        db.session.add(uniprot_entry)
        for function_obj in functions:
            function_obj.db_entry = db_entry
            function_obj.uniprot_entry = uniprot_entry
            db.session.add(function_obj)
        for disease_association_obj in disease_associations:
            disease_association_obj.db_entry = db_entry
            disease_association_obj.uniprot_entry = uniprot_entry
            db.session.add(disease_association_obj)
        for subcellular_location_obj in subcellular_locations:
            subcellular_location_obj.db_entry = db_entry
            subcellular_location_obj.uniprot_entry = uniprot_entry
            db.session.add(subcellular_location_obj)
        for isoform_data in isoforms:
            isoform_obj = isoform_data[0]
            notes = isoform_data[1]
            isoform_obj.db_entry = db_entry
            isoform_obj.uniprot_entry = uniprot_entry
            db.session.add(isoform_obj)
            for note_obj in notes:
                note_obj.uniprotisoform = isoform_obj
                db.session.add(note_obj)
        for domain_obj in domain_objs:
            domain_obj.db_entry = db_entry
            domain_obj.uniprot_entry = uniprot_entry
            db.session.add(domain_obj)
        for pdb_data_dict in pdb_data:
            pdb_obj = pdb_data_dict['pdb_obj']
            chain_data_dicts = pdb_data_dict['chain_data_dicts']
            pdb_obj.db_entry = db_entry
            db.session.add(pdb_obj)
            for chain_data_dict in chain_data_dicts:
                chain_obj = chain_data_dict['chain_obj']
                domain_obj = chain_data_dict['domain_obj']
                chain_obj.pdb_entry = pdb_obj
                chain_obj.uniprot_domain = domain_obj
                db.session.add(chain_obj)
        for gene_name_obj in gene_name_data:
            gene_name_obj.db_entry = db_entry
            db.session.add(gene_name_obj)
        for NCBIGeneEntry in ncbi_gene_entries:
            NCBIGeneEntry.db_entry = db_entry
            db.session.add(NCBIGeneEntry)
        for HGNCEntry in hgnc_entries:
            HGNCEntry.db_entry = db_entry
            db.session.add(HGNCEntry)
        for ensembl_transcript_id in ensembl_data:
            ensembl_gene_id = ensembl_data[ensembl_transcript_id]['gene']
            ensembl_gene_row = models.EnsemblGene(
                crawl_number=self.current_crawl_number,
                gene_id=ensembl_gene_id,
                db_entry=db_entry,
            )
            db.session.add(ensembl_gene_row)

            ensembl_transcript_row = models.EnsemblTranscript(
                crawl_number=self.current_crawl_number,
                transcript_id=ensembl_transcript_id,
                ensembl_gene=ensembl_gene_row,
            )
            ensembl_transcript_uniprot_isoform_ac = ensembl_data[ensembl_transcript_id]['uniprot_isoform_ac']
            if ensembl_transcript_uniprot_isoform_ac is not None:
                matching_uniprot_isoform_obj = [
                    isoform[0] for isoform in isoforms
                    if isoform[0].ac == ensembl_transcript_uniprot_isoform_ac
                ]
                if len(matching_uniprot_isoform_obj) != 0:
                    ensembl_transcript_row.uniprot_isoform = matching_uniprot_isoform_obj[0]
            db.session.add(ensembl_transcript_row)

            ensembl_protein_id = ensembl_data[ensembl_transcript_id]['protein']
            ensembl_protein_row = models.EnsemblProtein(
                crawl_number=self.current_crawl_number,
                protein_id=ensembl_protein_id,
                ensembl_gene=ensembl_gene_row,
                ensembl_transcript=ensembl_transcript_row,
            )
            db.session.add(ensembl_protein_row)

    def commit(self):
        # update db UniProt datestamp
        current_crawl_datestamp_row = models.DateStamps.query.filter_by(
            crawl_number=self.current_crawl_number
        ).first()
        current_crawl_datestamp_row.uniprot_datestamp = self.now
        if self.commit_to_db:
            db.session.commit()


def parse_uniprot_pdbref_chains(chains_span_str):
    """
    Examples of pdbref chains entries to be parsed:
    A=65-119             => {'A':[65,119]}
    A/C/E/G=64-121       => {'A/C/E/G':[64,121]}
    A=458-778, B=764-778 => {'A':[458,778],'B':[764,778]}
    """
    comma_sep = chains_span_str.split(',')
    chains_span = dict()
    for s in comma_sep:
        span = s.split('=')[1]
        begin = int(span.split('-')[0])
        end = int(span.split('-')[1])
        chainids = s.split('=')[0].strip().split('/')
        for c in chainids:
            chains_span[c] = [begin, end]
    return chains_span


def retrieve_uniprot(search_string_query, maxreadlength=100000000):
    """
    Searches the UniProt database given a search string, and retrieves an XML
    file, which is returned as a string.
    maxreadlength is the maximum size in bytes which will be read from the website
    (default 100MB)

    The function also removes the xmlns attribute from <uniprot> tag, as this
    makes xpath searching annoying
    """

    base_url = 'http://www.uniprot.org/uniprot/?query='
    url_request_string = base_url + urllib.quote_plus(search_string_query) + '&force=yes&format=xml'
    response = urllib2.urlopen(url_request_string)
    page = response.read(maxreadlength)
    page = page.replace('xmlns="http://uniprot.org/uniprot" ', '', 1)

    return page


def get_uniprot_mapping(query_data_type, retrieve_data_type, query_data):
    """
NOTE: one gene_id may return multiple uniprotACs.
Use a list of data to query uniprot ID mapping service
and retrieve other types of ID.
query_data can be a list of strings
e.g. get_uniprot_mapping('gene_id','uniprotAC',gene_ids)

Mapping codes can be found here:
http://www.uniprot.org/faq/28#id_mapping_examples
e.g. ACC+ID (from)
     AC (to)
     ID (to)
     P_ENTREZGENEID (both)
    """

    url = 'http://www.uniprot.org/mapping/'

    query_params = {
    'from':query_data_type,
    'to':retrieve_data_type,
    'format':'tab',
    'query':query_data,
    'reviewed':'yes'
    }

    url_query = urllib.urlencode(query_params)
    request = urllib2.Request(url, url_query)
    request.add_header('User-Agent', 'Python contact')
    response = urllib2.urlopen(request)
    page = response.read(200000)

    # page is in two-column format, with one header line. Extract the info we need into a list
    retrieved_data = page.split('\n')[1:-1]
    for l in range(len(retrieved_data)):
        retrieved_data[l] = retrieved_data[l].split('\t')[1]

    return retrieved_data


def retrieve_uniprotACs(query_data):
    """
    Gathers UniProt primary accession IDs by querying UniProt with a given set of GeneIDs.
    Pass it a list of geneID strings
    Returns a list of uniprotACs
    """
    
    url = 'http://www.uniprot.org/uniprot/'
    retrieved_data = []
    
    for q in query_data:
        print q
        query_params = {
        'query':'geneid:%s AND reviewed:yes' % q,
        'format':'tab',
        'columns':'id'
        }
        
        url_query = urllib.urlencode(query_params)
        request = urllib2.Request(url, url_query)
        request.add_header('User-Agent', 'Python contact')
        response = urllib2.urlopen(request)
        page = response.read(200000)
        
        # each page has one header line. Strip this and put the remaining info into a list
        page_stripped = page.split('\n')[1:-1]
        if len(page_stripped) > 1:
            raise Exception , 'got more than one uniprotAC back'
        # if all is well, then the first value in page_stripped was the only uniprotAC returned
        retrieved_data.append(page_stripped[0])
        
    return retrieved_data


def query_uniprot_multiple(query_params):
    """
    Queries UniProt with multiple sets of url params.
    Pass: a list of dicts containing params for urllib.urlencode
    Returns: a list of raw response strings.
    """

    url = 'http://www.uniprot.org/uniprot/'
    retrieved_data = []
    
    for q in query_params:
        
        url_query = urllib.urlencode(q)
        request = urllib2.Request(url, url_query)
        request.add_header('User-Agent', 'Python contact')
        response = urllib2.urlopen(request)
        page = response.read(200000)

        # each page has one header line. Strip this and put the remaining info into a list
        #page_stripped = page.split('\n')[1:-1]
        #retrieved_data.append(page_stripped)
        retrieved_data.append(page)
        
    return retrieved_data


def query_uniprot(query_params):
    """
    Queries UniProt with a set of url params.
    Pass: a dict of params for urllib.urlencode
    Returns: the raw response string
    """

    url = 'http://www.uniprot.org/uniprot/'

    url_query = urllib.urlencode(query_params)
    request = urllib2.Request(url, url_query)
    request.add_header('User-Agent', 'Python contact')
    response = urllib2.urlopen(request)
    page = response.read(200000)

    return page


def retrieve_uniprot_xml(uniprotAC):
    """
    Retrieves a UniProt entry in .xml format
    Pass it a uniprotAC string
    Returns the .xml file as a string
    """
    
    url = 'http://www.uniprot.org/uniprot/'
    
    response = urllib2.urlopen(url+uniprotAC+'.xml')
    page = response.read(200000)
    
    xml_file = page
        
    return xml_file
