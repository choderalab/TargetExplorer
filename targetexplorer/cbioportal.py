import urllib2
import re
import os
import datetime
import json
import gzip
from lxml import etree
from targetexplorer.core import int_else_none, xml_parser, external_data_dirpath, logger
from targetexplorer.utils import json_dump_pretty, set_loglevel
from targetexplorer.oncotator import retrieve_oncotator_mutation_data_as_json, build_oncotator_search_string
from targetexplorer.flaskapp import models, db


external_data_dir = os.path.join(external_data_dirpath, 'cBioPortal')
external_cbioportal_data_filepath = os.path.join(external_data_dir, 'cbioportal-mutations.xml')
external_oncotator_data_filepath = os.path.join(external_data_dir, 'oncotator-data.json.gz')

ensembl_transcript_id_regex = re.compile('(ENS[A-Z]{0,3}T[0-9]{11})')


class GatherCbioportalData(object):
    def __init__(self,
                 use_existing_cbioportal_data=False,
                 use_existing_oncotator_data=False,
                 write_extended_mutation_txt_files=False,
                 run_main=True,
                 commit_to_db=True
                 ):
        self.commit_to_db = commit_to_db
        self.use_existing_cbioportal_data = use_existing_cbioportal_data
        self.use_existing_oncotator_data = use_existing_oncotator_data
        self.write_extended_mutation_txt_files = write_extended_mutation_txt_files
        if os.path.exists(external_oncotator_data_filepath):
            with gzip.open(external_oncotator_data_filepath) as oncotator_data_file:
                self.existing_oncotator_data = json.load(oncotator_data_file)
        else:
            self.existing_oncotator_data = {}

        if not os.path.exists(external_data_dir):
            os.mkdir(external_data_dir)

        self.now = datetime.datetime.utcnow()

        crawldata_row = models.CrawlData.query.first()
        self.current_crawl_number = crawldata_row.current_crawl_number
        print('Current crawl number: {0}'.format(self.current_crawl_number))

        if run_main:
            self.get_hgnc_gene_symbols_from_db()
            self.get_mutation_data_as_xml()
            self.extract_mutation_data()
            self.write_oncotator_data()
            self.finish()

    def get_hgnc_gene_symbols_from_db(self):
        # HGNC gene symbols are necessary to query cBioPortal
        self.db_uniprot_acs = [
            value_tuple[0] for value_tuple
            in models.UniProtEntry.query.filter_by(
                crawl_number=self.current_crawl_number
            ).values(models.UniProtEntry.ac)
        ]
        self.hgnc_gene_symbols = [
            value_tuple[0] for value_tuple
            in models.HGNCEntry.query.filter_by(
                crawl_number=self.current_crawl_number
            ).values(models.HGNCEntry.approved_symbol)
        ]

    def get_mutation_data_as_xml(self):
        if os.path.exists(external_cbioportal_data_filepath) and self.use_existing_cbioportal_data:
            print 'cBioPortal data file found at:', external_cbioportal_data_filepath
        else:
            print 'Retrieving new cBioPortal data file from server...'
            self.cancer_studies = get_cancer_studies()
            retrieve_mutants_xml(
                external_cbioportal_data_filepath,
                self.cancer_studies,
                self.hgnc_gene_symbols,
                write_extended_mutation_txt_files=self.write_extended_mutation_txt_files,
            )

        # import shutil
        # shutil.copy(external_data_filepath, '/Users/partond/tmp')

        self.xmltree = etree.parse(external_cbioportal_data_filepath, xml_parser).getroot()

    def extract_mutation_data(self):
        # for gene_node in self.xmltree.findall('gene'):
        #     cases = gene_node.findall('mutant')
        case_nodes = self.xmltree.findall('gene/case')
        case_rows_by_case_id = {}
        for case_node in case_nodes:
            case_id = case_node.get('case_id')
            print('Extracting mutation data for case {0}'.format(case_id))
            if case_id not in case_rows_by_case_id:
                case_row = models.CbioportalCase(
                    crawl_number=self.current_crawl_number,
                    study=case_node.get('study'),
                    case_id=case_id,
                )
                db.session.add(case_row)
                case_rows_by_case_id[case_id] = case_row
            else:
                case_row = case_rows_by_case_id[case_id]

            mutation_nodes = case_node.findall('mutation')
            for mutation_node in mutation_nodes:
                mutation_type = mutation_node.get('mutation_type')
                chromosome_index = int_else_none(mutation_node.get('chromosome_index'))
                chromosome_startpos = int_else_none(mutation_node.get('chromosome_startpos'))
                chromosome_endpos = int_else_none(mutation_node.get('chromosome_endpos'))
                reference_allele = mutation_node.get('reference_allele')
                variant_allele = mutation_node.get('variant_allele')

                mutation_row = models.CbioportalMutation(
                    crawl_number=self.current_crawl_number,
                    type=mutation_type,
                    cbioportal_aa_change_string=mutation_node.get('aa_change'),
                    mutation_origin=mutation_node.get('mutation_origin'),
                    validation_status=mutation_node.get('validation_status'),
                    functional_impact_score=mutation_node.get('functional_impact_score'),
                    chromosome_index=chromosome_index,
                    chromosome_startpos=chromosome_startpos,
                    chromosome_endpos=chromosome_endpos,
                    reference_dna_allele=reference_allele,
                    variant_dna_allele=variant_allele,
                    cbioportal_case=case_row,
                    in_uniprot_domain=False,
                )

                if mutation_type == 'Missense_Mutation' and None not in [chromosome_index, chromosome_startpos, chromosome_endpos]:
                    oncotator_data = self.get_oncotator_data(
                        chromosome_index,
                        chromosome_startpos,
                        chromosome_endpos,
                        reference_allele,
                        variant_allele
                    )
                    if oncotator_data is not None:
                        reference_aa = oncotator_data['reference_aa']
                        aa_pos = oncotator_data['aa_pos']
                        variant_aa = oncotator_data['variant_aa']
                        ensembl_transcript_id = oncotator_data['ensembl_transcript_id']
                        mutation_row.oncotator_reference_aa = reference_aa
                        mutation_row.oncotator_aa_pos = aa_pos
                        mutation_row.oncotator_variant_aa = variant_aa
                        mutation_row.oncotator_ensembl_transcript_id = ensembl_transcript_id

                        matching_ensembl_transcript_row = models.EnsemblTranscript.query.filter_by(
                            transcript_id=ensembl_transcript_id
                        ).first()
                        if ((matching_ensembl_transcript_row is not None) and
                                (matching_ensembl_transcript_row.uniprot_isoform is not None) and
                                matching_ensembl_transcript_row.uniprot_isoform.is_canonical):
                            mutation_row.db_entry = matching_ensembl_transcript_row.ensembl_gene.db_entry

                            # is mutation within a uniprot domain?
                            matching_uniprot_domains = matching_ensembl_transcript_row.ensembl_gene.db_entry.uniprot_domains.all()
                            for domain in matching_uniprot_domains:
                                if aa_pos >= domain.begin and aa_pos <= domain.end:
                                    if mutation_row.oncotator_reference_aa != mutation_row.cbioportal_aa_change_string[0]:
                                        continue
                                    mutation_row.in_uniprot_domain = True
                                    mutation_row.uniprot_domain = domain

                db.session.add(mutation_row)

    def get_oncotator_data(
            self,
            chromosome_index,
            chromosome_startpos,
            chromosome_endpos,
            reference_allele,
            variant_allele
            ):
        oncotator_search_string = build_oncotator_search_string(
            chromosome_index,
            chromosome_startpos,
            chromosome_endpos,
            reference_allele,
            variant_allele
        )

        if (
                self.use_existing_oncotator_data and self.existing_oncotator_data
                and oncotator_search_string in self.existing_oncotator_data
                ):
                oncotator_data = self.existing_oncotator_data[oncotator_search_string]
        else:
            oncotator_data = retrieve_oncotator_mutation_data_as_json(
                chromosome_index,
                chromosome_startpos,
                chromosome_endpos,
                reference_allele,
                variant_allele
            )
            self.existing_oncotator_data[oncotator_search_string] = oncotator_data

        oncotator_ensembl_transcript_id = oncotator_data.get('transcript_id')
        if oncotator_ensembl_transcript_id is None:
            return None
        if oncotator_ensembl_transcript_id != oncotator_data.get('annotation_transcript'):
            print(
                'WARNING: oncotator transcript_id {0} does not match annotation_transcript {1}'.format(
                    oncotator_data['transcript_id'],
                    oncotator_data['annotation_transcript'],
                )
            )

        # example: 'ENST00000318560.5'
        ensembl_transcript_regex_match = re.match(
            ensembl_transcript_id_regex, oncotator_ensembl_transcript_id
        )
        if ensembl_transcript_regex_match is None:
            return None
        protein_change = oncotator_data['protein_change']
        protein_change_regex_match = re.match('p.([A-Z])([0-9]+)([A-Z])', protein_change)
        if protein_change_regex_match is None:
            return None
        reference_aa, aa_pos, variant_aa = protein_change_regex_match.groups()
        return {
            'ensembl_transcript_id': ensembl_transcript_regex_match.groups(0)[0],
            'reference_aa': reference_aa,
            'aa_pos': int(aa_pos),
            'variant_aa': variant_aa
        }

    def write_oncotator_data(self):
        with gzip.open(external_oncotator_data_filepath, 'w') as external_oncotator_data_file:
            json_dump_pretty(self.existing_oncotator_data, external_oncotator_data_file)

    def finish(self):
        # update db datestamps
        datestamp_row = models.DateStamps.query.filter_by(crawl_number=self.current_crawl_number).first()
        datestamp_row.cbioportal_datestamp = self.now
        if self.commit_to_db:
            db.session.commit()
        print 'Done.'


def get_cancer_studies():
    """
    Get list of all cancer studies available via cBioPortal
    """
    cancer_studies_url = 'http://www.cbioportal.org/public-portal/webservice.do?cmd=getCancerStudies'

    response = urllib2.urlopen(cancer_studies_url)
    page = response.read(100000000000)
    lines = page.splitlines()

    cancer_studies = []
    for line in lines[1:]:
        words = line.split('\t')
        cancer_studies.append(words[0])

    return cancer_studies


def get_genetic_profile_ids(cancer_study):
    """
    For a given cancer study, return a list of available genetic_profile_ids (e.g. 'laml_tcga_pub_mutations')
    """
    genetic_profile_url = 'http://www.cbioportal.org/public-portal/webservice.do?cmd=getGeneticProfiles&cancer_study_id={0}'.format(cancer_study)
    response = urllib2.urlopen(genetic_profile_url)
    page = response.read(100000000000)
    lines = page.splitlines()
    genetic_profile_ids = []
    for line in lines[1:]:
        words = line.split('\t')
        genetic_profile_ids.append(words[0])

    return genetic_profile_ids


def get_profile_data(case_set_id, genetic_profile_id, entrez_gene_ids):
    """
    For each cancer study and each genetic_profile_id, get mutation data for all genes

    ___Note on case_set_ids___
    These comprise cancer_study + '_' + subset.
    subset may be 'all', 'sequenced' (i.e. only include samples which have undergone sequencing), 'acgh' (i.e. only include samples which have undergone aCGH), 'cnaseq' (i.e. only include samples which have undergone sequencing and aCGH) or various others.
    For more options, see the drop-down menu under "Select Patient/Case Set" here: http://www.cbioportal.org/public-portal/index.do
    The subset type is important to consider when calculating the '% in cohort' statistic.

    ___Note on genetic_profile_ids___
    These comprise cancer_study + '_' + profile_type
    profile_type is e.g. 'mutations', 'gistic', 'methylation'
    Some cancer studies may only include a limited set of profile_types. For CNAs, some carry the 'gistic' type, some the 'cna_rae' type or 'cna_consensus'. TODO look more at these. Other types include 'log2CNA'
    """
    # TODO probably deprecate this once get_mutations and get_CNAs are ready

    gene_list_string = '+'.join( [ str(gene_id) for gene_id in entrez_gene_ids ] )
    mutation_url = 'http://www.cbioportal.org/public-portal/webservice.do?cmd=getProfileData&case_set_id=%(case_set_id)s&genetic_profile_id=%(genetic_profile_id)s&gene_list=%(gene_list_string)s' % vars()
    response = urllib2.urlopen(mutation_url)
    page = response.read(100000000000)
    lines = page.splitlines()
    gene_mutations = {}
    for line in lines[3:]:
        words = line.split('\t')
        gene = words[1]
        mutations = words[2:]
        gene_mutations[gene] = mutations

    return gene_mutations


def retrieve_mutants_xml(output_xml_filepath, cancer_studies, gene_ids,
                         write_extended_mutation_txt_files=False,
                         verbose=False
                         ):
    """
    Given a list of cBioPortal cancer studies (typically all those available) and a list of HGNC
    Approved gene Symbols, downloads all mutation data and writes as an XML file.
    IMPORTANT: residue positions are not yet mapped to the canonical UniProt sequence.

    Schema for returned XML:

    <CBPmuts>
      <gene gene_symbol= >
        <mutant source= study= case_id= >
          <mutation mutation_type= aa_change= ... >

    """

    # XML root node
    results_node = etree.Element('CBPmuts')

    # Use this dict to reference gene nodes by gene_id (may be quicker than using results_node.find())
    gene_nodes_dict = {}
    for gene_id in gene_ids:
        gene_node = etree.SubElement(results_node, 'gene')
        gene_node.set('gene_symbol', gene_id)
        gene_nodes_dict[gene_id] = gene_node

    # For each cancer_study
    for cancer_study in cancer_studies:

        # ==============
        # First get "ExtendedMutation" data
        # --------------
        case_set_id = cancer_study + '_sequenced'
        genetic_profile_id = cancer_study + '_mutations'
        print 'Retrieving ExtendedMutation data from cBioPortal for study %s...' % cancer_study

        if write_extended_mutation_txt_files:
            txt_output_filepath = os.path.join(external_data_dir, cancer_study+'.txt')
        else:
            txt_output_filepath = False

        lines = retrieve_extended_mutation_datatxt(
            case_set_id,
            genetic_profile_id,
            gene_ids,
            write_to_filepath=txt_output_filepath,
        )
        if lines == ['Error: Problem when identifying a cancer study for the request.']:
            print 'WARNING: case_set_id "%s" not available - probably means that sequencing data from the underlying cancer study is not yet available. Skipping this case set.' % case_set_id
            continue
        if lines[0][0:25] == '# Warning:  Unknown gene:':
            print lines[0]
            raise Exception
        print 'Done retrieving ExtendedMutation data from cBioPortal.'

        # These dicts will be used later to assign percent_in_cohort values for a given mutation by matching case_id and aa_change
        mutant_nodes_by_case_id = {}
        mutation_nodes_by_case_id_and_aa_change = {}

        # First line is a header
        # Second line describes the tab-separated data columns:
        # entrez_gene_id  gene_symbol case_id sequencing_center   mutation_status mutation_type   validation_status   amino_acid_change   functional_impact_score xvar_link   xvar_link_pdb   xvar_link_msa   chr start_position  end_position    reference_allele    variant_allele  genetic_profile_id
        # From third line onwards, each line is associated with a single mutant, e.g.:
        # 1956  EGFR    TCGA-16-1048    broad.mit.edu   Somatic Missense_Mutation   NA  F254I   M   [getma.org/...snip]   [getma.org/...snip] [getma.org/...snip] 7   55221716    55221716    T   A   gbm_tcga_mutations
        # EXCEPT for dual mutations of consecutive aas, e.g.:
        # 51231 VRK3    TCGA-AA-A01V    broad.mit.edu   Somatic Missense_Mutation   Unknown 324_325LA>FS    NA  NA  NA  NA  19  50493019    50493020    CC  AA  coadread_tcga_pub_mutations
        for line in lines[2:]:
            words = line.split('\t')
            returned_gene_id = words[1]
            case_id = words[2]
            mutation_status = words[4] # 'Somatic', ?'Germline'
            mutation_type = words[5] # 'Missense_Mutation', ...
            validation_status = words[6] # 'NA', ???
            aa_changes = [ words[7] ] # 'F254I', '324_325LA>FS'
            functional_impact_score = words[8]
            chromosome_index = words[12]
            chromosome_startpos = words[13]
            chromosome_endpos = words[14]
            reference_allele = words[15]
            variant_allele = words[16]
            #returned_cancer_study = words[-1][ 0 : words[-1].rfind('_mutations') ]

            # _______________
            # Some exceptions
            # ---------------

            if aa_changes[0] == 'MUTATED':
                aa_changes[0] = 'Unknown'

            # Dual mutations of consecutive aas
            if mutation_type == 'Missense_Mutation' and match_dual_consecutive_aa_change(aa_changes[0]):
                aa_changes = split_dual_consecutive_aa_change(aa_changes[0])

            if returned_gene_id == 'C9ORF96':
                returned_gene_id = 'C9orf96'

            # _______________

            for aa_change in aa_changes:
                # Look for a mutant node (corresponding to this study and case_id). This will only exist if a mutation has already been added for this mutant.
                case_node = gene_nodes_dict[returned_gene_id].find('case[@study="%s"][@case_id="%s"]' % (cancer_study, case_id))
                # If not found, create it.
                if case_node == None:
                    case_node = etree.SubElement(gene_nodes_dict[returned_gene_id], 'case')
                    case_node.set('source', 'cBioPortal')
                    case_node.set('study', cancer_study)
                    case_node.set('case_id', case_id)
                    mutant_nodes_by_case_id[case_id] = case_node

                # Each mutation gets a mutation_mode
                mutation_node = etree.SubElement(case_node, 'mutation')
                mutation_node.set('mutation_origin', mutation_status)
                mutation_node.set('mutation_type', mutation_type)
                mutation_node.set('validation_status', validation_status)
                mutation_node.set('aa_change', aa_change)
                mutation_node.set('functional_impact_score', functional_impact_score)
                mutation_node.set('chromosome_index', chromosome_index)
                mutation_node.set('chromosome_startpos', chromosome_startpos)
                mutation_node.set('chromosome_endpos', chromosome_endpos)
                mutation_node.set('reference_allele', reference_allele)
                mutation_node.set('variant_allele', variant_allele)

                # This dict is used later to assign percent_in_cohort values to mutations by matching case_id and aa_change
                mutation_nodes_by_case_id_and_aa_change[case_id + aa_change] = mutation_node

        # ============

        # ============
        # Now get the non-extended "Mutation" format data - this includes non-mutated samples, thus allowing calculation of percent_in_cohort values
        # ------------

        print 'Retrieving Mutation data from cBioPortal for study %s...' % cancer_study
        lines = retrieve_mutation_datatxt(case_set_id, genetic_profile_id, gene_ids)
        print 'Done retrieving Mutation data from cBioPortal.'

        # First two lines are header info
        # Third line contains the case_ids, tab-separated
        case_ids_returned = lines[2].split('\t')[2:]

        num_in_cohort = len(case_ids_returned)

        # Fourth line onwards: first column is Entrez Gene ID, second is gene_id, third onwards are the aa_change strings (in order of case_ids from the third line)
        # 1956    EGFR    NaN NaN C620Y ...
        # Multiple mutations for a given case are split by ',', except for dual consecutive mutations which are annotated in the form '324_325LA>FS'
        for line in lines[3:]:
            words = line.split('\t')
            returned_gene_id = words[1]
            aa_change_strings = words[2:]

            if verbose:
                print returned_gene_id

            # Make a flat list of aa changes (with multiple-mutation aa_changes split into separate elements)
            all_aa_changes_this_cohort = []
            # Also require a list of case_ids for each individual aa_change
            case_ids_returned_by_aa_change = []
            for a, aa_change_string in enumerate(aa_change_strings):
                aa_changes_split_comma = aa_change_string.split(',') # May result in ['R23K', '56_57EG>AH']
                aa_changes = []
                for b in range(len(aa_changes_split_comma)):
                    if match_dual_consecutive_aa_change(aa_changes_split_comma[b]):
                        aa_changes += split_dual_consecutive_aa_change(aa_changes_split_comma[b])
                    else:
                        aa_changes.append( aa_changes_split_comma[b] )

                all_aa_changes_this_cohort += aa_changes
                case_ids_returned_by_aa_change += ( [case_ids_returned[a]] * len(aa_changes) )

            # Number of aas with any aa_change for this gene and cancer study
            num_in_cohort_any_aa_change = sum( [ 1 for x in aa_change_strings if x != 'NaN' ] )
            percent_in_cohort_any_aa_change = float(num_in_cohort_any_aa_change) / num_in_cohort * 100.

            for a, aa_change in enumerate(all_aa_changes_this_cohort):
                if aa_change not in ['NaN', 'MUTATED']:
                    num_in_cohort_this_aa_change = sum( [ 1 for x in all_aa_changes_this_cohort if x == aa_change ] )
                    # Calculate percent_in_cohort
                    percent_in_cohort_this_aa_change = float(num_in_cohort_this_aa_change) / num_in_cohort * 100.
                    mutation_nodes_by_case_id_and_aa_change[ case_ids_returned_by_aa_change[a] + aa_change ].set('percent_in_cohort_this_aa_change', '%.3f' % percent_in_cohort_this_aa_change)
                    mutation_nodes_by_case_id_and_aa_change[ case_ids_returned_by_aa_change[a] + aa_change ].getparent().set('percent_in_cohort_any_aa_change', '%.3f' % percent_in_cohort_any_aa_change)
                    mutation_nodes_by_case_id_and_aa_change[ case_ids_returned_by_aa_change[a] + aa_change ].getparent().set('num_in_cohort', str(num_in_cohort))

        # ============

    # ===========
    # Write XML to file
    # -----------

    with open(output_xml_filepath, 'w') as output_xml_file:
        output_xml_file.write(etree.tostring(results_node, pretty_print=True))

    # ===========


def retrieve_mutation_datatxt(case_set_id,
                              genetic_profile_id,
                              gene_ids,
                              portal_version='public-portal',
                              verbose=False,
                              ):
    """
    Queries cBioPortal for "Mutation" format data, given a list of cBioPortal cancer studies and a list of HGNC Approved gene Symbols.
    Returns the data file as a list of text lines.
    """
    gene_ids_string = '+'.join(gene_ids)
    mutation_url = 'http://www.cbioportal.org/{0}/' \
                   'webservice.do' \
                   '?cmd=getProfileData' \
                   '&case_set_id={1}' \
                   '&genetic_profile_id={2}' \
                   '&gene_list={3}'.format(
                       portal_version,
                       case_set_id,
                       genetic_profile_id,
                       gene_ids_string
                   )
    if verbose:
        set_loglevel('debug')
        logger.debug(mutation_url)
    response = urllib2.urlopen(mutation_url)
    page = response.read(1000000000)
    lines = page.splitlines()
    return lines


def retrieve_extended_mutation_datatxt(case_set_id,
                                       genetic_profile_id,
                                       gene_ids,
                                       portal_version='public-portal',
                                       write_to_filepath=False
                                       ):
    """
    Queries cBioPortal for "ExtendedMutation" format data, given a list of cBioPortal cancer studies and a list of HGNC Approved gene Symbols.
    Returns the data file as a list of text lines.

    Parameters
    ----------
    portal_version: str
        'public-portal': use only public cBioPortal data
        'private': use private cBioPortal data
    write_to_filepath: str (or False)
    """
    gene_ids_string = '+'.join(gene_ids)
    mutation_url = 'http://www.cbioportal.org/{0}/' \
                   'webservice.do' \
                   '?cmd=getMutationData' \
                   '&case_set_id={1}' \
                   '&genetic_profile_id={2}' \
                   '&gene_list={3}'.format(
                       portal_version,
                       case_set_id,
                       genetic_profile_id,
                       gene_ids_string
                   )
    response = urllib2.urlopen(mutation_url)
    page = response.read(1000000000)
    if write_to_filepath:
        with open(write_to_filepath, 'w') as ofile:
            ofile.write(page)
    lines = page.splitlines()
    return lines


def match_dual_consecutive_aa_change(aa_change_string):
    """
    Matches an aa_change string of the type: '324_325LA>FS'
    Returns True for a match, False for a non-match.
    Checks the aas are consecutive.
    """
    # searches for the 'LA>FS' type pattern - putting the regex in parentheses means that the matched pattern will also be returned in the list
    match = re.match('^[0-9]+_[0-9]+[A-Z][A-Z]>[A-Z][A-Z]', aa_change_string) # ['324_325', 'LA>FS', '']

    if not match:
        return False

    aa_change_string_split = re.split('([A-Z][A-Z]>[A-Z][A-Z])', aa_change_string) # ['324_325', 'LA>FS', '']
    positions = aa_change_string_split[0].split('_') # ['324', '325']

    # check positions are consecutive
    if int(positions[1]) == (int(positions[0]) + 1):
        return True
    else:
        return False


def split_dual_consecutive_aa_change(aa_change_string):
    """
    For an aa_change string of the type: '324_325LA>FS'
    Splits into two aa_changes, e.g. ['L324F', 'A325S']
    """
    # searches for the 'LA>FS' type pattern - putting the regex in parentheses means that the matched pattern will also be returned in the list
    aa_change_string_split = re.split('([A-Z][A-Z]>[A-Z][A-Z])', aa_change_string) # ['324_325', 'LA>FS', '']
    positions = aa_change_string_split[0].split('_') # ['324', '325']
    resnames = aa_change_string_split[1].split('>') # ['LA', 'FS']
    aa_changes = [resnames[0][0] + positions[0] + resnames[1][0], resnames[0][1] + positions[1] + resnames[1][1]]
    return aa_changes


def percent_cases_with_mutations(gene_node):
    """
    Calculates the percentage of cases with mutations for a given gene node (kinDB XML style), and across all cancer studies with data contained in the node.
    """
    cancer_studies = set([x.get('study') for x in gene_node.findall('mutants/mutant')])
    if len(cancer_studies) == 0:
        return None
    nmutants = 0
    num_in_cohort = 0
    for cancer_study in cancer_studies:
        try:
            nmutants += len( gene_node.findall('mutants/mutant[@study="%(cancer_study)s"]' % vars()) )
            num_in_cohort += int( gene_node.find('mutants/mutant[@study="%(cancer_study)s"]' % vars()).get('num_in_cohort') )
        except TypeError:
            continue

    percent_cases_with_mutations = (float(nmutants) / float(num_in_cohort)) * 100.
    return percent_cases_with_mutations

