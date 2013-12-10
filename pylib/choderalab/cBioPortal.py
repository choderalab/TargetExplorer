# Functions for retrieval of cBioPortal data via the web API
#
# Daniel L. Parton <partond@mskcc.org> - 9 Oct 2013

import urllib2, re, sys, os
from lxml import etree

def get_cancer_studies():
    '''
    Get list of all cancer studies available via cBioPortal
    '''
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
    '''
    For a given cancer study, return a list of available genetic_profile_ids (e.g. 'laml_tcga_pub_mutations')
    '''
    genetic_profile_ids = []

    genetic_profile_url = 'http://www.cbioportal.org/public-portal/webservice.do?cmd=getGeneticProfiles&cancer_study_id=%(cancer_study)s' % vars()
    response = urllib2.urlopen(genetic_profile_url)
    page = response.read(100000000000)
    lines = page.splitlines()
    genetic_profile_ids = []
    for line in lines[1:]:
        words = line.split('\t')
        genetic_profile_ids.append(words[0])

    return genetic_profile_ids


def get_profile_data(case_set_id, genetic_profile_id, entrez_gene_ids):
    '''
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
    '''
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


def retrieve_Mutation_datatxt(case_set_id, genetic_profile_id, gene_ids):
    '''
    Queries cBioPortal for "Mutation" format data, given a list of cBioPortal cancer studies and a list of HGNC Approved gene Symbols.
    Returns the data file as a list of text lines.
    '''
    gene_ids_string = '+'.join(gene_ids)
    mutation_url = 'http://www.cbioportal.org/public-portal/webservice.do?cmd=getProfileData&case_set_id=%(case_set_id)s&genetic_profile_id=%(genetic_profile_id)s&gene_list=%(gene_ids_string)s' % vars()
    response = urllib2.urlopen(mutation_url)
    page = response.read(1000000000)
    lines = page.splitlines()
    return lines


def retrieve_ExtendedMutation_datatxt(case_set_id, genetic_profile_id, gene_ids):
    '''
    Queries cBioPortal for "ExtendedMutation" format data, given a list of cBioPortal cancer studies and a list of HGNC Approved gene Symbols.
    Returns the data file as a list of text lines.
    '''
    gene_ids_string = '+'.join(gene_ids)
    mutation_url = 'http://www.cbioportal.org/public-portal/webservice.do?cmd=getMutationData&case_set_id=%(case_set_id)s&genetic_profile_id=%(genetic_profile_id)s&gene_list=%(gene_ids_string)s' % vars()
    response = urllib2.urlopen(mutation_url)
    page = response.read(1000000000)
    lines = page.splitlines()
    return lines


def retrieve_mutants_xml(output_xml_filepath, cancer_studies, gene_ids, debug=False):
    '''
    Given a list of cBioPortal cancer studies (typically all those available) and a list of HGNC Approved gene Symbols, downloads all mutation data and writes as an XML file.
    IMPORTANT: residue positions are not yet mapped to the canonical UniProt sequence.

    Schema for returned XML:

    <CBPmuts>
      <gene gene_symbol= >
        <mutant source= study= case_id= >
          <mutation mutation_type= aa_change= ... >

    '''

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
        lines = retrieve_ExtendedMutation_datatxt(case_set_id, genetic_profile_id, gene_ids)
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
                mutant_node = gene_nodes_dict[returned_gene_id].find('mutant[@study="%s"][@case_id="%s"]' % (cancer_study, case_id))
                # If not found, create it.
                if mutant_node == None:
                    mutant_node = etree.SubElement(gene_nodes_dict[returned_gene_id], 'mutant')
                    mutant_node.set('source', 'cBioPortal')
                    mutant_node.set('study', cancer_study)
                    mutant_node.set('case_id', case_id)
                    mutant_nodes_by_case_id[case_id] = mutant_node

                # Each mutation gets a mutation_mode
                mutation_node = etree.SubElement(mutant_node, 'mutation')
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
        lines = retrieve_Mutation_datatxt(case_set_id, genetic_profile_id, gene_ids)
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

            if debug:
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
        output_xml_file.write( etree.tostring(results_node, pretty_print=True) )

    # ===========

def match_dual_consecutive_aa_change(aa_change_string):
    '''
    Matches an aa_change string of the type: '324_325LA>FS'
    Returns True for a match, False for a non-match.
    Checks the aas are consecutive.
    '''
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
    '''
    For an aa_change string of the type: '324_325LA>FS'
    Splits into two aa_changes, e.g. ['L324F', 'A325S']
    '''
    # searches for the 'LA>FS' type pattern - putting the regex in parentheses means that the matched pattern will also be returned in the list
    aa_change_string_split = re.split('([A-Z][A-Z]>[A-Z][A-Z])', aa_change_string) # ['324_325', 'LA>FS', '']
    positions = aa_change_string_split[0].split('_') # ['324', '325']
    resnames = aa_change_string_split[1].split('>') # ['LA', 'FS']
    aa_changes = [resnames[0][0] + positions[0] + resnames[1][0], resnames[0][1] + positions[1] + resnames[1][1]]
    return aa_changes

def percent_cases_with_mutations(gene_node):
    '''
    Calculates the percentage of cases with mutations for a given gene node (kinDB XML style), and across all cancer studies with data contained in the node.
    '''
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

