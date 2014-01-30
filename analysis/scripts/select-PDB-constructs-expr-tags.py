# Same as select-PDB_constructs.py but selects only PDB construct sequences which include a 6His tag. This is aimed at avoiding improperly annotated sequences.
#
# Daniel L. Parton <partond@mskcc.org> - 3 Jan 2014
#

import sys, os, openpyxl, re
from lxml import etree
from lxml.builder import E
import TargetExplorer as clab
import Bio.Seq
import Bio.Alphabet
from openpyxl import Workbook

# ===========
# parameters
# ===========

try:
    override_target = sys.argv[ sys.argv.index('-target') + 1 ]
    targets = [override_target]
    ndesired_targets = 1
except ValueError:
    targets = 'All'
    ndesired_targets = 96

results_dir = os.path.join('analysis', 'PDB_construct_selection')

output_Excel_filename = 'PDB_constructs.xlsx'
output_selections_filename = 'PDB_constructs-data.txt'
manual_exceptions_filename = 'manual_exceptions.yaml'
output_Excel_filepath = os.path.join(results_dir, output_Excel_filename)
output_selections_filepath = os.path.join(results_dir, output_selections_filename)
manual_exceptions_filepath = os.path.join(results_dir, manual_exceptions_filename)

html_alignments_dir = os.path.join(results_dir, 'alignments')

css_filename = 'seqlib.cs'
css_filepath = os.path.join(html_alignments_dir, css_filename)

# ===========
# function definitions
# ===========

def generate_html_from_alignment(title, alignment, alignment_IDs, additional_data_fields=None, aa_css_class_list=None):
    '''
    additional_data_fields structure: [ { [PDB_ID]_[PDB_CHAIN_ID] : data }, { [PDB_ID]_[PDB_CHAIN_ID] : data }, ... ] with a separate dict for each data type, and where data is of len(alignment)
    aa_css_class_list can be used to override the CSS classes assigned to each residue. Should be given as a list of lists, with shape: (len(alignment), len(alignment[0])).
    '''
    html_body = E.body()
    html_body.append( E.h2(title) )
    html_table = E.table()
    html_body.append( html_table )

    for i in range(len(alignment)):
        #table = E.table( E.tr( E.td( E.div(alignment_IDs[i],CLASS='ali')), E.td( E.div(id='blank') ), E.td( E.div(id=alignment_IDs[i],CLASS='ali'),nowrap='' ) ) )
        # row
        row = E.tr()

        # alignment ID div
        row.append( E.td( E.div(alignment_IDs[i],CLASS='ali') ) )

        # add any additional data fields (also in divs)
        for d in range(len(additional_data_fields)):
            data = additional_data_fields[d][alignment_IDs[i]]
            if data != None:
                row.append( E.td( E.div(data, CLASS='ali'), nowrap='') )
            else:
                row.append( E.td( E.div('', CLASS='ali'), nowrap='') )

        # format sequence with css classes. Returned as a list of span objects
        if aa_css_class_list != None:
            if aa_css_class_list[i] != None:
                prettyseq = clab.core.seq2pretty_html(alignment[i], aa_css_class_list=aa_css_class_list[i])
            else:
                prettyseq = clab.core.seq2pretty_html(alignment[i])
        else:
            prettyseq = clab.core.seq2pretty_html(alignment[i])

        # set up sequence div
        #row.append( E.td( E.div(id='sequence',CLASS='ali'), nowrap='') )
        #seq_div = row.find('td/div[@id="sequence"]')
        seq_div = E.div(id='sequence',CLASS='ali')
        seq_div.set('style','background-color:#dddddd;letter-spacing:-5px')

        # add sequence to div
        for span in prettyseq:
            seq_div.append(span)

        row.append( E.td( seq_div, nowrap='' ) )

        # add the row to the table
        html_table.append(row)

    return html_body

def match_regex(context, attrib_values, xpath_argument):
    # If no attrib found
    if len(attrib_values) == 0:
        return False
    # If attrib found, then run match against regex
    else:
        return bool( re.match(xpath_argument, attrib_values[0]) )

def write_output(ofilename, html_tree):
    ofile = open(ofilename, 'w')
    ofile.write( etree.tostring(html_tree, pretty_print=True) )
    ofile.close()

def process_target(t):

    target = targets_data[t].keys()[0]
    target_results = {target : []} # results will be returned in this dict
    nmatching_PDB_structures = targets_data[t][target][0]
    print 'Working on target:', target
    # ===========
    # Generate html header
    # ===========
    output_html_tree = E.html(
        E.head(
            E.link()
        ),
        E.body(
        )
    )
    output_html_body = output_html_tree.find('body')
    css_link = output_html_tree.find('head/link')
    css_link.set('type','text/css')
    #css_path = os.path.join('..', '..', '..', 'TargetExplorer', 'seqlib.cs')
    css_path = os.path.join(css_filename)
    css_link.set('href',css_path)
    css_link.set('rel','stylesheet')

    # ===========
    # Get target info from DB
    # ===========

    # get DB entry
    DB_domain = DB_root.find('entry/UniProt/domains/domain[@target_id="%s"]' % target)
    target_domain_len = int(DB_domain.get('length'))
    if target_domain_len > 350 or target_domain_len < 191:
        print 'Target domain length %d. Skipping...' % target_domain_len
        target_results[target] = [0, None, [None, None, None, 0]]
        return target_results
    DB_entry = DB_domain.getparent().getparent().getparent()

    # get IDs
    target_UniProt_entry_name = DB_entry.find('UniProt').get('entry_name')
    target_NCBI_Gene_node = DB_entry.find('NCBI_Gene/entry[@ID]')
    if target_NCBI_Gene_node == None:
        print 'Gene ID not found for target %s' % target_UniProt_entry_name
        target_results[target] = [0, None, [None, None, None, 0]]
        return target_results
    target_NCBI_GeneID = int(target_NCBI_Gene_node.get('ID'))

    if target_NCBI_GeneID not in plasmid_NCBI_GeneIDs:
        print 'Gene ID %s (%s) not found in Harvard plasmid library.' % (target_NCBI_GeneID, target_UniProt_entry_name)
        target_results[target] = [0, None, [None, None, None, 0]]
        return target_results

    # get UniProt canonical isoform sequence
    UniProt_canonseq = clab.core.sequnwrap( DB_entry.findtext('UniProt/isoforms/canonical_isoform/sequence') )

    # ===========
    # Get PDB info from DB
    # ===========

    # get PDB sequences (only those which have the desired expression_system tag) and store in dict e.g. { '3GKZ_B' : 'MGYL...' }
    PDB_matching_seq_nodes = DB_entry.xpath( 'PDB/structure/expression_data[match_regex(@EXPRESSION_SYSTEM, "%s")]/../chain/experimental_sequence/sequence' % desired_expression_system_regex, extensions = { (None, 'match_regex'): match_regex } )
    if len(PDB_matching_seq_nodes) == 0:
        target_results[target] = [0, None, [None, None, None, 0]] # Skip targets with no PDB sequences with matching expression data
        return target_results
    # PDB_seqs structure: { [PDB_ID]_[PDB_CHAIN_ID] : sequence }
    # remove 'X' residues
    PDB_seqs = { seq_node.getparent().getparent().getparent().get('ID') + '_' + seq_node.getparent().getparent().get('ID') : clab.core.sequnwrap(seq_node.text.replace('X', '')) for seq_node in PDB_matching_seq_nodes }

    # PDB_seqs_aln_vs_UniProt_conflicts has same structure and contains the PDB construct seq "aligned" against the UniProt sequence, with conflicts in lower-case (derived using SIFTS info).
    # replace 'x' residues with '-'
    PDB_seqs_aln_vs_UniProt_conflicts = { seq_node.getparent().getparent().getparent().get('ID') + '_' + seq_node.getparent().getparent().get('ID') : clab.core.sequnwrap(seq_node.getparent().find('sequence_aln_conflicts').text.replace('x','-')) for seq_node in PDB_matching_seq_nodes }

    # html_additional_data structure: [ { [PDB_ID]_[PDB_CHAIN_ID] : data }, { [PDB_ID]_[PDB_CHAIN_ID] : data }, ... ] with a separate dict for each data type, and where data is of len(alignment)
    expression_system_data = { seq_node.getparent().getparent().getparent().get('ID') + '_' + seq_node.getparent().getparent().get('ID') : seq_node.getparent().getparent().getparent().find('expression_data').get('EXPRESSION_SYSTEM') for seq_node in PDB_matching_seq_nodes }
    html_additional_data = [{'UniProt': '', 'plasmid_seq': ''}] # No expression_system data to show for UniProt and plasmid sequences
    for key in expression_system_data.keys():
        html_additional_data[0][key] = expression_system_data[key]

    # ===========
    # run MSA (using ClustalO)
    # ===========
    alignment_IDs = ['UniProt', 'plasmid_seq']
    pre_alignment_seqs = [UniProt_canonseq, plasmid_aa_seqs[target_NCBI_GeneID]]
    for PDB_chain_ID in PDB_seqs.keys():
        alignment_IDs.append(PDB_chain_ID)
        pre_alignment_seqs.append(PDB_seqs[PDB_chain_ID])
    aligned_seqs = clab.align.run_clustalo(alignment_IDs, pre_alignment_seqs)
    alignment = [ [alignment_IDs[i], seq] for i, seq in enumerate(aligned_seqs) ] # convert aligned_seqs to this list of 2-element lists
    UniProt_canon_seq_aligned = alignment[0]

    # ===========
    # compare plasmid and PDB seqs with aligned UniProt canonical seq - put non-matching aas in lower case
    # ===========
    for seq_iter in range(len(alignment))[1:]:
        seq = list(alignment[seq_iter][1]) # sequence needs to be mutable so convert to list
        for aa_iter in range(len(seq)):
            if seq[aa_iter] != UniProt_canon_seq_aligned[1][aa_iter]:
                seq[aa_iter] = seq[aa_iter].lower()
        alignment[seq_iter][1] = ''.join(seq)

    # ===========
    # find the target domain within the aligned target domain seq
    # ===========
    target_domain_seq = clab.core.sequnwrap( DB_domain.findtext('sequence') )
    # to do this, we construct a regex which accounts for the possible presence of '-' chars within the target domain sequence.
    target_domain_seq_regex = ''.join([ aa + '-*' for aa in target_domain_seq ])[:-2] # ignore last '-*'
    aligned_UniProt_seq_target_domain_match = re.search(target_domain_seq_regex, UniProt_canon_seq_aligned[1])
    len_aligned_UniProt_seq_target_domain_match = aligned_UniProt_seq_target_domain_match.end() - aligned_UniProt_seq_target_domain_match.start()

    # ===========
    # set up data structures in preparation for sorting
    # ===========
    sorted_alignment = alignment[0:2] # the UniProt seq and plasmid seq stay at the beginning
    # PDB_construct_seqs_aligned structure: [ ['PDB_chain_ID', 'sequence'], ... ]
    PDB_construct_seqs_aligned = alignment[2:]

    # ===========
    # Detect extraneous expression tags - then calculate an "authenticity score", in an attempt to downrank misannotated construct sequences
    # ===========
    # Highest score for 'gh' within non-matching sequence (Abl1 E coli constructs use this at the N-term - comes after a His tag with TEV cleavage site)
    # Next search for His tag
    # Then look for non-matching sequence at the N- or C-term, of >3 aas.
    # Favor N-terminal tags, since QB3 MacroLab plasmids with His tags and TEV cleavage sites perform best in this configuration

    # TODO ?go through and change re.match to re.search and alter regex strings appropriately

    TEV_cleaved_Nterm_regex = '^g[has]m{0,1}g{0,1}s{0,1}[A-Z]+[A-Z]{30}'
    TEV_uncleaved_Nterm_regex = '.*[eE][nNvV][lL][yY]{0,1}[fF][qQ].*[A-Z]{30}'
    TEV_Cterm_regex = '.*[A-Z]{30}.*[eE][nN][lL][yY][fF][qQ]'
    histag_Nterm_regex = '.*[hH]{6}.*[A-Z]+'
    histag_Cterm_regex = '.*[A-Z]+.*[hH]{6}'
    other_extra_seq_Nterm_regex = '.*[a-z]{3}.*[A-Z]{30}'
    other_extra_seq_Cterm_regex = '.*[A-Z]{30}.*[a-z]{3}'

    authenticity_scores = [0] * len(PDB_construct_seqs_aligned)
    # construct_data structure: {alignment_ID : [expression_system, tag_type, terminus, authenticity_score], ...} where terminus can be 'Nterm' or 'Cterm'
    construct_data = { x[0] : None for x in alignment }
    expr_tag_strings = { x[0] : None for x in alignment }
    for i in range(len(PDB_construct_seqs_aligned)):
        ID = PDB_construct_seqs_aligned[i][0]
        PDB_entry_ID = ID.split('_')[0]
        seq = PDB_construct_seqs_aligned[i][1].replace('-', '') # remove '-' from sequence for regex searches
        construct_data[ID] = [ expression_system_data[ID] ]

        # first check for manual exceptions
        manual_exception_behavior = clab.core.parse_nested_dicts(manual_exceptions, [target, PDB_entry_ID, 'authenticity_score', 'behavior'])
        if manual_exception_behavior == 'deprioritize authenticity_score':
            authenticity_scores[i] = -10
            expr_tag_strings[ID] = 'manually deprioritized'
            construct_data[ID] += [None, None, authenticity_scores[i]]
            continue

        # now use regexes to check for the presence of expression tags, and use this information to set the authenticity_scores
        elif re.match(TEV_cleaved_Nterm_regex, seq):
            authenticity_scores[i] = 10
            expr_tag_strings[ID] = 'TEV_cleaved_Nterm'
            construct_data[ID] += ['TEV_cleaved', 'Nterm', authenticity_scores[i]]
            continue
        elif re.match(TEV_uncleaved_Nterm_regex, seq):
            authenticity_scores[i] = 9
            expr_tag_strings[ID] = 'TEV_uncleaved_Nterm'
            construct_data[ID] += ['TEV_uncleaved', 'Nterm', authenticity_scores[i]]
            continue
        elif re.match(TEV_Cterm_regex, seq):
            authenticity_scores[i] = 8
            expr_tag_strings[ID] = 'TEV_Cterm'
            construct_data[ID] += ['TEV', 'Cterm', authenticity_scores[i]]
            continue
        elif re.match(histag_Nterm_regex, seq):
            authenticity_scores[i] = 7
            expr_tag_strings[ID] = 'Histag_Nterm'
            construct_data[ID] += ['Histag', 'Nterm', authenticity_scores[i]]
            continue
        elif re.match(histag_Cterm_regex, seq):
            authenticity_scores[i] = 6
            expr_tag_strings[ID] = 'Histag_Cterm'
            construct_data[ID] += ['Histag', 'Cterm', authenticity_scores[i]]
            continue
        elif re.match(other_extra_seq_Nterm_regex, seq):
            authenticity_scores[i] = 4
            expr_tag_strings[ID] = 'other_extra_seq_Nterm'
            construct_data[ID] += ['other_extra_seq', 'Nterm', authenticity_scores[i]]
            continue
        elif re.match(other_extra_seq_Cterm_regex, seq):
            authenticity_scores[i] = 3
            expr_tag_strings[ID] = 'other_extra_seq_Cterm'
            construct_data[ID] += ['other_extra_seq', 'Cterm', authenticity_scores[i]]
            continue
        else:
            authenticity_scores[i] = 0
            expr_tag_strings[ID] = None
            construct_data[ID] += [None, None, authenticity_scores[i]]
            continue

    html_additional_data.append( expr_tag_strings )

    # ===========
    # calculate the number of aas outside the target domain sequence
    # ===========

    num_aas_outside_target_domain = [0] * len(PDB_construct_seqs_aligned)
    for i in range(len(PDB_construct_seqs_aligned)):
        PDB_chain_ID = PDB_construct_seqs_aligned[i][0]
        # to get the total number of aas in the PDB construct, we take the mapping of that sequence against the UniProt sequence (derived from SIFTS data, not an alignment)
        # get this from the PDB construct sequence mapped (using SIFTS data - not an alignment) against the UniProt seq coords
        PDB_construct_start_UniProt_coords = re.search('[A-Za-z]', PDB_seqs_aln_vs_UniProt_conflicts[PDB_chain_ID]).start()
        PDB_construct_end_UniProt_coords = len(PDB_seqs_aln_vs_UniProt_conflicts[PDB_chain_ID]) - re.search('[A-Z]', PDB_seqs_aln_vs_UniProt_conflicts[PDB_chain_ID][::-1]).start()
        # this calculation will include '-' and 'x' chars within the span
        PDB_construct_len = PDB_construct_end_UniProt_coords - PDB_construct_start_UniProt_coords + 1
        num_aas_outside_target_domain[i] = PDB_construct_len - target_domain_len

    # ===========
    # score the alignment quality for the target domain sequence
    # ===========
    aln_scores = [0] * len(PDB_construct_seqs_aligned)
    for i in range(len(PDB_construct_seqs_aligned)):
        # extract target domain sequences
        UniProt_canon_seq_target_domain_seq = UniProt_canon_seq_aligned[1][ aligned_UniProt_seq_target_domain_match.start() : aligned_UniProt_seq_target_domain_match.end() ]
        PDB_seq_target_domain_seq = PDB_construct_seqs_aligned[i][1][ aligned_UniProt_seq_target_domain_match.start() : aligned_UniProt_seq_target_domain_match.end() ]
        # compare using PAM matrix
        aln_scores[i] = 0 - clab.align.score_aln(UniProt_canon_seq_target_domain_seq, PDB_seq_target_domain_seq) # subtract from 0 to reverse ordering

    # ===========
    # sort the PDB constructs based firstly on the authenticity_score (construct authenticity likelihood), then the alignment score, and finally on the number of aas outside the target domain sequence
    # ===========
    # dict_for_sorting structure: { sequence : (authenticity_score, alignment_score, num_aas_outside_target_domain), ... }
    # negate values for reverse sorting
    dict_for_sorting = { x[0] : (-authenticity_scores[i], aln_scores[i], num_aas_outside_target_domain[i]) for i, x in enumerate(PDB_construct_seqs_aligned) }
    # PDB_construct_seqs_aligned structure: [ ['PDB_chain_ID', 'sequence'], ... ]
    PDB_construct_seqs_aligned = sorted( PDB_construct_seqs_aligned, key = lambda x: dict_for_sorting[x[0]])
    sorted_alignment += PDB_construct_seqs_aligned
    top_PDB_chain_ID = PDB_construct_seqs_aligned[0][0]

    # ===========
    # generate custom set of css classes which will be used to highlight target domain of UniProt sequence in red ('c4')
    # ===========
    aa_css_class_list_UniProt_seq = [None] * len(alignment)
    aa_css_class_list_UniProt_seq[0] = ['bl'] * len(UniProt_canon_seq_aligned[1])
    aa_css_class_list_UniProt_seq[0][ aligned_UniProt_seq_target_domain_match.start() : aligned_UniProt_seq_target_domain_match.end() ] = ['c4'] * len_aligned_UniProt_seq_target_domain_match

    # ===========
    # generate html for the alignment
    # ===========
    sorted_alignment_seqs = [ x[1] for x in sorted_alignment ]
    sorted_alignment_IDs = [ x[0] for x in sorted_alignment ]
    alignment_html = generate_html_from_alignment(target, sorted_alignment_seqs, sorted_alignment_IDs, additional_data_fields=html_additional_data, aa_css_class_list=aa_css_class_list_UniProt_seq)

    # add the alignment html to the main html tree
    for child in alignment_html.getchildren():
        output_html_body.append(child)

    # ===========
    # write to html file
    # ===========
    ofilename = os.path.join(html_alignments_dir, target + '.html')
    write_output(ofilename, output_html_tree)

    # ===========
    # find the "construct target region" (the central portion of the construct with sequence matching the target protein sequence, i.e. excluding expression tags etc.)
    # ===========

    # get this from the PDB construct sequence mapped (using SIFTS data - not an alignment) against the UniProt seq coords
    # desired span is from the first upper-case character to the last

    top_PDB_seq_aln_vs_UniProt_conflicts = PDB_seqs_aln_vs_UniProt_conflicts[top_PDB_chain_ID]
    regex_match_forward = re.search('[A-Z]', top_PDB_seq_aln_vs_UniProt_conflicts)
    top_PDB_seq_aln_vs_UniProt_conflicts_reverse = top_PDB_seq_aln_vs_UniProt_conflicts[::-1]
    regex_match_backward = re.search('[A-Z]', top_PDB_seq_aln_vs_UniProt_conflicts_reverse)
    construct_target_region_start_UniProt_coords = regex_match_forward.start()
    construct_target_region_end_UniProt_coords = len(top_PDB_seq_aln_vs_UniProt_conflicts) - regex_match_backward.end()
    construct_target_region_seq = top_PDB_seq_aln_vs_UniProt_conflicts[ construct_target_region_start_UniProt_coords : construct_target_region_end_UniProt_coords + 1]

    #print construct_target_region_seq
    #print sorted_alignment_seqs[0]
    #print sorted_alignment_seqs[2]

    # now get the construct target region span in the coordinates of the alignment
    # do this by constructing a regex which accounts for the presence of '-' chars, and searching it against the aligned construct sequence
    # ignore '-' chars existing within construct_target_region_seq (which shouldn't be there according to the PDB standard, but frequently are, as SEQRES records often contain the observed sequence rather than the experimental sequence)

    construct_target_region_regex = ''.join([ aa + '-*' for aa in construct_target_region_seq if aa != '-' ])[:-2] # ignore last '-*'
    regex_match = re.search(construct_target_region_regex, sorted_alignment_seqs[2])
    try:
        construct_target_region_start_aln_coords = regex_match.start()
    except Exception as e:
        print UniProt_canon_seq_target_domain_seq
        print construct_target_region_seq
        print sorted_alignment_seqs[2]
        print sorted_alignment_seqs[2].replace('-', '')
        raise e
    construct_target_region_end_aln_coords = regex_match.end() - 1

    # now get the aligned plasmid sequence for the same span
    construct_target_region_aln_plasmid_seq = sorted_alignment_seqs[1][ construct_target_region_start_aln_coords : construct_target_region_end_aln_coords + 1 ]
    # now get the (non-aligned) plasmid sequence; just need to remove '-' chars and convert conflicts back to upper case
    construct_target_region_plasmid_seq = construct_target_region_aln_plasmid_seq.replace('-', '').upper()
    # now get the start and end for the construct target region in the coords of the original plasmid seq
    regex_match = re.search(construct_target_region_plasmid_seq, str(plasmid_aa_seqs[target_NCBI_GeneID]))
    construct_target_region_start_plasmid_coords = regex_match.start()
    construct_target_region_end_plasmid_coords = regex_match.end() - 1

    # ===========
    # add data to targets_results
    # ===========

    # targets_results structure:[ { targetID : data } , ... ] where data is constructed as [nmatching_PDB_structures, top_PDB_chain_ID, construct_data ], where construct_data is constructed as [expression_system, tag_type, tag_loc, authenticity_score]

    # And the top PDB ID to targets_data
    target_results[target].append(nmatching_PDB_structures)
    target_results[target].append(top_PDB_chain_ID)
    # And the construct_data for the top_PDB_chain_ID
    target_results[target].append(construct_data[top_PDB_chain_ID])
    # And the Gene ID
    target_results[target].append(target_NCBI_GeneID)
    # Add the plasmid sequence corresponding to the construct target region, and the start and end aa coordinates
    target_results[target].append(construct_target_region_start_plasmid_coords)
    target_results[target].append(construct_target_region_end_plasmid_coords)
    target_results[target].append(construct_target_region_plasmid_seq)

    return target_results







# ===========
# Main
# ===========

if __name__ == '__main__':

    # ===========
    # Set up directories
    # ===========

    analysis_dir = 'analysis'
    if not os.path.exists(results_dir):
        os.mkdir(results_dir)
    if not os.path.exists(html_alignments_dir):
        os.mkdir(html_alignments_dir)

    # ===========
    # Parse DB
    # ===========

    DB_path = os.path.join('database', 'database-stage.xml')
    DB_root = etree.parse(DB_path).getroot()
    nentries = len(DB_root)

    # ===========
    # Target selection
    # ===========

    if targets == 'All':
        targets = [ domain_node.get('target_id') for domain_node in DB_root.findall('entry/UniProt/domains/domain[@target_id]') ]

    desired_expression_system_regex = 'E.* COLI'

    # ===========
    # Parse manual exceptions file
    # ===========
    import yaml
    manual_exceptions = yaml.load( open(manual_exceptions_filepath, 'r').read() )

    # ===========
    # Parse manual exceptions file
    # ===========
    clab.core.write_css_stylesheet(css_filepath)

    # ===========
    # Sort targets based on number of PDB structures with matching expression systems
    # ===========

    # targets_data and targets_results structure:[ { targetID : data } , ... ] where data will be constructed as [nmatching_PDB_structures, top_PDB_chain_ID, construct_data ]
    targets_data = []

    for target in targets:
        # get DB entry
        DB_domain = DB_root.find('entry/UniProt/domains/domain[@target_id="%s"]' % target)
        DB_entry = DB_domain.getparent().getparent().getparent()

        # get PDB structures
        matching_PDB_structures = DB_entry.xpath( 'PDB/structure/expression_data[match_regex(@EXPRESSION_SYSTEM, "%s")]' % desired_expression_system_regex, extensions = { (None, 'match_regex'): match_regex } )
        targets_data.append( { target : [ len(matching_PDB_structures) ] } )

    # ===========
    # Get Harvard plasmid library sequences
    # ===========

    plasmid_library_dir = os.path.join('external-data', 'Harvard-DFHCC-HIP_human_kinase_collection-pJP1520')
    Harvard_plasmid_library_filepath = os.path.join(plasmid_library_dir, 'Mehle_Kinase_VS1_pJP1520_new_plates.xlsx')

    wb = openpyxl.load_workbook(Harvard_plasmid_library_filepath)
    sheet_ranges = wb.get_sheet_by_name(name = 'Kinase_VS1_pJP1520_new_plates')
    nrows = sheet_ranges.get_highest_row()
    plasmid_aa_seqs = {}
    plasmid_dna_seqs = {}
    for row in range(2,nrows):
        NCBI_GeneID = sheet_ranges.cell('H%d' % row).value    # type(int)
        Symbol = sheet_ranges.cell('I%d' % row).value
        dna_seq = sheet_ranges.cell('L%d' % row).value 
        if len(dna_seq) % 3 != 0:
            print 'WARNING: length of DNA sequence not divisible by 3, for plasmid with Gene Symbol %s and NCBI Gene ID %s' % (Symbol, NCBI_GeneID)
        # Translate to DNA sequence (don't include stop codon)
        aa_seq = Bio.Seq.Seq(dna_seq, Bio.Alphabet.generic_dna).translate(to_stop=True)
        #if len(dna_seq) % 3 != 0:
        #    print aa_seq
        plasmid_aa_seqs[NCBI_GeneID] = aa_seq
        plasmid_dna_seqs[NCBI_GeneID] = dna_seq

    plasmid_NCBI_GeneIDs = plasmid_aa_seqs.keys()

    print ''

    # ===========
    # iterate through targets
    # ===========

    from multiprocessing import Pool
    pool = Pool()
    targets_results = pool.map(process_target, range(len(targets_data)))

    # ===========
    # sort targets based firstly on the PDB construct authenticity score, and secondly on the number of PDB constructs with the desired expression system
    # ===========
    # targets_results structure:[ { targetID : data } , ... ] where data is constructed as [nmatching_PDB_structures, top_PDB_chain_ID, construct_data ], where construct_data is constructed as [expression_system, tag_type, tag_loc, authenticity_score]
    # negate values for reverse sorting
    dict_for_sorting = { target_dict.keys()[0] : [ -target_dict.values()[0][2][3], -target_dict.values()[0][0] ] for target_dict in targets_results }
    # PDB_construct_seqs_aligned structure: [ ['PDB_chain_ID', 'sequence'], ... ]
    targets_results = sorted( targets_results, key = lambda x: dict_for_sorting[x.keys()[0]] )

    # ===========
    # print and write file containing sorted targets with details on PDB structure selection
    # ===========

    print ''

    PDB_selections_text = ''
    PDB_selections_text += '= Targets sorted by number of PDB structures with matching expression system tag =\n\n'
    PDB_selections_text += '%18s  %24s  %16s  %23s  %14s  %30s\n' % ('target', 'nmatching_PDB_structures', 'top_PDB_chain_ID', 'detected_expression_tag', 'authenticity_score', 'expression_system')
    PDB_selections_text += '%18s  %24s  %16s  %23s  %14s  %30s\n' % ('______', '________________________', '________________', '_______________________', '__________________', '_________________')

    ntargets_zero_or_more_PDB = 0
    ntargets_one_or_more_PDB = 0
    ntargets_zero_or_more_PDB_expr_tag = 0
    ntargets_one_or_more_PDB_expr_tag = 0
    for t, target_dict in enumerate(targets_results):
        # Add a line after the 96th target
        if t == 96:
            PDB_selections_text += ('=' * len(PDB_selections_text.split('\n')[2])) + '\n'

        # ignore targets with no matching PDB structures or plasmid sequence, or authenticity_score of 0
        if target_dict == None:
            continue
        # (Note: ror some reason, target_dict.values()[0][0] (nmatching_PDB_structures) needs to be referenced directly - assigning to a variable beforehand results in a value of 0. No idea why.)
        if target_dict.values()[0][0] == 0 or len(target_dict.values()[0]) < 2:
            continue

        target = target_dict.keys()[0]
        top_PDB_chain_ID = target_dict.values()[0][1]
        construct_data = target_dict.values()[0][2]
        if construct_data[1] == None:
            expr_tag_string = 'None'
        else:
            expr_tag_string = construct_data[1] + '_' + construct_data[2]
        PDB_selections_text += '%18s  %24d  %16s  %23s  %14d  %30s\n' % (target, target_dict.values()[0][0], top_PDB_chain_ID, expr_tag_string, construct_data[3], construct_data[0])
        if target_dict.values()[0][0] > 0:
            ntargets_zero_or_more_PDB += 1
            if construct_data[1] != None:
                ntargets_zero_or_more_PDB_expr_tag += 1
        if target_dict.values()[0][0] > 1:
            ntargets_one_or_more_PDB += 1
            if construct_data[1] != None:
                ntargets_one_or_more_PDB_expr_tag += 1

    PDB_selections_text += '\nTotal targets with plasmid and > 0 matching PDB structures: %d\n' % ntargets_zero_or_more_PDB
    PDB_selections_text += 'Total targets with plasmid and > 1 matching PDB structures: %d\n' % ntargets_one_or_more_PDB
    PDB_selections_text += 'Total targets with plasmid and > 0 matching PDB structures and a detected expr_tag: %d\n' % ntargets_zero_or_more_PDB_expr_tag
    PDB_selections_text += 'Total targets with plasmid and > 1 matching PDB structures and a detected expr_tag: %d\n' % ntargets_one_or_more_PDB_expr_tag

    print PDB_selections_text

    with open(output_selections_filepath, 'w') as output_selections_file:
        output_selections_file.write(PDB_selections_text)

    # ===========
    # write .xlsx spreadsheet containing the top 96 targets and necessary details for expression testing
    # ===========

    wb = Workbook()
    ws = wb.get_active_sheet()
    ws.title = output_Excel_filename[0 : output_Excel_filename.find('.xlsx')]

    # Also print data to a csv file
    #csv_filepath = output_Excel_filepath[0 : output_Excel_filepath.find('.xlsx')] + '.txt'
    #with open(csv_filepath, 'w') as csv_file:
    #csv_file.write('Construct index, construct start residue (1-based aa), construct end residue, construct aa sequence, construct DNA sequence (from Harvard DF/HCC library of pJP1520 kinase plasmids)\n')

    headings = ['targetID', 'GeneID', 'expression tag location', 'aa_start', 'aa_end', 'aa_seq', 'dna_start', 'dna_end', 'dna_seq']
    for h in range(len(headings)):
        heading_cell = ws.cell(row=0, column=h)
        heading_cell.value = headings[h]

    t_iter = 0
    ntargets_selected = 0
    while True:

        targetID = targets_results[t_iter].keys()[0]

        #print t_iter, ndesired_targets
        #print targets_results[t_iter]

        if targets_results[t_iter][targetID][0] == 0:
            t_iter += 1
            if t_iter + 1 > ndesired_targets:
                break
            continue

        top_PDB_chain_ID = targets_results[t_iter][targetID][1]

        top_exp_tag_type = targets_results[t_iter][targetID][2]

        target_NCBI_GeneID = targets_results[t_iter][targetID][3]
        construct_aa_start = targets_results[t_iter][targetID][4]
        construct_aa_end = targets_results[t_iter][targetID][5] # this refers to the last aa in 0-based aa coordinates
        construct_aa_seq = targets_results[t_iter][targetID][6]

        construct_dna_start = (construct_aa_start * 3)
        construct_dna_end = (construct_aa_end * 3) + 2 # this refers to the last nucleotide in 0-based nucleotide coordinates

        # Build spreadsheeet

        ID = targetID + '_' + top_PDB_chain_ID
        ID_cell = ws.cell(row=ntargets_selected+1, column=0)
        ID_cell.value = ID

        GeneID_cell = ws.cell(row=ntargets_selected+1, column=1)
        GeneID_cell.value = target_NCBI_GeneID
        
        exp_tag_loc = targets_results[t_iter][targetID][2][2]
        exp_tag_loc_cell = ws.cell(row=ntargets_selected+1, column=2)
        exp_tag_loc_cell.value = exp_tag_loc

        # residue span: 1-based inclusive aa coordinates
        construct_aa_start_cell = ws.cell(row=ntargets_selected+1, column=3)
        construct_aa_end_cell = ws.cell(row=ntargets_selected+1, column=4)
        construct_aa_start_cell.value = construct_aa_start + 1
        construct_aa_end_cell.value = construct_aa_end + 1
        
        # aa_seq
        construct_aa_seq_cell = ws.cell(row=ntargets_selected+1, column=5)
        construct_aa_seq_cell.value = construct_aa_seq

        # DNA span: 1-based inclusive dna nucleotide coordinates
        construct_dna_start_cell = ws.cell(row=ntargets_selected+1, column=6)
        construct_dna_end_cell = ws.cell(row=ntargets_selected+1, column=7)
        construct_dna_start_cell.value = construct_dna_start + 1
        construct_dna_end_cell.value = construct_dna_end + 1
        
        # dna_seq
        orig_plasmid_dna_seq = plasmid_dna_seqs[target_NCBI_GeneID]
        construct_dna_seq = orig_plasmid_dna_seq[ construct_dna_start : construct_dna_end + 1 ]
        construct_dna_seq_cell = ws.cell(row=ntargets_selected+1, column=8)
        construct_dna_seq_cell.value = construct_dna_seq

        if len(construct_dna_seq) % 3 != 0:
            raise Exception, 'modulo 3 of DNA sequence length should be 0. Instead was %d' % (len(construct_dna_seq) % 3)

        t_iter += 1
        ntargets_selected += 1
        if t_iter + 1 > ndesired_targets:
            break

    # Save spreadsheet
    wb.save(output_Excel_filepath)

