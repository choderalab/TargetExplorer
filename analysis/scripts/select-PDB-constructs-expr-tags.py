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

try:
    override_target = sys.argv[ sys.argv.index('-target') + 1 ]
    targets = [override_target]
except ValueError:
    targets = 'All'

output_selections_filepath = os.path.join('analysis', 'PDB_construct_selection', 'PDB_construct_selections.txt')

# ===========
# Function definitions
# ===========

def generate_html_from_alignment(title, alignment, alignment_IDs, additional_data_fields=None, aa_css_class_list=None):
    '''
    additional_data_fields structure: [ { [PDB_ID]_[PDB_CHAIN_ID] : data }, { [PDB_ID]_[PDB_CHAIN_ID] : data }, ... ]
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
# Sort targets based on number of PDB structures with matching expression systems
# ===========

targets_nmatching_PDB_structures = []

for target in targets:
    # get DB entry
    DB_domain = DB_root.find('entry/UniProt/domains/domain[@target_id="%s"]' % target)
    DB_entry = DB_domain.getparent().getparent().getparent()

    # get PDB structures
    matching_PDB_structures = DB_entry.xpath( 'PDB/structure/expression_data[match_regex(@EXPRESSION_SYSTEM, "%s")]' % desired_expression_system_regex, extensions = { (None, 'match_regex'): match_regex } )
    targets_nmatching_PDB_structures.append( { target : [ len(matching_PDB_structures) ] } )

# targets_data structure:[ { targetID : data } , ... ] where data will be constructed as [nmatching_PDB_structures, TODO ]
targets_data = sorted( targets_nmatching_PDB_structures, key = lambda x:x.values()[0][0], reverse=True )

# Print the sorted targets and the number of PDB structures with matching expression system
print ''
print '= Targets sorted by number of PDB structures with matching expression system tag =\n'
print '%18s %24s' % ('target', 'nmatching_PDB_structures')
print '%18s %24s' % ('______', '________________________')
for target_dict in targets_data:
    target = target_dict.keys()[0]
    nmatching_PDB_structures = target_dict.values()[0][0]
    if nmatching_PDB_structures == 0:
        continue
    print '%18s %24d' % (target, nmatching_PDB_structures)

print ''
print 'Total targets with > 0 matching PDB structures:', sum([1 for target_dict in targets_data if target_dict.values()[0][0] > 0])
print 'Total targets with > 1 matching PDB structures', sum([1 for target_dict in targets_data if target_dict.values()[0][0] > 1])

# ===========
# Get Harvard plasmid library sequences
# ===========

plasmid_library_dir = os.path.join('external-data', 'Harvard-DFHCC-HIP_human_kinase_collection-pJP1520')
Harvard_plasmid_library_filepath = os.path.join(plasmid_library_dir, 'Mehle_Kinase_VS1_pJP1520_new_plates.xlsx')

wb = openpyxl.load_workbook(Harvard_plasmid_library_filepath)
sheet_ranges = wb.get_sheet_by_name(name = 'Kinase_VS1_pJP1520_new_plates')
nrows = sheet_ranges.get_highest_row()
plasmid_aa_seqs = {}
for row in range(2,nrows):
    NCBI_GeneID = sheet_ranges.cell('H%d' % row).value    # type(int)
    Symbol = sheet_ranges.cell('I%d' % row).value
    DNA_seq = sheet_ranges.cell('L%d' % row).value 
    # Translate to DNA sequence (don't include stop codon)
    aa_seq = Bio.Seq.Seq(DNA_seq, Bio.Alphabet.generic_dna).translate(to_stop=True)
    plasmid_aa_seqs[NCBI_GeneID] = aa_seq

plasmid_NCBI_GeneIDs = plasmid_aa_seqs.keys()

# ===========
# Iterate through targets
# ===========

print ''
print '= Carrying out sequence alignments ='

analysis_dir = 'analysis'
output_toplevel_dir = os.path.join(analysis_dir, 'PDB_construct_selection')
html_alignments_dir = os.path.join(output_toplevel_dir, 'alignments')
if not os.path.exists(output_toplevel_dir):
    os.mkdir(output_toplevel_dir)
if not os.path.exists(html_alignments_dir):
    os.mkdir(html_alignments_dir)

for t in range(len(targets_data)):
    target = targets_data[t].keys()[0]
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
    css_path = os.path.abspath( os.path.join('TargetExplorer', 'seqlib.cs') )
    css_link.set('href',css_path)
    css_link.set('rel','stylesheet')

    # get DB entry
    DB_domain = DB_root.find('entry/UniProt/domains/domain[@target_id="%s"]' % target)
    DB_entry = DB_domain.getparent().getparent().getparent()

    # get IDs
    target_UniProt_entry_name = DB_entry.find('UniProt').get('entry_name')
    target_NCBI_Gene_node = DB_entry.find('NCBI_Gene/entry[@ID]')
    if target_NCBI_Gene_node == None:
        print 'Gene ID not found for target %s' % target_UniProt_entry_name
        continue
    target_NCBI_GeneID = int(target_NCBI_Gene_node.get('ID'))

    if target_NCBI_GeneID not in plasmid_NCBI_GeneIDs:
        print 'Gene ID %s (%s) not found in Harvard plasmid library.' % (target_NCBI_GeneID, target_UniProt_entry_name)
        continue

    # get UniProt canonical isoform sequence
    UniProt_canonseq = clab.core.sequnwrap( DB_entry.findtext('UniProt/isoforms/canonical_isoform/sequence') )

    # get PDB sequences (only those which have the desired expression_system tag) and store in dict e.g. { '3GKZ_B' : 'MGYL...' }
    PDB_matching_seq_nodes = DB_entry.xpath( 'PDB/structure/expression_data[match_regex(@EXPRESSION_SYSTEM, "%s")]/../chain/experimental_sequence/sequence' % desired_expression_system_regex, extensions = { (None, 'match_regex'): match_regex } )
    if len(PDB_matching_seq_nodes) == 0:
        continue # Skip targets with no PDB sequences with matching expression data
    # PDB_seqs structure: { [PDB_ID]_[PDB_CHAIN_ID] : sequence }
    PDB_seqs = { seq_node.getparent().getparent().getparent().get('ID') + '_' + seq_node.getparent().getparent().get('ID') : clab.core.sequnwrap(seq_node.text) for seq_node in PDB_matching_seq_nodes }

    # additional_data structure: [ { [PDB_ID]_[PDB_CHAIN_ID] : data }, { [PDB_ID]_[PDB_CHAIN_ID] : data }, ... ]
    expression_system_data = { seq_node.getparent().getparent().getparent().get('ID') + '_' + seq_node.getparent().getparent().get('ID') : seq_node.getparent().getparent().getparent().find('expression_data').get('EXPRESSION_SYSTEM') for seq_node in PDB_matching_seq_nodes }
    additional_data = [{'UniProt': '', 'plasmid_seq': ''}] # No expression_system data to show for UniProt and plasmid sequences
    for key in expression_system_data.keys():
        additional_data[0][key] = expression_system_data[key]

    # ===========
    # run MSA (using ClustalO)
    # ===========
    alignment_IDs = ['UniProt', 'plasmid_seq']
    alignment_seqs = [UniProt_canonseq, plasmid_aa_seqs[target_NCBI_GeneID]]
    for PDB_ID in PDB_seqs.keys():
        alignment_IDs.append(PDB_ID)
        alignment_seqs.append(PDB_seqs[PDB_ID])
    aligned_seqs = clab.align.run_clustalo(alignment_IDs, alignment_seqs)
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
    # PDB_construct_seqs_aligned structure: [ ['PDB_ID', 'sequence'], ... ]
    PDB_construct_seqs_aligned = alignment[2:]

    # ===========
    # Calculate construct sequence "authenticity score" - attempt to downrank misannotated sequences
    # ===========
    # Highest score for 'ghm' within non-matching sequence (Abl1 E coli constructs use this at the N-term - comes after a His tag with TEV cleavage site)
    # Next search for His tag
    # Then look for non-matching sequence at the N- or C-term, of >3 aas.
    ghm_regex = '.*ghm[A-Z]+|.*[A-Z]+mhg'
    histag_regex = '.*(h|H){6,999}'
    nonmatching_seq_regex = '.*[a-z]{3,9999}[A-Z]+|.*[A-Z]+[a-z]{3,9999}'

    # more negative scores rank higher
    authenticity_scores = [0] * len(PDB_construct_seqs_aligned)
    for i in range(len(PDB_construct_seqs_aligned)):
        seq = PDB_construct_seqs_aligned[i][1]
        if re.match(ghm_regex, seq):
            authenticity_scores[i] = -10
            continue
        elif re.match(histag_regex, seq):
            authenticity_scores[i] = -8
            continue
        elif re.match(nonmatching_seq_regex, seq):
            authenticity_scores[i] = -4
            continue

    # ===========
    # calculate the number of aas outside the target domain sequence
    # ===========
    num_aas_outside_target_domain = [0] * len(PDB_construct_seqs_aligned)
    for i in range(len(PDB_construct_seqs_aligned)):
        # calculate number of aas present outside the target domain
        for aa_iter in range(len(PDB_construct_seqs_aligned[i][1])):
            if PDB_construct_seqs_aligned[i][1][aa_iter] != '-' and ( aa_iter < aligned_UniProt_seq_target_domain_match.start() or aa_iter > aligned_UniProt_seq_target_domain_match.end() ):
                num_aas_outside_target_domain[i] += 1

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
    # rank the PDB constructs based firstly on the authenticity likelihood, then the alignment score, and finally on the number of aas outside the target domain sequence
    # ===========
    # dict_for_sorting structure: { 'sequence' : (authenticity_score, alignment_score, num_aas_outside_target_domain), ... }
    dict_for_sorting = { x[1] : (authenticity_scores[i], aln_scores[i], num_aas_outside_target_domain[i]) for i, x in enumerate(PDB_construct_seqs_aligned) }
    # PDB_construct_seqs_aligned structure: [ ['PDB_ID', 'sequence'], ... ]
    PDB_construct_seqs_aligned = sorted( PDB_construct_seqs_aligned, key = lambda x: dict_for_sorting[x[1]])
    sorted_alignment += PDB_construct_seqs_aligned
    top_PDB_ID = PDB_construct_seqs_aligned[0][0]

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
    alignment_html = generate_html_from_alignment(target, sorted_alignment_seqs, sorted_alignment_IDs, additional_data_fields=additional_data, aa_css_class_list=aa_css_class_list_UniProt_seq)

    # add the alignment html to the main html tree
    for child in alignment_html.getchildren():
        output_html_body.append(child)

    # ===========
    # write to html file
    # ===========
    ofilename = os.path.join(html_alignments_dir, target + '.html')
    write_output(ofilename, output_html_tree)

    # Add top PDB ID to targets_data dicts
    targets_data[t][target].append(top_PDB_ID)

# ===========
# print and write file containing sorted targets with details on PDB structure selection
# ===========

print ''

PDB_selections_text = ''
PDB_selections_text += '= Targets sorted by number of PDB structures with matching expression system tag =\n\n'
PDB_selections_text += '%18s  %24s  %16s\n' % ('target', 'nmatching_PDB_structures', 'top_PDB_chain_ID')
PDB_selections_text += '%18s  %24s  %16s\n' % ('______', '________________________', '________________')

for target_dict in targets_data:
    target = target_dict.keys()[0]
    # For some reason, target_dict.values()[0][0] (nmatching_PDB_structures) needs to be referenced directly - assigning to a variable beforehand results in a value of 0. No idea why.
    if target_dict.values()[0][0] == 0 or len(target_dict.values()[0]) < 2:
        continue
    top_PDB_ID = target_dict.values()[0][1]
    PDB_selections_text += '%18s  %24d  %16s\n' % (target, target_dict.values()[0][0], top_PDB_ID)

PDB_selections_text += '\nTotal targets with > 0 matching PDB structures: %d\n' % sum([1 for target_dict in targets_data if target_dict.values()[0][0] > 0])
PDB_selections_text += 'Total targets with > 1 matching PDB structures: %d\n' % sum([1 for target_dict in targets_data if target_dict.values()[0][0] > 1])

print PDB_selections_text

with open(output_selections_filepath, 'w') as output_selections_file:
    output_selections_file.write(PDB_selections_text)

