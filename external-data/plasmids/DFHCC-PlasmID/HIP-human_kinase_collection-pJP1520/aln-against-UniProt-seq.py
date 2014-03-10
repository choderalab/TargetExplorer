import sys, os, re
import openpyxl
from openpyxl import Workbook
import Bio, Bio.Seq, Bio.Alphabet
from lxml import etree
from lxml.builder import E
import TargetExplorer as clab
import pandas as pd

DB_path = os.path.join('..', '..', '..', '..', 'database', 'database.xml')
DB_root = etree.parse(DB_path).getroot()

plasmid_df = pd.read_csv('plasmid_insert_data.csv')


DB_targets = DB_root.findall('entry/UniProt/domains/domain[@targetID]')

# To be used to construct a pandas DataFrame
output_data = {
'targetID':[],
'cloneID':[],
'target_seq_aln':[],
'plasmid_seq_aln':[],
'nconflicts':[],
'len_seq_aln':[],
'pctidentity':[]
}

# ===========
# Iterate through plasmids
# ===========

for p in plasmid_df.index:
    cloneID = plasmid_df['cloneID'][p]
    UniProtAC = plasmid_df['UniProtAC'][p]
    UniProt_entry_name = plasmid_df['UniProt_entry_name'][p]
    #UniProt_canonseq = plasmid_df['UniProt_canonseq'][p]
    insert_aa_seq = plasmid_df['insert_aa_seq'][p]

    if UniProtAC == 'None':
        continue

    DB_UniProt_node = DB_root.find('entry/UniProt[@AC="%s"]' % UniProtAC)
    DB_entry = DB_UniProt_node.getparent()

    DB_domains = DB_UniProt_node.findall('domains/domain[@targetID]')

    # ===========
    # Iterate through target domains
    # ===========

    for d in range(len(DB_domains)):

        targetID = DB_domains[d].get('targetID')
        domain_seq = clab.core.sequnwrap( DB_domains[d].find('sequence').text )

        # =====
        # Run alignment
        # =====

        alnIDs = [targetID, cloneID]
        pre_aln_seqs = [domain_seq, insert_aa_seq]

        aligned_seqs = clab.align.run_clustalo(alnIDs, pre_aln_seqs)
        alignment = [ [alnIDs[i], seq] for i, seq in enumerate(aligned_seqs) ]

        # ===========
        # convert conflicts in plasmid seq to lower case and find the target domain region in the alignment coords
        # ===========
        seq = list(alignment[1][1]) # sequence needs to be mutable so convert to list
        for aa_iter in range(len(seq)):
            if seq[aa_iter] != alignment[0][1][aa_iter]:
                seq[aa_iter] = seq[aa_iter].lower()
        alignment[1][1] = ''.join(seq)

        regex_start = re.search('[A-Z]', alignment[0][1])
        regex_end = re.search('[A-Z]', alignment[0][1][::-1])
        target_domain_start_aln_coords = regex_start.start()
        target_domain_end_aln_coords = len(alignment[0][1]) - regex_end.start() - 1  # 0-based index of the last residue

        # ===========
        # calculate number of conflicts and percentage
        # ===========
        aln_seqs_target_domain_region = [aln_seq[1][target_domain_start_aln_coords : target_domain_end_aln_coords + 1] for aln_seq in alignment]
        nconflicts = len( [ aa_iter for aa_iter in range(len(aln_seqs_target_domain_region[0])) if aln_seqs_target_domain_region[0][aa_iter] == '-' or aln_seqs_target_domain_region[1][aa_iter] == '-' or re.match('[a-z]', aln_seqs_target_domain_region[1][aa_iter]) ] )
        pctidentity = (len(domain_seq) - nconflicts) * 100. / float(len(domain_seq))

        # ===========
        # append to data dict, to be used later to construct pandas DataFrame
        # ===========
        output_data['targetID'].append(targetID)
        output_data['cloneID'].append(cloneID)
        output_data['target_seq_aln'].append(alignment[0][1])
        output_data['plasmid_seq_aln'].append(alignment[1][1])
        output_data['nconflicts'].append(nconflicts)
        output_data['len_seq_aln'].append(str(len(domain_seq)))
        output_data['pctidentity'].append('%.2f' % pctidentity)

    #break

# construct pandas DataFrame and write to csv
output_data = pd.DataFrame(output_data)
#df.to_csv('aln-against-UniProt-seq.csv', index=False)

output_data = output_data.sort('nconflicts', ascending=True)
output_data.index = range(len(output_data)) # redo indices, as these will have been sorted in the previous step

# write data txt file
with open('aln.txt', 'w') as otxtfile:
    columnwidths = {}
    for key in output_data.keys():
        columnwidths[key] = max( [ len(str(x)) for x in output_data[key] ] )

    IDcolumnwidth = max( [columnwidths['targetID'], columnwidths['cloneID']] )

    for i in range(len(output_data)):
        otxtfile.write('%-*s' % (IDcolumnwidth, output_data['targetID'][i]))
        otxtfile.write('  ')
        otxtfile.write('%*s' % (columnwidths['pctidentity'], '') )
        otxtfile.write('  ')
        otxtfile.write('%*s' % (columnwidths['nconflicts'], '') )
        otxtfile.write(' ')
        otxtfile.write('%-*s' % (columnwidths['len_seq_aln'], '') )
        otxtfile.write('  ')
        otxtfile.write('%s' % output_data['target_seq_aln'][i])
        otxtfile.write('\n')

        otxtfile.write('%-*s' % (IDcolumnwidth, output_data['cloneID'][i]))
        otxtfile.write('  ')
        otxtfile.write('%*s' % (columnwidths['pctidentity'], output_data['pctidentity'][i]) )
        otxtfile.write('  ')
        otxtfile.write('%*s' % (columnwidths['nconflicts'], str(output_data['nconflicts'][i])) )
        otxtfile.write('/')
        otxtfile.write('%-*s' % (columnwidths['len_seq_aln'], output_data['len_seq_aln'][i]) )
        otxtfile.write('  ')
        otxtfile.write('%s' % output_data['plasmid_seq_aln'][i])
        otxtfile.write('\n\n')

#with open('aln.html', 'w') as ohtmlfile:
#    ohtmlfile.write( etree.tostring(output_html_tree, pretty_print=True) )

