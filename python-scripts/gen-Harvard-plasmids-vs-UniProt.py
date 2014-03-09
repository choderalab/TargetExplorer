import os, shutil

main_branch_path = '..'
html_page_path_main_branch = os.path.join(main_branch_path, 'external-data', 'plasmids', 'DFHCC-PlasmID', 'HIP_human_kinase_collection-pJP1520', 'aln-against-UniProt-seq.html')

shutil.copy(html_page_path_main_branch, 'Harvard-plasmids-vs-UniProt.html')

with open('Harvard-plasmids-vs-UniProt.html', 'r') as aln_file:
    aln_text = aln_file.read()

with open('Harvard-plasmids-vs-UniProt.html', 'w') as aln_file:
    aln_file.write(aln_text.replace('<link type="text/css" href="seqlib.cs" rel="stylesheet"/>', '<link type="text/css" href="stylesheets/seqlib.css" rel="stylesheet"/>'))


