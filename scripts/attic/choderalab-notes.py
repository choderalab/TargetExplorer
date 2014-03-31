# Various notes on kinases. Produces choderalab-notes.xml
#
# Daniel L. Parton <partond@mskcc.org> - 3 May 2013
#

import os
from lxml import etree

ofile_path = os.path.join('database', 'choderalab-notes.xml')

root = etree.Element('kinome-notes')

# ===================================
# Put notes here
# ===================================

kinase = etree.Element('kinase')
kinase.set('uniprotAC','Q96GX5')
kinase.set('entry_name','GWL_HUMAN')
etree.SubElement(kinase,'note').text = '''
Long insert of ~540 aa, toward the C-term of the activation loop (checked via comparison with 3CQU structure of AKT1).

From __REF: Blake-Hodex Goldberg Mol Cell Biol 2012 10.1128/MCB.06525-11__
* The kinase Gwl appears to be essential for M phase entry and maintenance in drosophila, xenopus laevis, and human tissue culture cells
* The authors refer to the ~500 aa insert as NCMR, for nonconserved middle region
* The NCMR is poorly conserved between species
* Appears to be largely dispensible for Gwl functionality, at least in vitro
* However, many mitotic phosphorylation sites are found throughout the NCMR

PSIPRED indicates the NCMR is mostly unstructured, with a few short structured regions.

Looks like NCMR is truly a part of the transcribed protein (sequencing error seemed a possibility), but its function remains clear.
'''
root.append(kinase)

# ===================================

kinase = etree.Element('kinase')
kinase.set('uniprotAC','Q96GX5')
kinase.set('entry_name','GWL_HUMAN')
etree.SubElement(kinase,'note').text = '''
Long insert of ~540 aa, toward the C-term of the activation loop (checked via comparison with 3CQU structure of AKT1).

From __REF: Blake-Hodex Goldberg Mol Cell Biol 2012 10.1128/MCB.06525-11__
* The kinase Gwl appears to be essential for M phase entry and maintenance in drosophila, xenopus laevis, and human tissue culture cells
* The authors refer to the ~500 aa insert as NCMR, for nonconserved middle region
* The NCMR is poorly conserved between species
* Appears to be largely dispensible for Gwl functionality, at least in vitro
* However, many mitotic phosphorylation sites are found throughout the NCMR

PSIPRED indicates the NCMR is mostly unstructured, with a few short structured regions.

Looks like NCMR is truly a part of the transcribed protein (sequencing error seemed a possibility), but its function remains clear.
'''
root.append(kinase)

# ===================================
# Output the .xml document
# ===================================

ofile = open(ofile_path, 'w')
ofile.write( etree.tostring(root, pretty_print=True) )
ofile.close()

