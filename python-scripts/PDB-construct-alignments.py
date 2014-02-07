#
# Daniel L. Parton <partond@mskcc.org> - 6 Feb 2014
#

import os
from lxml import etree
from lxml.builder import E
import locallib

# ============
# parameters
# ============

output_html_index_filepath = 'kinase-PDB-constructs.html'

PDB_constructs_xml_filepath = os.path.join('..', 'analysis', 'PDB_construct_selection', 'PDB_constructs-data.xml')

parser = etree.XMLParser(remove_blank_text=True)

aln_dir = 'PDB-construct-alignments'

# ============
# definitions
# ============

class PDBCnstrct_Table(object):

    def __init__(self):
        self.table_node = E.table(CLASS='table')
        self.add_table_headers()

    def add_table_headers(self):
        row = E.tr()
        row.append( E.th ( E.div('index'), nowrap='' ) )
        row.append( E.th ( E.div('nreps'), nowrap='' ) )
        row.append( E.th ( E.div('targetID'), nowrap='' ) )
        row.append( E.th ( E.div('nPDBs'), nowrap='' ) )
        row.append( E.th ( E.div('top_PDB'), nowrap='' ) )
        row.append( E.th ( E.div('expr_tag'), nowrap='' ) )
        row.append( E.th ( E.div('expr_tag_score'), nowrap='' ) )
        row.append( E.th ( E.div('expr_sys'), nowrap='' ) )
        row.append( E.th ( E.div('DB_target_score'), nowrap='' ) )
        row.append( E.th ( E.div('DB_target_rank'), nowrap='' ) )
        row.append( E.th ( E.div('UniProt_AC'), nowrap='' ) )
        self.table_node.append(row)

    def add_construct_row(self, index, construct_xml):
        targetID = construct_xml.get('targetID')
        nreplicates = construct_xml.get('nreplicates')
        nPDBs = construct_xml.get('nmatching_PDB_structures')
        top_PDB = construct_xml.get('top_PDB_chain_ID')
        PDBID = top_PDB.split('_')[0]
        expr_tag = construct_xml.get('expr_tag_string')
        auth_score = construct_xml.get('authenticity_score')
        expr_sys = construct_xml.get('expression_system')
        targ_score = construct_xml.get('DB_target_score')
        targ_rank = construct_xml.get('DB_target_rank')
        AC = construct_xml.get('UniProt_AC')

        row = E.tr()

        row.append( E.td( E.div(str(index+1)), nowrap='' ) )
        row.append( E.td( E.div(nreplicates), nowrap='' ) )

        # targetID, href to alignment html
        targetID_html = E.a(targetID)
        targetID_html.set('href', os.path.join(aln_dir, targetID + '.html'))
        row.append( E.td( targetID_html, nowrap='' ) )

        row.append( E.td( E.div(nPDBs), nowrap='' ) )

        # top PDB, href to PDB page
        top_PDB_html = E.a(top_PDB)
        top_PDB_html.set('href', 'http://www.rcsb.org/pdb/explore.do?structureId=%s' % PDBID)
        row.append( E.td( top_PDB_html, nowrap='' ) )

        row.append( E.td( E.div(expr_tag), nowrap='' ) )
        row.append( E.td( E.div(auth_score), nowrap='' ) )
        row.append( E.td( E.div(expr_sys), nowrap='' ) )
        row.append( E.td( E.div(targ_score), nowrap='' ) )
        row.append( E.td( E.div(targ_rank), nowrap='' ) )

        # UniProt AC, href to UniProt page
        AC_html = E.a(AC)
        AC_html.set('href', 'http://www.uniprot.org/uniprot/%s' % AC)
        row.append( E.td( AC_html, nowrap='' ) )

        self.table_node.append(row)


# ============
# main
# ============

# read in XML

PDB_constructs_xml = etree.parse(PDB_constructs_xml_filepath, parser).getroot()

# Generate main content section node
main_content_node = locallib.gen_main_content_node()

# TODO add info on expr_tag_score etc.

# Title for main content section
h2_node = etree.SubElement(main_content_node, 'h2')
a_node = etree.SubElement(h2_node, 'a')
a_node.set('name', 'PDB construct alignments')
a_node.set('class', 'anchor')
a_node.set('href', '#PDBalns')
span_node = etree.SubElement(a_node, 'span')
span_node.set('class', 'octicon octicon-link')
h2_node.text = 'PDB construct alignments'

# Now construct the table of PDB construct data

pdbcnstrct_table = PDBCnstrct_Table()

for t, target_xml in enumerate(PDB_constructs_xml):
    pdbcnstrct_table.add_construct_row(t, target_xml)

#main_content_node.append(pdbcnstrct_table.table_node)
table_div = E.div()
table_div.append(pdbcnstrct_table.table_node)
#table_div.set('style', 'position:absolute;margin-left: -200px;')
main_content_node.append(table_div)





# TODO add footer info




webpage = locallib.Webpage()
webpage.add_standard_header()
webpage.add_main_content(main_content_node)
webpage.write_html(output_html_index_filepath)


