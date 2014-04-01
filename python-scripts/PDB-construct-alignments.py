#
# Daniel L. Parton <partond@mskcc.org> - 6 Feb 2014
#

import os, shutil
from lxml import etree
from lxml.builder import E
import locallib
import pandas as pd
import numpy as np

# ============
# parameters
# ============

output_html_index_filepath = 'kinase-PDB-constructs.html'

parser = etree.XMLParser(remove_blank_text=True)

aln_dir = 'PDB-construct-alignments'

manual_exceptions_yaml_filepath = os.path.join(aln_dir, 'manual_exceptions.yaml')
manual_exceptions_html_filepath = os.path.join(aln_dir, 'manual_exceptions.html')
manual_exceptions_yaml_filepath_main_branch = os.path.join('..', 'PDB_construct_selection', 'manual_exceptions.yaml')

PDB_constructs_csv_filepath = os.path.join(aln_dir, 'PDB_constructs-data.csv')
PDB_constructs_csv_filepath_main_branch = os.path.join('..', 'PDB_construct_selection', 'PDB_constructs-data.csv')

aln_dir_main_branch = os.path.join('..', 'PDB_construct_selection', 'alignments')

# ============
# definitions
# ============

class PDBCnstrct_Table(object):

    def __init__(self):
        self.table_node = E.table(CLASS='PDBCnstrctTable')
        self.add_table_headers()

    def add_table_headers(self):
        row = E.tr()
        row.append( E.th ( E.div('index'), nowrap='' ) )
        #row.append( E.th ( E.div('nreps'), nowrap='' ) )
        row.append( E.th ( E.div('targetID'), nowrap='' ) )
        row.append( E.th ( E.div('nPDBs'), nowrap='' ) )
        row.append( E.th ( E.div('top_PDB'), nowrap='' ) )
        row.append( E.th ( E.div('expr_tag'), nowrap='' ) )
        row.append( E.th ( E.div('expr_tag_score'), nowrap='' ) )
        row.append( E.th ( E.div('expr_sys'), nowrap='' ) )
        row.append( E.th ( E.div('construct_organism'), nowrap='' ) )
        row.append( E.th ( E.div('DB_target_score'), nowrap='' ) )
        row.append( E.th ( E.div('DB_target_rank'), nowrap='' ) )
        row.append( E.th ( E.div('UniProt_AC'), nowrap='' ) )
        self.table_node.append(row)

    def add_construct_row(self, index, construct_df):
        targetID = construct_df.get('targetID')
        #nreplicates = construct_df.get('nreplicates')
        nPDBs = str( construct_df.get('nmatching_PDB_structures') )
        top_PDB = construct_df.get('top_PDB_chain_ID')
        PDBID = top_PDB.split('_')[0]
        expr_tag = construct_df.get('top_cnstrct_expr_tag')
        auth_score = str( construct_df.get('top_cnstrct_auth_score') )
        expr_sys = construct_df.get('top_cnstrct_expr_sys')
        taxname = construct_df.get('top_cnstrct_taxname')
        targ_score = str( construct_df.get('DB_target_score') )
        targ_rank = str( int( construct_df.get('DB_target_rank') ) )
        AC = construct_df.get('UniProtAC')

        row = E.tr()

        row.append( E.td( E.div(str(index+1)), nowrap='' ) )
        #row.append( E.td( E.div(nreplicates), nowrap='' ) )

        # targetID, href to alignment html
        targetID_html = E.a(targetID)
        targetID_html.set('href', os.path.join(aln_dir, targetID + '.html'))
        row.append( E.td( targetID_html, nowrap='' ) )

        row.append( E.td( E.div(nPDBs), nowrap='' ) )

        # top PDB, href to PDB page
        top_PDB_html = E.a(top_PDB)
        top_PDB_html.set('href', 'http://www.rcsb.org/pdb/explore.do?structureId=%s' % PDBID)
        row.append( E.td( top_PDB_html, nowrap='' ) )

        try:
            row.append( E.td( E.div(expr_tag), nowrap='' ) )
        except TypeError:
            row.append( E.td( E.div('None'), nowrap='' ) )
        row.append( E.td( E.div(auth_score), nowrap='' ) )
        row.append( E.td( E.div(expr_sys), nowrap='' ) )
        row.append( E.td( E.div(taxname), nowrap='' ) )
        row.append( E.td( E.div(targ_score), nowrap='' ) )
        row.append( E.td( E.div(targ_rank), nowrap='' ) )

        # UniProt AC, href to UniProt page
        AC_html = E.a(AC)
        AC_html.set('href', 'http://www.uniprot.org/uniprot/%s' % AC)
        row.append( E.td( AC_html, nowrap='' ) )

        self.table_node.append(row)

def manual_exceptions_yaml2html(input_yaml_filepath, output_html_filepath):
    '''
    Just a quick/dirty solution for now.
    '''
    with open(output_html_filepath, 'w') as html_file:
        with open(input_yaml_filepath, 'r') as yaml_file:
            yaml_text = yaml_file.read()
            body_node = E.body(yaml_text, style='white-space: pre-wrap;')
            html_file.write(etree.tostring(body_node, pretty_print=True))


# ============
# main
# ============

# copy PDB_constructs-data.csv from main branch
shutil.copy(PDB_constructs_csv_filepath_main_branch, PDB_constructs_csv_filepath)

# copy HTML alignment files and seqlib.css from main branch
shutil.copy(manual_exceptions_yaml_filepath_main_branch, manual_exceptions_yaml_filepath)

# copy HTML alignment files and seqlib.css from main branch
for filename in os.walk(aln_dir_main_branch).next()[2]:
    if filename[0] != '.':
        mainbranch_filepath = os.path.join(aln_dir_main_branch, filename)
        ghpages_filepath = os.path.join(aln_dir, filename)
        shutil.copy(mainbranch_filepath, ghpages_filepath)

# read in CSV

PDB_constructs_df = pd.DataFrame.from_csv(PDB_constructs_csv_filepath)

# Generate main section section node
PDBCnstrct_section = locallib.gen_PDBCnstrct_section_node()

# Title for main section
h2_node = etree.SubElement(PDBCnstrct_section, 'h2')
h2_node.text = 'Selection of kinase constructs from PDB data'
a_node = etree.SubElement(h2_node, 'a')
a_node.set('name', 'PDB construct alignments')
a_node.set('class', 'anchor')
a_node.set('href', '#PDBalns')
span_node = etree.SubElement(a_node, 'span')
span_node.set('class', 'octicon octicon-link')

# Information subsection
h3_node = etree.SubElement(PDBCnstrct_section, 'h3')
h3_node.text = 'Information'

p_node = etree.SubElement(PDBCnstrct_section, 'p')
p_node.text = 'The aim is to select 96 kinase catalytic domain constructs suitable for E coli expression, using the PDB as a source of construct sequences which have been expressed successfully. The results are shown in the table below. '
p_node.set('style', 'max-width: 1000px; text-align: justify;')

#bullet = E.li( E.strong('nreps: '), 'Number of replicates chosen for expression testing' )
#PDBCnstrct_section.append(bullet)
bullet = E.li( E.strong('targetID: '), '[UniProt entry name]_[domain index], e.g. KS6A1 (aka RPS6KA1) has two kinase catalytic domains, which would be named KS6A1_HUMAN_D0 and KS6A1_HUMAN_D1. Click the targetID to see the PDB construct alignments' )
PDBCnstrct_section.append(bullet)
bullet = E.li( E.strong('nPDBs: '), 'Number of PDBs which list E. coli as the expression system' )
PDBCnstrct_section.append(bullet)
bullet = E.li( E.strong('top_PDB: '), 'Top-ranked PDB construct' )
PDBCnstrct_section.append(bullet)
bullet = E.li( E.strong('expr_tag: '), 'Type of expression tag detected from the construct sequence' )
PDBCnstrct_section.append(bullet)
bullet = E.li( E.strong('expr_tag_score: '), 'Arbitrary score assigned to the expression tag type. TEV-cleaved N-terminal tags are given the highest score' )
PDBCnstrct_section.append(bullet)
bullet = E.li( E.strong('expr_sys: '), 'Expression system listed for the PDB construct' )
PDBCnstrct_section.append(bullet)
bullet = E.li( E.strong('construct_organism: '), 'Organism source of the PDB construct' )
PDBCnstrct_section.append(bullet)
bullet = E.li( E.strong('DB_target_score: '), 'Score of kinase (subjective) importance based on metrics such as: number of disease associations; number of publications; percent mutations in a given cohort of sequenced tumor samples; number of bioassays listed in BindingDB' )
PDBCnstrct_section.append(bullet)
bullet = E.li( E.strong('DB_target_rank: '), 'Ranking based on DB_target_score' )
PDBCnstrct_section.append(bullet)
bullet = E.li( E.strong('UniProt_AC: '), 'Uniprot accession code' )
PDBCnstrct_section.append(bullet)
etree.SubElement(PDBCnstrct_section, 'br')

p_node = etree.SubElement(PDBCnstrct_section, 'p')
p_node.text = 'For a given kinase target, all PDB constructs expressed in E coli are run through a multiple sequence alignment, along with the UniProt canonical isoform sequence, and the sequence from a plasmid library (we will be subcloning our constructs from this library). The alignments can be viewed by clicking on the targetIDs in the table below (the kinase catalytic domain is highlighted in red in the UniProt sequence). The PDB constructs are ranked firstly based on their expr_tag_score, which denotes whether an expression tag has been detected in the sequence. The reason for this is to avoid PDB entries with erroneously annotated SEQRES records, which often record the structurally resolved residues rather than the input sequence (as they should). We consider construct sequences which include expression tags more likely to be authentic. After expr_tag_score, the construct sequences are ranked by alignment score (to avoid constructs with mutations), and then again by their sequence length compared to the kinase catalytic domain (since we only want to express this isolated domain).'
p_node.set('style', 'max-width: 1000px; text-align: justify;')

p_node = etree.SubElement(PDBCnstrct_section, 'p')
#p_node.text = 'Kinase targets are ranked firstly by expr_tag_score for the top PDB construct, then by nPDBs. Only constructs for which an expression tag has been detected (expr_tag_score > 0) are selected for expression testing. Currently this process results in 69 unique kinase constructs. For the remaining 27 well positions, replicate experiments will be run for the kinases with the highest DB_target_score. Running replicate experiments will help to allay cloning errors (error rate ~70%).'
p_node.text = 'Kinase targets are ranked firstly by expr_tag_score for the top PDB construct, then by nPDBs. Only constructs for which an expression tag has been detected (expr_tag_score > 0) are selected for expression testing.'
p_node.set('style', 'max-width: 1000px; text-align: justify;')

p_node = etree.SubElement(PDBCnstrct_section, 'p')
p_node = E.p('A series of ', E.a('manual exceptions', href=manual_exceptions_html_filepath), ' have also been made to further refine the results, which can be badly affected by misannotated PDB entries.')
p_node.set('style', 'max-width: 1000px; text-align: justify;')
PDBCnstrct_section.append(p_node)

# Now construct the table of PDB construct data
pdbcnstrct_table = PDBCnstrct_Table()

#for t, target_xml in enumerate(PDB_constructs_xml):
for t in range(len(PDB_constructs_df)):
    PDB_constructs_row = PDB_constructs_df.loc[t]
    if PDB_constructs_row['nmatching_PDB_structures'] > 0 and PDB_constructs_row['top_cnstrct_auth_score'] >= 0:
        pdbcnstrct_table.add_construct_row(t, PDB_constructs_row)

table_div = E.div()
table_div.append(pdbcnstrct_table.table_node)
PDBCnstrct_section.append(table_div)





webpage = locallib.PDBCnstrct_Webpage()
webpage.add_standard_header()
webpage.add_content(PDBCnstrct_section)
webpage.write_html(output_html_index_filepath)

manual_exceptions_yaml2html(manual_exceptions_yaml_filepath, manual_exceptions_html_filepath)


