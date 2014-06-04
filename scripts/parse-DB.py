#!/usr/bin/env python
#
# Simple script to print kinDB_ids of kinases returned by a given xpath query
#
# Daniel L. Parton <partond@mskcc.org> - 7 May 2013
#

import os
from lxml import etree
import TargetExplorer as clab

analysis_dir = 'analysis'
DB_version = 'database.xml'
DB_filepath = os.path.join('database', DB_version)
DB_root = etree.parse(DB_filepath).getroot()
nentries = len(DB_root)
nhuman_entries = len(DB_root.findall('entry/UniProt[@NCBI_taxID="9606"]'))
print 'Number of entries:', nentries
print 'Number of human entries:', nhuman_entries

# XXX
npdb_structures = len( DB_root.findall('entry/UniProt/../PDB/structure') )
npdb_structures_hum = len( DB_root.findall('entry/UniProt[@NCBI_taxID="9606"]/../PDB/structure') )
print 'Number of PDB structures:', npdb_structures
print 'Number of human PDB structures:', npdb_structures_hum

npdb=0
npdbhum=0
nscop=0
ndisease=0
npdb_disease = 0
nscop_disease = 0
npdb_no_scop = 0
nscop_no_pdb = 0
ntargets = 0
nfamily_kinases=0
for i in range(nentries):
    taxID = DB_root[i].find('UniProt').get('NCBI_taxID')
    AC = DB_root[i].find('UniProt').get('AC')
    entry_name = DB_root[i].find('UniProt').get('entry_name')
    pdb = DB_root[i].findall('PDB/structure')
    #scop = DB_root[i].findall('scop/scop_domain')
    disease = DB_root[i].findall('UniProt/disease_associations/disease_association')
    target_domains = DB_root[i].findall('UniProt/domains/domain')
    ntargets += len(target_domains)
    family = DB_root[i].find('UniProt').get('family')
    if pdb != []:
        npdb+=1
        if taxID == '9606':
            npdbhum+=1
    if disease != [] and taxID == '9606':
        ndisease+=1
    if (pdb != []) and (disease != []) and taxID == '9606':
        npdb_disease+=1
    #if (scop != []) & (disease != []):
    #    nscop_disease+=1
    #if (len(scop) == 0) & (len(pdb) != 0):
    #    npdb_no_scop+=1
    #if (len(pdb) == 0) & (len(scop) != 0):
    #    nscop_no_pdb+=1
    if family in ['AGC','CAMK','CMGC','CK1','STE','TKL','TK'] and taxID == '9606':
        nfamily_kinases+=1
    #    print entry_name

print 'Number of human kinases belonging to one of the major families:', nfamily_kinases
print 'Number of human target domains:', ntargets
#print 'Number of kinases with a disease association:', ndisease
#print 'Number of kinases with a SCOP PK domain:', nscop
#print 'Number of kinases with both a SCOP PK domain and a disease association:', nscop_disease
print 'Number of kinases with a PDB structure which includes the PK domain:', npdb
print 'Number of human kinases with a PDB structure which includes the PK domain:', npdbhum
#print 'Number of kinases with both a PK domain PDB structure and a disease association:', npdb_disease
#print 'Number of kinases with a PK domain PDB structure and no SCOP entry:', npdb_no_scop
#print 'Number of kinases with a SCOP entry and no PK domain PDB structure:', nscop_no_pdb

nmutants = len(DB_root.findall('entry/UniProt[@NCBI_taxID="9606"]/../mutants/mutant'))
nbioassays = len(DB_root.findall('entry/UniProt[@NCBI_taxID="9606"]/../bioassays/bioassay'))
npubs = len(DB_root.findall('entry/UniProt[@NCBI_taxID="9606"]/../NCBI_Gene/entry/publications/publication'))
print 'Total number of mutants (human):', nmutants
print 'Total number of bioassays (human):', nbioassays
print 'Total number of publications (human):', npubs

print ''

# ==========================
# For a given type of data node, print all values of a given attrib
# ==========================
def print_attribs(query_data_node, query_attrib):
    data_nodes = DB_root.findall(query_data_node)
    attrib_values = []
    for data_node in data_nodes:
        attrib_values.append(data_node.get(query_attrib))
    attrib_values = set(attrib_values)
    for attrib_value in attrib_values:
        print attrib_value

# ==========================
# Print pk_pdb experimental_aln_conflicts for pk_domains > 350 in span
# ==========================

# print '\n= Length of PK domains annotated in UniProt =\n'
# for k in range(nentries):
#     length = int( DB_root[k].find('uniprot/pk_domain').get('length') )
#     print length
#     if length > 350:
#         entry_name = DB_root[k].find('uniprot').get('entry_name')
#         AC = DB_root[k].find('uniprot').get('AC')
#         print entry_name, AC
#         print DB_root[k].findtext('uniprot/function')
#         pk_pdbs = DB_root[k].findall('pk_pdb')
#         for p in pk_pdbs:
#             pdbid = p.get('id')
#             chains = p.findall('chain')
#             for c in chains:
#                 chainid = c.get('id')
#                 exp_seq = c.findtext('experimental/sequence_aln_conflicts')
#                 print pdbid, chainid
#                 print ''.join(exp_seq.split('\n'))


# ==========================
# Print pk_pdb experimental sequences in fasta format
# ==========================

# for k in range(nentries):
#     entry_name = DB_root[k].find('uniprot').get('entry_name')
#     pk_pdbs = DB_root[k].findall('pk_pdb')
#     pk_domain = DB_root[k].find('uniprot/pk_domain')
#     pk_begin = int(pk_domain.get('begin'))
#     pk_end = int(pk_domain.get('end'))
#     for p in pk_pdbs:
#         chains = p.findall('chain')
#         for c in chains:
#             pdbid = p.get('id')
#             whole_exp_seq = ''.join(c.findtext('experimental/sequence_aln_conflicts').split('\n'))
#             pk_exp_seq = whole_exp_seq[pk_begin-1 : pk_end]
#             print '>' + entry_name + '_' + pdbid
#             print pk_exp_seq + '\n'

# ==========================
# Print kinDB_ids in order of target_priority
# ==========================

def write_targets_prioritized():
    ofilepath = os.path.join(analysis_dir, 'targets-prioritized.txt')

    print 'Writing target priorities to %s' % ofilepath

    targets = []
    for k in range(nentries):
        DB_entry = DB_root[k]
        taxid = DB_entry.find('UniProt').get('NCBI_taxID')
        if taxid != '9606':
            continue
        target_domains = DB_entry.findall('UniProt/domains/domain[@targetID]')
        target_score_node = DB_entry.find('target_score')
        pubs_score = int(float(target_score_node.get('publications')))
        cbioportal_mutations_score = target_score_node.get('cBioPortal_mutations')
        disease_assocs = len(DB_entry.findall('UniProt/disease_associations/disease_association'))

        nbioassays = len(DB_entry.findall('bioassays/bioassay'))

        for target_domain in target_domains:    
            targetID = target_domain.get('targetID')
            target_domain_score_node = DB_entry.find('target_score/domain[@targetID="%s"]' % targetID)
            target_score = float( target_domain_score_node.get('target_score') )
            family = DB_entry.find('UniProt').get('family')
            target_domain_length = target_domain.get('length')
            nPDBs = target_domain_score_node.get('nPDBs')
            pseudogene_score = float(target_domain_score_node.get('pseudogene'))
            if pseudogene_score == 0.:
                pseudogene = ''
            else:
                pseudogene = 'pseudo'

            # Add info the target_info dict
            target_info = {}
            target_info['targetID'] = targetID
            target_info['target_score'] = target_score
            target_info['family'] = family
            target_info['target_domain_length'] = target_domain_length
            target_info['pubs'] = pubs_score
            target_info['muts'] = cbioportal_mutations_score
            target_info['disease'] = disease_assocs
            target_info['nbioassays'] = nbioassays
            target_info['nPDBs'] = nPDBs
            target_info['pseudogene'] = pseudogene
            targets.append(target_info)

    sorted_targets = sorted(targets, key=lambda x: x['target_score'], reverse=True)
    with open(ofilepath, 'w') as ofile:
        ofile.write( '{:<25}  {:<12}    {:<10}  {:<13}  {:<6}    {:<8}  {:<8}  {:<15}  {:<10}  {:<10}\n'.format('targetID', 'target_score', 'family', 'domain_length', 'npubs', '%muts', 'nPDBs', 'ndisease_assocs', 'nbioassays', 'pseudogene') )
        for t, target in enumerate(sorted_targets):
            ofile.write( '{0[targetID]:<25}  {0[target_score]: 12.1f}    {0[family]:<10}  {0[target_domain_length]:<13}  {0[pubs]:<6}    {0[muts]:<8}  {0[nPDBs]:<8}  {0[disease]:<15}  {0[nbioassays]:<10}  {0[pseudogene]:<10}\n'.format(target) )
            if t == 9:
                ofile.write( '-'*130 + '\n' )

    print 'Done.'

# ==========================
# Print kinases with pk_pdb entries, but no scop entries
# ==========================

#for k in range(nentries):
#    kinase = DB_root[k]
#    name = kinase.find('uniprot').get('entry_name')
#    pk_pdbs = kinase.findall('pk_pdb')
#    scops = kinase.findall('scop/scop_domain')
#    #if len(pk_pdbs) > 0 and len(scops) == 0:
#    #    print name
#    if len(scops) > 0:
#        print name

# ==========================
# Print total number of pk_pdb entries
# ==========================

#n_pk_pdb = len(DB_root.findall('kinase/pk_pdb'))
#n_chains = len(DB_root.findall('kinase/pk_pdb/chain'))
#print 'Total number of pk_pdbs', n_pk_pdb
#print 'Total number of pk_pdb chains', n_chains

# ==========================
# Print NCBI Gene IDs
# ==========================

def print_GeneIDs():
    for k in range(nentries):
        kinase = DB_root[k]
        entry_name = kinase.find('uniprot').get('entry_name')
        print entry_name
        IDs = [ id.get('ID') for id in kinase.findall('NCBI_Gene/entry') ]
        if len(IDs) != 1:
            print '####                        ', IDs

# ==========================
# Print Ensembl gene IDs
# ==========================

#for k in range(nentries):
#    kinase = DB_root[k]
#    entry_name = kinase.find('uniprot').get('entry_name')
#    print entry_name
#    ensemblgeneIDs = [ id.text for id in kinase.findall('EnsemblGeneID') ]
#    if len(ensemblgeneIDs) != 1:
#        print '####                        ', ensemblgeneIDs

# ==========================
# Print HGNC gene Approved Symbols
# ==========================

#for k in range(nentries):
#    kinase = DB_root[k]
#    entry_name = kinase.find('uniprot').get('entry_name')
#    print entry_name.split('_')[0],
#    HGNC_entries = kinase.findall('HGNC/entry')
#    for HGNC_entry in HGNC_entries:
#        HGNC_Approved_Symbol = HGNC_entry.get('Approved_Symbol')
#        print HGNC_Approved_Symbol
#    if len(HGNC_entries) != 1:
#        print '\n####                        ', HGNC_entries

# ==========================
# Look for gatekeeper mutants
# ==========================

#for entry_name in ['SRC_HUMAN', 'ABL1_HUMAN']:
#    kinase_node = DB_root.find('kinase/uniprot[@entry_name="%(entry_name)s"]/..' % vars())
#    if entry_name == 'SRC_HUMAN':
#        mutation = kinase_node.find('mutants/mutant/mutation[aa_change="T341I"]')
#    elif entry_name == 'ABL1_HUMAN':
#        mutation = kinase_node.find('mutants/mutant/mutation[aa_change="T315I"]')
#    print mutation

# ==========================
# Print number of publications retrieved from NCBI Gene for each kinase
# ==========================

def print_pubs():
    for k in range(nentries):
        kinase = DB_root[k]
        entry_name = kinase.find('uniprot').get('entry_name')
        print entry_name
        pubs = [ id.text for id in kinase.findall('NCBI_Gene/entry/publications/publication') ]
        print len(pubs)

# ==========================
# Print number of bioassays retrieved from BindingDB for each kinase
# ==========================

def print_bioassays(bioassay_type=None):
    for k in range(nentries):
        kinase = DB_root[k]
        entry_name = kinase.find('uniprot').get('entry_name')
        if bioassay_type == None:
            bioassays = kinase.findall('bioassays/bioassay')
        else:
            bioassays = kinase.xpath('bioassays/bioassay[@%s!="n/a"]' % bioassay_type)
            if len(bioassays) == 0:
                continue
        print entry_name
        print len(bioassays)

# ==========================
# Print example bioassays data for a given kinase
# ==========================

def print_bioassays_kinase(entry_name):
    kinase_node = DB_root.find('kinase/uniprot[@entry_name="%s"]/..' % entry_name)
    bioassays = kinase_node.findall('bioassays/bioassay[@ligand_BindingDB_ID="50187566"]')
    for bioassay in bioassays:
        print bioassay.attrib

# ==========================
# Print total number of BindingDB bioassays of each type (Kd, IC50 etc.)
# ==========================

def print_bioassay_types(bioassay_type=None):
    nKd = len( DB_root.findall('kinase/bioassays/bioassay[@Kd]') )
    nKi = len( DB_root.findall('kinase/bioassays/bioassay[@Ki]') )
    nIC50 = len( DB_root.findall('kinase/bioassays/bioassay[@IC50]') )
    nEC50 = len( DB_root.findall('kinase/bioassays/bioassay[@EC50]') )
    nkon = len( DB_root.findall('kinase/bioassays/bioassay[@kon]') )
    nkoff = len( DB_root.findall('kinase/bioassays/bioassay[@koff]') )
    print '{:<10} {:<8}'.format('Kd:', nKd)
    print '{:<10} {:<8}'.format('Ki:', nKi)
    print '{:<10} {:<8}'.format('IC50:', nIC50)
    print '{:<10} {:<8}'.format('EC50:', nEC50)
    print '{:<10} {:<8}'.format('kon:', nkon)
    print '{:<10} {:<8}'.format('koff:', nkoff)

# ==========================
# Print sequences recorded in BindingDB and compare with UniProt sequence
# ==========================

def print_bioassay_seqs(kinase_entry_names='all'):
    if kinase_entry_names == 'all':
        kinase_nodes = DB_root.findall('kinase')
    else:
        kinase_nodes = []
        for entry_name in kinase_entry_names:
            kinase_nodes.append( DB_root.find('kinase/uniprot[@entry_name="%s"]/..' % entry_name) )
    for kinase_node in kinase_nodes:
        entry_name = kinase_node.find('uniprot').get('entry_name')
        uniprot_seq = clab.core.sequnwrap( kinase_node.findtext('uniprot/sequence') )
        #bioassays = kinase_node.findall('bioassays/bioassay[@BindingDB_source="ChEMBL"][@ligand_ChEMBL_ID="CHEMBL1221411"]')
        bioassays = kinase_node.findall('bioassays/bioassay')
        for bioassay in bioassays:
            try:
                bindingdbseq = bioassay.get('target_sequence')
                if bindingdbseq == 'n/a':
                    continue
                clustalo_result = clab.align.run_clustalo(['BindingDB', entry_name], [bindingdbseq, uniprot_seq]) 
                print entry_name, bioassay.get('ligand_BindingDB_ID')
                comparison = clab.align.seq_comparator(clustalo_result[0], clustalo_result[1])
                if comparison != '*' * len(comparison):
                    print bindingdbseq
                    print comparison
            except Exception as e:
                print bindingdbseq
                raise e

# ==========================
# Print kinases with PDB structures which were expressed in a given expression system
# ==========================

def print_kinases_expressed_in(expression_system):
    all_target_nodes = DB_root.findall('entry/UniProt/domains/domain[@targetID]')
    target_nodes = []
    for t in all_target_nodes:
        target_domain_id = t.get('domainID')
        expression_data = t.xpath('../../../PDB/structure/chain[@domainID="%s"]/../expression_data[contains(@EXPRESSION_SYSTEM,"%s")]' % (target_domain_id, expression_system))
        if expression_data != []:
            target_nodes.append(t)

    #target_nodes = DB_root.xpath('kinase/uniprot/pk_domain/../../pk_pdb/expression_data[contains(@EXPRESSION_SYSTEM,"%s")]' % expression_system)

    print 'Number of PK domains with E coli expression host listed in PDB:' , len(target_nodes)
    target_info_dicts = []
    for target_node in target_nodes:
        target_info_dict = clab.DB.target_info_dict(target_node)
        target_info_dicts.append( target_info_dict )
        example_pdbs = [ x.get('ID') for x in target_node.xpath('../../../PDB/structure/chain[@domain_ID="%d"]/../expression_data[contains(@EXPRESSION_SYSTEM,"%s")]/..' % (target_info_dict['pk_domain_id'], expression_system)) ]
        example_pdbs_string = '%6d' % len(example_pdbs) + '  ' + '  '.join(example_pdbs[0:10])
        target_info_dicts[-1]['example_pdbs_string'] = example_pdbs_string

    #print len(target_info)
    #print len(target_info[0])

    # Sort
    print clab.DB.format_target_info(return_format_string=True) + '      npdbs with expression_system "*%s*", and first 10 PDB codes' % expression_system
    sorted_target_info = sorted(target_info_dicts, key=lambda x: x['target_score'], reverse=True)
    for t in sorted_target_info:
        print clab.DB.format_target_info(t), t['example_pdbs_string']

# ==========================
# Print number of kinases with PDB structures which were expressed in a given expression system
# ==========================

def print_nkinases_expressed_in(expression_systems):
    expr_nodes = DB_root.xpath('entry/UniProt[@family="TK"]/../PDB/structure/expression_data[@EXPRESSION_SYSTEM]')
    expr_systems = [node.get('EXPRESSION_SYSTEM') for node in expr_nodes]
    set_expr_systems = set(expr_systems)

    for expression_system in set_expr_systems:
        expression_data = DB_root.xpath('entry/UniProt[@family="TK"]/../PDB/structure/expression_data[contains(@EXPRESSION_SYSTEM,"%s")]' % (expression_system))
        print expression_system, len(expression_data)

    # for expression_system in expression_systems:

# ==========================
# Print non-human PDBs expressed in
# ==========================

def print_nonhuman_pdbs_expressed_in(expression_system):
    PDBs = DB_root.xpath('entry/UniProt[@NCBI_taxID!="9606"]/../PDB/structure/expression_data[contains(@EXPRESSION_SYSTEM,"%s")]/..' % expression_system)
   #  PDBs = DB_root.xpath('entry/PDB/structure/expression_data[contains(@EXPRESSION_SYSTEM,"%s")]' % expression_system)
    pdbids = ['%20s  %s' % (pdb.getparent().getparent().find('UniProt').get('entry_name'), pdb.get('ID')) for pdb in PDBs]
    print '\n'.join(pdbids)

# ==========================
# Prints PDB info for a given DB entry
# ==========================

def print_pdbs_by_UniProt_entry_name(entry_name, print_expr_data=False):
    DB_entry = DB_root.find('entry/UniProt[@entry_name="%s"]' % entry_name).getparent()
    pdb_structures = DB_entry.findall('PDB/structure')
    pdb_chains = DB_entry.findall('PDB/structure/chain')
    print 'Number of PDB structures:', len(pdb_structures)
    print 'Number of PDB chains:', len(pdb_chains)
    print 'PDB structures:'
    for pdb_structure in pdb_structures:
        PDB_ID = pdb_structure.get('ID')
        print_str = PDB_ID

        if print_expr_data != False:
            expr_data = pdb_structure.find('expression_data')
            if expr_data != None:
                if print_expr_data == 'all':
                    for key in expr_data.keys():
                        print_str += '  ' + key + '=' + expr_data.get(key)
                else:
                    print_str += '  ' + print_expr_data + '=' + expr_data.get(print_expr_data)

        print print_str

# ==========================
# Put function calls here
# ==========================
#print_attribs('entry/PDB/structure/expression_data', 'EXPRESSION_SYSTEM')
#write_targets_prioritized()
#print_GeneIDs()
#print_pubs()
#print_bioassays_kinase('ABL1_HUMAN')
#print_bioassay_types()
#print_bioassay_seqs(['SRC_HUMAN'])
#print_kinases_expressed_in('ESCHERICHIA')
#print_nkinases_expressed_in(['ESCHERICHIA', 'SPODOPTERA'])
print_nonhuman_pdbs_expressed_in('ESCHERICHIA')
#print_pdbs_by_UniProt_entry_name('ABL1_HUMAN', print_expr_data='EXPRESSION_SYSTEM')

