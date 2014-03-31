# Assign scores to each PK domain in kinDB, based on various properties such as
# number of PDB structures, references, disease associations etc.
#
# Daniel L. Parton <partond@mskcc.org> - 29 April 2013
#

import sys, os
from lxml import etree
from numpy import *
import choderalab as clab

# ==============
# Parameters
# ==============

kinDB_path = os.path.join('database', 'kinDB-pdb-scop-mut-ncbi-bindingdb.xml')
ofile = os.path.join('database', 'kinDB-complete.xml')

parser = etree.XMLParser(remove_blank_text=True)
kinDB = etree.parse(kinDB_path, parser).getroot()
nkinases = len(kinDB)
print 'Number of kinases:', nkinases

# ==============
# Scoring functions
# ==============

def publications_score_fn(kinase_node):
    '''
    Scoring function based on the number of PubMed publications retrieved from NCBI Gene.
    '''
    pubs_node = kinase_node.find('NCBI_Gene/entry/publications')
    if pubs_node == None:
        return 0.
    else:
        return float(len(pubs_node))

def cbioportal_mutations_score_fn(kinase_node):
    '''
    Scoring function based on the number of mutations for a given kinase/mutants node
    '''
    percent_cases_with_mutations = clab.cBioPortal.percent_cases_with_mutations(kinase_node)
    return percent_cases_with_mutations

def disease_score_fn(kinase_node):
    '''
    Scoring function based on the number of disease associations for a certain kinase.
    '''
    ndisease_assocs = len( kinase_node.findall('uniprot/disease_association') )
    score = log(ndisease_assocs + 1)
    return score

def bioassays_score_fn(kinase_node):
    '''
    Scoring function based on the number of disease associations for a certain kinase.
    '''
    nbioassays = len( kinase_node.findall('bioassays/bioassay') )
    score = log(nbioassays + 1)
    return score

def pk_pdb_score_fn(x):
    '''
    Scoring function based on the number of PDB entries for a certain PK domain.
    '''
    y = log(x+1)
    return y

def pseudogene_score_fn(kinDB_id):
    '''
    Deprioritizes certain genes thought to be pseudogenes.

    PRKY_HUMAN_O43930: see comments on UniProt and Ensembl page
    '''    
    if kinDB_id in ['PRKY_HUMAN_O43930_PK0', 'PDPK2_HUMAN_Q6A1A2_PK0']:
        return -100.
    else:
        return 0.

def pk_length_score_fn(pk_domain):
    '''
    Deprioritizes certain genes thought to be pseudogenes.

    PRKY_HUMAN_O43930: see comments on UniProt and Ensembl page
    '''    
    pk_length_score = 0.
    pk_begin = int( pk_domain.get('begin') )
    pk_end = int( pk_domain.get('end') )
    pk_length = pk_end - pk_begin + 1
    if pk_length > 350 or pk_length < 191:
        pk_length_score -= 100.
    return pk_length_score

def overall_target_score_fn(pubs_score, pubs_score_max, cbioportal_mutations_score, cbioportal_mutations_score_max, disease_score, disease_score_max, bioassays_score, bioassays_score_max, pk_pdb_score, pk_pdb_score_max, pseudogene_score, pk_length_score):
    return (pubs_score / pubs_score_max) + (0.5 * (cbioportal_mutations_score / cbioportal_mutations_score_max)) + (0.5 * (disease_score / disease_score_max)) + (bioassays_score / bioassays_score_max) + (0.5 * (pk_pdb_score / pk_pdb_score_max)) + pseudogene_score + pk_length_score


# ==============
# Go through each kinase pk_domain, conduct initial scoring, and store as attribs
# ==============

for k in range(nkinases):
    kinase_node = kinDB[k]
    name = kinase_node.find('uniprot').get('entry_name')

    target_score_node = etree.SubElement(kinase_node, 'target_score')

    pk_domains = kinase_node.findall('uniprot/pk_domain')

    # Publications score
    pubs_score = publications_score_fn(kinase_node)
    target_score_node.set('publications', '%.3f' % pubs_score)

    # Mutants score
    cbioportal_mutations_score = cbioportal_mutations_score_fn(kinase_node)
    if cbioportal_mutations_score == None:
        target_score_node.set('cBioPortal_mutations', 'No data')
    else:
        target_score_node.set('cBioPortal_mutations', '%.3f' % cbioportal_mutations_score)

    # Disease association score
    disease_score = disease_score_fn(kinase_node)
    target_score_node.set('disease_association', '%.3f' % disease_score)

    # Bioassays score
    bioassays_score = bioassays_score_fn(kinase_node)
    target_score_node.set('bioassays', '%.3f' % bioassays_score)

    pk_pdbs = kinase_node.findall('pk_pdb')

    for pk_domain in pk_domains:
        pk_domain_id = pk_domain.get('id')
        kinDB_id = pk_domain.get('kinDB_id')

        # Set up a target_score/pk_domain node
        pk_domain_score_node = etree.SubElement(target_score_node, 'pk_domain')
        pk_domain_score_node.set('id', pk_domain_id)
        pk_domain_score_node.set('kinDB_id', kinDB_id)

        # Score for number of PDB structures including the pk_domain (only count the first chain)
        num_matching_pk_pdbs = 0
        for pk_pdb in pk_pdbs:
            matching_pk_pdb_chain = pk_pdb.find('chain[@pk_domain_id="%s"]' % pk_domain_id)
            if matching_pk_pdb_chain != None:
                num_matching_pk_pdbs += 1
        pk_pdb_score = pk_pdb_score_fn(num_matching_pk_pdbs)
        pk_domain_score_node.set('pk_pdbs', '%.3f' % pk_pdb_score)
        pk_domain_score_node.set('npk_pdbs', '%d' % num_matching_pk_pdbs)

        # Deprioritize pseudogenes
        pseudogene_score = pseudogene_score_fn(kinDB_id)
        pk_domain_score_node.set('pseudogene', '%.3f' % pseudogene_score)

        # Score for length of pk_domain (deprioritize unusually long or short pk_domains) - in the future this may be replaced with something more sophisticated which deprioritizes kinases with long inserts
        pk_length_score = pk_length_score_fn(pk_domain)
        pk_domain_score_node.set('pk_length', '%.3f' % pk_length_score)


# ==============
# Calculate score maxima (required for normalization)
# ==============

pubs_score_max = max( [float(x.get('publications')) for x in kinDB.findall('kinase/target_score')] )

cbioportal_mutations_score_max = max( [float(x.get('cBioPortal_mutations')) for x in kinDB.findall('kinase/target_score') if x.get('cBioPortal_mutations') != 'No data'] )

disease_score_max = max( [float(x.get('disease_association')) for x in kinDB.findall('kinase/target_score')] )

bioassays_score_max = max( [float(x.get('bioassays')) for x in kinDB.findall('kinase/target_score')] )

pk_pdb_score_max = max( [float(x.get('pk_pdbs')) for x in kinDB.findall('kinase/target_score/pk_domain')] )

# ==============
# Iterate through kinases again, calculating overall target_scores
# ==============

for k in range(nkinases):
    kinase_node = kinDB[k]
    name = kinase_node.find('uniprot').get('entry_name')
    pk_domains = kinase_node.findall('uniprot/pk_domain')
    kinase_score_node = kinase_node.find('target_score')
    for pk_domain in pk_domains:
        pk_domain_id = pk_domain.get('id')
        kinDB_id = pk_domain.get('kinDB_id')
        pk_domain_score_node = kinase_node.find('target_score/pk_domain[@kinDB_id="%(kinDB_id)s"]' % vars())

        # Get scores for individual parameters
        pubs_score = float(kinase_score_node.get('publications'))
        cbioportal_mutations_score = kinase_score_node.get('cBioPortal_mutations')
        if cbioportal_mutations_score == 'No data':
            cbioportal_mutations_score = 0.
        else:
            cbioportal_mutations_score = float(cbioportal_mutations_score)
        disease_score = float(kinase_score_node.get('disease_association'))
        bioassays_score = float(kinase_score_node.get('bioassays'))
        pk_pdb_score = float(pk_domain_score_node.get('pk_pdbs'))

        pseudogene_score = float(pk_domain_score_node.get('pseudogene'))
        pk_length_score = float(pk_domain_score_node.get('pk_length'))

        # Calculate overall target score and set value in XML DB
        overall_target_score = overall_target_score_fn(pubs_score, pubs_score_max, cbioportal_mutations_score, cbioportal_mutations_score_max, disease_score, disease_score_max, bioassays_score, bioassays_score_max, pk_pdb_score, pk_pdb_score_max, pseudogene_score, pk_length_score)
        pk_domain_score_node.set('target_score', str(overall_target_score))


# ==============
# Normalize the overall target scores
# ==============

max_score = max([ float(t.get('target_score')) for t in kinDB.findall('kinase/target_score/pk_domain') ])

for k in range(nkinases):
    kinase_node = kinDB[k]
    pk_domains = kinase_node.findall('uniprot/pk_domain')
    for pk_domain in pk_domains:
        kinDB_id = pk_domain.get('kinDB_id')
        pk_domain_score_node = kinase_node.find('target_score/pk_domain[@kinDB_id="%(kinDB_id)s"]' % vars())

        overall_target_score = float(pk_domain_score_node.get('target_score'))
        normalized_target_score_percentage = overall_target_score / max_score * 100.
        pk_domain_score_node.set('target_score', '%.1f' % normalized_target_score_percentage)


# ==============
# write the XML DB
# ==============

ofile = open(ofile, 'w')
ofile.write( etree.tostring(kinDB, pretty_print=True) )
ofile.close()


