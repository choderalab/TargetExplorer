#!/usr/bin/env python
#
# Assign scores to each target domain in the DB, based on various properties such as
# number of PDB structures, references, disease associations etc.
#
# Daniel L. Parton <partond@mskcc.org> - 29 April 2013
#

import sys, os, datetime
from lxml import etree
from numpy import *
import TargetExplorer as clab

# ==============
# Parameters
# ==============

if '-stage' in sys.argv:
    run_mode = 'stage'
elif '-dev' in sys.argv:
    run_mode = 'dev'
else:
    run_mode = 'nowrite'

print 'Running in mode: %s' % run_mode

database_dir = 'database'

DBstage_filepath = os.path.join(database_dir, 'database-stage.xml')
if not os.path.exists(DBstage_filepath):
    raise Exception, '%s not found.' % DBstage_filepath

if run_mode != 'nowrite':
    DB_out_filename = 'database-%(run_mode)s.xml' % vars()
    DB_out_filepath = os.path.join(database_dir, DB_out_filename)

now = datetime.datetime.utcnow()
datestamp = now.strftime(clab.DB.datestamp_format_string)

parser = etree.XMLParser(remove_blank_text=True)

target_taxID = 9606

# ==============
# Scoring functions
# ==============

def publications_score_fn(DBentry_node):
    '''
    Scoring function based on the number of PubMed publications retrieved from NCBI Gene.
    '''
    pubs_node = DBentry_node.find('NCBI_Gene/entry/publications')
    if pubs_node == None:
        return 0.
    else:
        return float(len(pubs_node))

def cbioportal_mutations_score_fn(DBentry_node):
    '''
    Scoring function based on the number of mutations for a given kinase/mutants node
    '''
    percent_cases_with_mutations = clab.cBioPortal.percent_cases_with_mutations(DBentry_node)
    return percent_cases_with_mutations

def disease_score_fn(DBentry_node):
    '''
    Scoring function based on the number of disease associations for a certain kinase.
    '''
    ndisease_assocs = len( DBentry_node.findall('UniProt/disease_associations/disease_association') )
    score = log(ndisease_assocs + 1)
    return score

def bioassays_score_fn(DBentry_node):
    '''
    Scoring function based on the number of disease associations for a certain kinase.
    '''
    nbioassays = len( DBentry_node.findall('bioassays/bioassay') )
    score = log(nbioassays + 1)
    return score

def PDB_score_fn(x):
    '''
    Scoring function based on the number of PDB entries for a certain PK domain.
    '''
    y = log(x+1)
    return y

def pseudogene_score_fn(targetID):
    '''
    Deprioritizes certain genes thought to be pseudogenes.

    PRKY_HUMAN_D0: see comments on UniProt and Ensembl page
    '''    
    if targetID in ['PRKY_HUMAN_D0', 'PDPK2_HUMAN_D0']:
        return -100.
    else:
        return 0.

def target_domain_length_score_fn(target_domain):
    '''
    Deprioritizes certain genes thought to be pseudogenes.

    PRKY_HUMAN_O43930: see comments on UniProt and Ensembl page
    '''    
    target_domain_length_score = 0.
    target_domain_begin = int( target_domain.get('begin') )
    target_domain_end = int( target_domain.get('end') )
    target_domain_length = target_domain_end - target_domain_begin + 1
    if target_domain_length > 350 or target_domain_length < 191:
        target_domain_length_score -= 100.
    return target_domain_length_score

def overall_target_score_fn(pubs_score, pubs_score_max, cbioportal_mutations_score, cbioportal_mutations_score_max, disease_score, disease_score_max, bioassays_score, bioassays_score_max, PDB_score, PDB_score_max, pseudogene_score, target_domain_length_score):
    # print pubs_score_max, cbioportal_mutations_score_max, disease_score_max, bioassays_score_max, PDB_score_max
    target_score = (pubs_score / pubs_score_max) + (0.5 * (cbioportal_mutations_score / cbioportal_mutations_score_max)) + (0.5 * (disease_score / disease_score_max)) + (bioassays_score / bioassays_score_max) + (0.5 * (PDB_score / PDB_score_max)) + pseudogene_score + target_domain_length_score
    return target_score


# =================
# Read in existing database
# =================

DB_root = etree.parse(DBstage_filepath, parser).getroot()
nentries = len(DB_root)
target_entries = DB_root.findall('entry/UniProt[@NCBI_taxID="%d"]/..' % target_taxID)
ntarget_entries = len(target_entries)

# ==============
# Go through each target domain, conduct initial scoring, and store as attribs
# ==============

for e in range(ntarget_entries):
    DBentry_node = target_entries[e]
    name = DBentry_node.find('UniProt').get('entry_name')

    target_score_node = etree.SubElement(DBentry_node, 'target_score')

    target_domains = DBentry_node.findall('UniProt/domains/domain[@targetID]')

    # Publications score
    pubs_score = publications_score_fn(DBentry_node)
    target_score_node.set('publications', '%.3f' % pubs_score)

    # Mutants score
    cbioportal_mutations_score = cbioportal_mutations_score_fn(DBentry_node)
    if cbioportal_mutations_score == None:
        target_score_node.set('cBioPortal_mutations', 'No data')
    else:
        target_score_node.set('cBioPortal_mutations', '%.3f' % cbioportal_mutations_score)

    # Disease association score
    disease_score = disease_score_fn(DBentry_node)
    target_score_node.set('disease_association', '%.3f' % disease_score)

    # Bioassays score
    bioassays_score = bioassays_score_fn(DBentry_node)
    target_score_node.set('bioassays', '%.3f' % bioassays_score)

    PDBs = DBentry_node.findall('PDB/structure')

    for target_domain in target_domains:
        target_domainID = target_domain.get('domainID') # TODO this will change to domainID at some point
        targetID = target_domain.get('targetID')

        # Set up a target_score/domain node
        target_domain_score_node = etree.SubElement(target_score_node, 'domain')
        target_domain_score_node.set('domainID', target_domainID)
        target_domain_score_node.set('targetID', targetID)

        # Score for number of PDB structures including the pk_domain (only count the first chain)
        num_matching_PDBs = 0
        for PDB_structure in PDBs:
            matching_PDB_chain = PDB_structure.find('chain[@domainID="%s"]' % target_domainID)
            if matching_PDB_chain != None:
                num_matching_PDBs += 1
        PDB_score = PDB_score_fn(num_matching_PDBs)
        target_domain_score_node.set('PDBs', '%.3f' % PDB_score)
        target_domain_score_node.set('nPDBs', '%d' % num_matching_PDBs)

        # Deprioritize pseudogenes
        pseudogene_score = pseudogene_score_fn(targetID)
        target_domain_score_node.set('pseudogene', '%.3f' % pseudogene_score)

        # Score for length of target_domain (deprioritize unusually long or short target_domains) - in the future this may be replaced with something more sophisticated which deprioritizes kinases with long inserts
        target_domain_length_score = target_domain_length_score_fn(target_domain)
        target_domain_score_node.set('target_domain_length', '%.3f' % target_domain_length_score)


# ==============
# Calculate score maxima (required for normalization)
# ==============

pubs_score_max = max( [float(x.get('publications')) for x in DB_root.findall('entry/target_score')] )

cbioportal_mutations_score_max = max( [float(x.get('cBioPortal_mutations')) for x in DB_root.findall('entry/target_score') if x.get('cBioPortal_mutations') != 'No data'] )

disease_score_max = max( [float(x.get('disease_association')) for x in DB_root.findall('entry/target_score')] )

bioassays_score_max = max( [float(x.get('bioassays')) for x in DB_root.findall('entry/target_score')] )

PDB_score_max = max( [float(x.get('PDBs')) for x in DB_root.findall('entry/target_score/domain')] )

# ==============
# Iterate through kinases again, calculating overall target_scores
# ==============

for e in range(ntarget_entries):
    DBentry_node = target_entries[e]
    name = DBentry_node.find('UniProt').get('entry_name')

    target_score_node = DBentry_node.find('target_score')

    target_domains = DBentry_node.findall('UniProt/domains/domain[@targetID]')

    for target_domain in target_domains:
        target_domainID = target_domain.get('domainID') # TODO this will change to domainID at some point
        targetID = target_domain.get('targetID')
        target_domain_score_node = target_score_node.find('domain[@targetID="%s"]' % targetID)

        # Get scores for individual parameters
        pubs_score = float(target_score_node.get('publications'))
        cbioportal_mutations_score = target_score_node.get('cBioPortal_mutations')
        if cbioportal_mutations_score == 'No data':
            cbioportal_mutations_score = 0.
        else:
            cbioportal_mutations_score = float(cbioportal_mutations_score)
        disease_score = float(target_score_node.get('disease_association'))
        bioassays_score = float(target_score_node.get('bioassays'))
        pk_pdb_score = float(target_domain_score_node.get('PDBs'))

        pseudogene_score = float(target_domain_score_node.get('pseudogene'))
        target_domain_length_score = float(target_domain_score_node.get('target_domain_length'))

        # Calculate overall target score and set value in XML DB
        overall_target_score = overall_target_score_fn(pubs_score, pubs_score_max, cbioportal_mutations_score, cbioportal_mutations_score_max, disease_score, disease_score_max, bioassays_score, bioassays_score_max, PDB_score, PDB_score_max, pseudogene_score, target_domain_length_score)
        target_domain_score_node.set('target_score', str(overall_target_score))

# ==============
# Normalize the overall target scores and add ranks
# ==============

max_score = max([ float(t.get('target_score')) for t in DB_root.findall('entry/target_score/domain') ])

# ranks will be generated by sorting this list of targetIDs
all_targetIDs = [ t.get('targetID') for t in DB_root.findall('entry/target_score/domain') ]
dict_for_sorting = { t.get('targetID') : t.get('target_score') for t in DB_root.findall('entry/target_score/domain') }
all_targetIDs_sorted = sorted(all_targetIDs, key = lambda x: dict_for_sorting[x], reverse=True)

for e in range(ntarget_entries):
    DBentry_node = target_entries[e]
    target_domains = DBentry_node.findall('UniProt/domains/domain[@targetID]')
    for target_domain in target_domains:
        targetID = target_domain.get('targetID')
        target_domain_score_node = DBentry_node.find('target_score/domain[@targetID="%s"]' % targetID)

        overall_target_score = float(target_domain_score_node.get('target_score'))
        normalized_target_score_percentage = overall_target_score / max_score * 100.
        target_domain_score_node.set('target_score', '%.1f' % normalized_target_score_percentage)

        # ranks
        target_domain_rank = all_targetIDs_sorted.index(targetID) + 1
        target_domain_score_node.set('target_rank', '%d' % target_domain_rank)




# =======================
# If staging, update date_run
# =======================

if run_mode == 'stage':
    DB_root.set('prioritization_last_run', datestamp)

#==============================================================================
# If staging, compare new and old DBs
#==============================================================================

data_modified = False

if run_mode == 'stage':
    # Parse the old DB
    DBold_root = etree.parse(DB_out_filepath, parser).getroot()

    # Just check if the overall scores match
    DB_target_scores = [float(score_node.get('target_score')) for score_node in DB_root.findall('entry/target_score/domain')]
    DBold_target_scores = [float(score_node.get('target_score')) for score_node in DBold_root.findall('entry/target_score/domain')]
    if DB_target_scores != DBold_target_scores:
        print 'Comparison of latest prioritization data with data in %s indicates changes. DB will be re-written with new data.' % DB_out_filepath
        data_modified = True

    # If the numbers of PDB chain nodes match, then proceeed to a full comparison of PDB-derived data using diff
    else:
        print 'Comparison of latest prioritization data with data in %s indicates no changes. File will be rewritten with updated prioritization_last_run attrib, but other data will not be modified.' % DB_out_filepath
        data_modified = False


#==============================================================================
# If staging and there have been modifications, update date_modified
#==============================================================================

if run_mode == 'stage' and data_modified:
    DB_root.set('prioritization_last_modif', datestamp)

# =======================
# write the XML DB
# =======================
if run_mode == 'stage' and data_modified:
    clab.DB.writeDB(DB_root, DB_out_filepath)

elif run_mode == 'stage' and not data_modified:
    DBold_root = etree.parse(DB_out_filepath, parser).getroot()
    DBold_root.set('prioritization_last_run', datestamp)
    clab.DB.writeDB(DBold_root, DB_out_filepath)

elif run_mode == 'dev':
    clab.DB.writeDB(DB_root, DB_out_filepath)

print ''
print 'Done.'






