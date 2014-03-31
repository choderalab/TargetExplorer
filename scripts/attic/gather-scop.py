# Retrieve kinase domain structure info from SCOP, and merge with existing XML
# database. PDB IDs are matched to UniProtACs via sifts.
# UniProtAC is obtained by searching www.rcsb.org/pdb
#
# Daniel L. Parton <partond@mskcc.org> - 31 Mar 2013

#==============================================================================
# IMPORTS
#==============================================================================

import sys,os.path,urllib2,gzip,ast
from lxml import etree
from Bio import SCOP
from choderalab.pdb import retrieve_pdb
from choderalab.uniprot import get_uniprot_mapping

#==============================================================================
# PARAMETERS
#==============================================================================

scop_sccs_search_term = 'd.144.1.7' # Protein kinases, catalytic subunit
kinDB_path = os.path.join('database', 'kinDB-pdb.xml') # Path for existing kinDB XML database
structures_dir = os.path.join('..','structures')
sifts_dir = os.path.join(structures_dir,'sifts')
sifts_map_pdb_uniprot_path = os.path.join(sifts_dir,'pdb_chain_uniprot.csv')
okinDB_path = os.path.join('database', 'kinDB-pdb-scop.xml')

#==============================================================================
# SCOP file sources.
#==============================================================================

scop_dir = os.path.join(structures_dir,'scop')
scop_version = '1.75B'
scop_cla_filename = os.path.join(scop_dir,'dir.cla.scop.%(scop_version)s.txt' % vars())
scop_des_filename = os.path.join(scop_dir,'dir.des.scop.%(scop_version)s.txt' % vars())

#==============================================================================
# Build a list of SCOP domain info for entries which match the sccs search
# term, and a list of associated UniProtACs
#==============================================================================

print 'Retrieving data from local SCOP database files...\n'

# First open sifts file for mapping PDB IDs to UniProt ACs
# sifts_map_pdb_uniprot will be of form:
# {'101m' : {'A' : ['P02185','1','154','0','153','1','154']}}
with open(sifts_map_pdb_uniprot_path, 'r') as sifts_map_pdb_uniprot_file:
    sifts_map_pdb_uniprot = dict()
    for line in sifts_map_pdb_uniprot_file.readlines():
        if line[0:10] != 'PDB,CHAIN,':
            words = line.strip().split(',')
            try:
                sifts_map_pdb_uniprot[words[0]]
            except KeyError:
                sifts_map_pdb_uniprot[words[0]] = dict()
            sifts_map_pdb_uniprot[words[0]][words[1]] = words[2:]

scop_info = list()
scop_cla_file = open(scop_cla_filename)
for scop_entry in SCOP.Cla.parse(scop_cla_file):
    if scop_entry.sccs == scop_sccs_search_term:
        # If the sccs term matches, first grab some data from the line
        pdb_id = scop_entry.residues.pdbid
        species = scop_entry.hierarchy['sp'] # this returns the sunid for this species node - note that species nodes are organized as descendents of subfamilies, so "Human" will have many different 'sp' nodes
        domain_def = scop_entry.residues.fragments  # tuple of form (('A','241','531'),)
        chain_id = domain_def[0][0]
        # In many cases there is no residue start/stop info present, but in those cases res_start and res_stop will be set to an empty string
        scop_res_start = domain_def[0][1]
        scop_res_stop = domain_def[0][2]
        sunid = scop_entry.sunid  # unique SCOP identifier
        sid = scop_entry.sid

        # Use sifts mapping file to get UniProtAC, as well as res starts and stops which use UniProt numbering
        uniprotAC = sifts_map_pdb_uniprot[pdb_id][chain_id][0]
        res_start = sifts_map_pdb_uniprot[pdb_id][chain_id][5]
        res_stop = sifts_map_pdb_uniprot[pdb_id][chain_id][6]

        # Uncomment to debug
        #print {'sunid':sunid,'pdb_id':pdb_id,'chain_id':chain_id,'res_start':res_start,'res_stop':res_stop}

        scop_info.append({'sunid':sunid,'sid':sid,'pdb_id':pdb_id,'chain_id':chain_id,'res_start':res_start,'res_stop':res_stop,'uniprotAC':uniprotAC,'species':species})

nscop = len(scop_info)
print 'Number of SCOP domain entries which have been retrieved:', nscop

scop_cla_file.close()

#==============================================================================
# Retrieve descriptions from 'sp' (species) nodes in SCOP Des file.
# This is used later to check we are not missing any human kinase domains.
#==============================================================================

scop_sp_descriptions = dict()
scop_des_file = open(scop_des_filename)
for scop_des in SCOP.Des.parse(scop_des_file):
    if scop_des.nodetype == 'sp':
        scop_sp_descriptions[scop_des.sunid] = scop_des.description

#==============================================================================
# Retrieve a list of UniProtACs from the currently existing XML database
#==============================================================================

print 'Retrieving UniProtACs from the existing XML database...'
parser = etree.XMLParser(remove_blank_text=True)
kinDB = etree.parse(kinDB_path, parser).getroot()
kinDB_kinases = kinDB.xpath('kinase')
nkinDB = len(kinDB_kinases)

#==============================================================================
# Iterate through each kinase in the XML database. If the kinase uniprotAC
# matches that associated with one of the SCOP domain entries, then the
# pertinent info from that SCOP domain entry will be stored in the XML
# database.
#==============================================================================

print 'Adding SCOP data to XML database...'

for s in range(nscop):
    # Search for the uniprotAC within kinDB
    matching_kinase = kinDB.xpath('./kinase/uniprot[@AC="%s"]/..' % scop_info[s]['uniprotAC'])
    # Should either find the kinase or not, so matching_kinase should be a list of len 0 or 1.
    if len(matching_kinase) == 0:
        # All of these are non-human, except for:
        # d1k9aa3 and related - mislabelled as Human in the SCOP Des file (it is actually rattus norvegicus)
        # d3o50a_ and related - the PDB entry links to an unreviewed TrEMBL entry (B4DX16; deposited 2008). Not much information, so can probably ignore.
        # d3fpqa_ and related - mislabelled as Human in the SCOP Des file (it is actually rattus norvegicus)
        # d3og7a_ and related - this is the kinase domain of the AKAP9-BRAF fusion protein. Only recently discovered. Can cause thyroid cancers. Still working out what to do about fusion proteins.
        print 'No match found in %s for scop sid: %s; species description: %s' % ( kinDB_path, scop_info[s]['sid'], scop_sp_descriptions[ scop_info[s]['species'] ] )
        continue
    elif len(matching_kinase) > 1:
        print "This really shouldn't happen", scop_info[s]
    else:
        matching_kinase = matching_kinase[0]
    scopnode = etree.SubElement(matching_kinase, 'scop')
    scopdomain = etree.SubElement(scopnode, 'scop_domain')
    scopdomain.set('sunid', str(scop_info[s]['sunid']))
    scopdomain.set('sid', scop_info[s]['sid'])
    scopdomain.set('pdb_id', scop_info[s]['pdb_id'])
    scopdomain.set('chain_id', scop_info[s]['chain_id'])
    scopdomain.set('res_start', str(scop_info[s]['res_start']))
    scopdomain.set('res_stop', str(scop_info[s]['res_stop']))

#==============================================================================
# Checking
#==============================================================================

print 'Number of SCOP domains added to database:', len(kinDB.xpath('./kinase/scop'))
nscop_disease = 0
for k in range(nkinDB):
    kinase = kinDB_kinases[k]
    #print kinase.find('scop'), kinase.find('uniprot/disease_association'), ((kinase.find('scop') != None) and (kinase.find('uniprot/disease_association') != None))
    if ((kinase.find('scop') != None) and (kinase.find('uniprot/disease_association') != None)):
        nscop_disease += 1
print 'Number of kinases with both SCOP domains and disease_annotations:', nscop_disease

#==============================================================================
# Output new kinDB
#==============================================================================

okinDB = open(okinDB_path, 'w')
okinDB.write( etree.tostring(kinDB, pretty_print=True) )
okinDB.close()

