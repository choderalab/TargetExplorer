# For each PDB in kinDB.xml, download the SIFTS residue-mapping .xml file
# Add experimental and resolved sequences to kinDB (numbered according to the UniProt sequence)
# Also add alignments of these sequences against the UniProt sequence
#
# Daniel L. Parton <partond@mskcc.org> - 10 Apr 2013
#
# Perhaps get DSSP info from here: http://www.rcsb.org/pdb/rest/das/pdbchainfeatures/features?segment=5pti.A
#
# TODO get secondary structure info from SIFTS files

#==============================================================================
# IMPORTS
#==============================================================================

import sys,os,gzip
from lxml import etree
from choderalab.pdb import retrieve_sifts
from choderalab.core import seqwrap
from Bio.PDB import to_one_letter_code

#==============================================================================
# PARAMETERS
#==============================================================================

kinDB_path = 'kinDB.xml'
okinDB_path = 'kinDB-pdb.xml'
structures_dir = os.path.join('..', 'structures')
local_sifts_path = os.path.join(structures_dir, 'sifts')

#==============================================================================
# MAIN
#==============================================================================

# Read in the kinDB XML document
print 'Reading', kinDB_path
parser = etree.XMLParser(remove_blank_text=True)
kinDB = etree.parse(kinDB_path, parser).getroot()
nkinases = len(kinDB)
print 'Number of kinases:', nkinases

for k in range(nkinases):
    # For each PDB in kinDB.xml, download the SIFTS residue-mapping .xml file if it is not already present
    pdb_nodes = kinDB[k].findall('pk_pdb')
    uniprot_sequence = kinDB[k].findtext('uniprot/sequence').strip()
    uniprot_sequence = ''.join(uniprot_sequence.split('\n'))
    uniprotAC = kinDB[k].find('uniprot').get('AC')
    entry_name = kinDB[k].find('uniprot').get('entry_name')
    #if uniprotAC != 'P31751':
    #if uniprotAC != 'Q00532':
    #if uniprotAC != 'P07333':
    #if uniprotAC != 'Q16539':
    #if uniprotAC != 'Q16539':
    #if uniprotAC != 'O14965':
    #    continue
    for pdb_node in pdb_nodes:
        pdbid = pdb_node.get('id')
        #if pdbid == '1O6L':
        #    sys.exit()
        local_sifts_file_path = os.path.join(local_sifts_path, pdbid+'.xml.gz')
        if os.path.exists(local_sifts_file_path):
            pass
        else:
            print 'Downloading SIFTS file and saving as (compressed):', local_sifts_file_path
            page = retrieve_sifts(pdbid)
            with gzip.open(local_sifts_file_path, 'wb') as local_sifts_file:
                local_sifts_file.write(page + '\n')

        # Parse the sifts XML document
        with gzip.open(local_sifts_file_path,'rb') as local_sifts_file:
            sifts = etree.parse(local_sifts_file, parser).getroot()

        # Get the chains to be searched from kinDB
        kinDB_chain_nodes = pdb_node.findall('chain')
        for chain_node in kinDB_chain_nodes:
            chainid = chain_node.get('id')
            # First check whether the first residue with matching chainid and a UniProt crossref has the same UniProt AC as was picked up from UniProt (by gather-uniprot.py).
            # 3O50 and 3O51 are picked up by gather-uniprot.py from uniprot AC O14965. But these have uniprot AC B4DX16 in the sifts .xml files, which is a TrEMBL entry. Sequences are almost identical except for deletion of ~70 residues prior to PK domain of B4DX16. This means that experimental_sequence_aln and related sequences are not added by gather-pdb.py. Need to sort out a special case for these pdbs. Should check for similar cases in other kinases.
            # 3O50 and 3O51 can be ignored. (Plenty of other PDBs for that protein)
            # 3OG7 is picked up from uniprot AC P15056, but the PDB entry links to Q5IBP5 - this is the AKAP9-BRAF fusion protein.
            # XXX TODO XXX 3OG7 will be ignored for now, but at some point should make separate entries for fusion proteins, and add the PDB files accordingly.
            first_matching_uniprot_resi = sifts.find('entity[@type="protein"]/segment/listResidue/residue/crossRefDb[@dbSource="PDB"][@dbChainId="%s"]/../crossRefDb[@dbSource="UniProt"]' % chainid)
            sifts_uniprotAC = first_matching_uniprot_resi.get('dbAccessionId')
            if uniprotAC != sifts_uniprotAC:
                print 'PDB %s chain %s picked up from UniProt entry %s %s. Non-matching UniProtAC in sifts: %s. This pk_pdb entry will be deleted when outputting %s' %  (pdbid, chainid, entry_name, uniprotAC, sifts_uniprotAC, okinDB_path)
                chain_node.set('DELETE_ME','')

            #
            #
            # TODO check if there are any PDBs where two proteins share the same chainid (I seem to remember that there are - check previous scripts)
            #
            #

            # Now extract the sequence data
            # These are the sifts residues which include a PDB crossref with matching chainid
            chain_residues = sifts.findall('entity[@type="protein"]/segment/listResidue/residue/crossRefDb[@dbSource="PDB"][@dbChainId="%s"]/..' % chainid)
            experimental_sequence = ''
            experimental_sequence_pdb_resids = []
            experimental_sequence_uniprot_res_indices = []
            observed_sequence = ''
            experimental_sequence_aln = ['-'] * len(uniprot_sequence) # This will contain the alignment of the experimental sequence against the full UniProt sequence. Conflicting residues will be added if they are contiguous with non-conflicting segments
            experimental_sequence_aln_conflicts = ['-'] * len(uniprot_sequence) # Same, but conflicting residues are added as 'x' if they have no UniProt crossref, or 'c' if they have a UniProt crossref but also have 'Conflict' residueDetail
            n_crossref_uniprot_matches = 0
            for r in chain_residues:
                residue_details = [ detail.text.strip() for detail in r.findall('residueDetail') ] # list of strings
                resname = r.attrib['dbResName'] 
                if resname == None:
                    print 'WARNING: UniProt crossref not found for conflicting residue!', k, pdbid, chainid, r.attrib
                    sys.exit()
                try:
                    # Note that this BioPython dict converts a modified aa to the single-letter code of its unmodified parent (e.g. "TPO":"T")
                    single_letter = to_one_letter_code[ resname ]
                except KeyError:
                    if r.attrib['dbResName'] == 'ACE': # Just ignore N-terminal ACE
                        continue
                # Add residue to experimental_sequence
                experimental_sequence += single_letter

                # Also save the pdb resids, which we will use later
                pdb_resid = r.find('crossRefDb[@dbSource="PDB"]').attrib['dbResNum']
                # Some pdb resids are e.g. '464A'
                if pdb_resid.isdigit() == False:
                    if pdbid in ['1O6L','2JDO','2JDR','2UW9','2X39','2XH5']: # These pdbs include three residues with pdb resids 464A, 464B, 464C, (all with UniProt crossrefs) then continue from 465. We will change this so that the pdb resids continue to iterate
                        corrected_pdb_resids = {'464A':465, '464B':466, '464C':467}
                        if pdb_resid in corrected_pdb_resids.keys():
                            pdb_resid = corrected_pdb_resids[pdb_resid]
                        elif int(pdb_resid[0:3]) > 464:
                            pdb_resid = int(pdb_resid) + 3
                    # Otherwise just extract the number (this will also detect negative numbers)
                    else:
                        pdb_resid = ''.join([char for char in pdb_resid if (char.isdigit() or char == '-')])
                try:
                    experimental_sequence_pdb_resids.append( int(pdb_resid) )
                except:
                    print uniprotAC, pdbid, chainid, pdb_resid
                    sys.exit()

                # Also add residue to experimental_sequence_aln. Residues which do not match the uniprot sequence (and thus do not have a uniprot crossref) will be added later
                crossref_uniprot = r.find('crossRefDb[@dbSource="UniProt"][@dbAccessionId="%s"]' % uniprotAC)
                if crossref_uniprot != None:
                    n_crossref_uniprot_matches += 1
                    index = int(crossref_uniprot.attrib['dbResNum']) - 1
                    experimental_sequence_aln[index] = single_letter
                    if 'Conflict' in residue_details:
                        experimental_sequence_aln_conflicts[index] = 'c'
                    else:
                        experimental_sequence_aln_conflicts[index] = single_letter
                    experimental_sequence_uniprot_res_indices.append(index)
                else:
                    experimental_sequence_uniprot_res_indices.append(None)
                    pass
                # Add residue to observed_sequence if it is observed, otherwise '-'
                if 'Not_Observed' in residue_details:
                    observed_sequence += '-'
                else:
                    observed_sequence += single_letter

            # ======
            # Now we add the residues which do not have a uniprot crossref
            # ======
            #print k, uniprotAC, pdbid, chainid
            #print experimental_sequence
            #print experimental_sequence_pdb_resids
            #print ''.join(experimental_sequence_aln_conflicts)
            i = 0

            # But first we have to deal with cases where residues have been added at the N-terminus which extend before the start of the uniprot sequence.
            # Get the uniprot residue index of the first residue with a uniprot crossref
            for s in range(len(experimental_sequence_uniprot_res_indices)):
                UP_res_index = experimental_sequence_uniprot_res_indices[s]
                if UP_res_index != None:
                    first_exp_seq_uniprot_res_index = UP_res_index
                    # And the corresponding pdb resid
                    corresponding_pdb_resid = experimental_sequence_pdb_resids[s]
                    exp_seq_first_uniprot_res_index = s
                    break
            # And get the pdb resid of the first residue in the experimental sequence
            for s in experimental_sequence_pdb_resids:
                if s != None:
                    first_exp_seq_pdb_resid = s
                    break
            ignore_excess_Nterm_residues_flag = False
            # If the experimental sequence includes the first residue of the full uniprot sequence
            if first_exp_seq_uniprot_res_index == 0:
                # And if the value of the first pdb resid is lower than that of the pdb resid corresponding to the first uniprot residue
                if first_exp_seq_pdb_resid < corresponding_pdb_resid:
                    # Then we will ignore the excess residues
                    ignore_excess_Nterm_residues_flag = True

            while i < len(experimental_sequence):
                resname_i = experimental_sequence[i]
                uniprot_res_index_i = experimental_sequence_uniprot_res_indices[i]
                pdb_resid_i = experimental_sequence_pdb_resids[i]

                if (ignore_excess_Nterm_residues_flag == True) and (pdb_resid_i < corresponding_pdb_resid):
                    pass # we ignore these residues

                # If this residue does not have a uniprot crossref
                elif uniprot_res_index_i == None:
                    # Start a list of residues with no uniprot crossref
                    contiguous_noUP_residues = [ resname_i ]
                    # Then check the next residue
                    j = i + 1
                    while j < len(experimental_sequence):
                        resname_j = experimental_sequence[j]
                        uniprot_res_index_j = experimental_sequence_uniprot_res_indices[j]
                        pdb_resid_j = experimental_sequence_pdb_resids[j]
                        #print 'len, i, j:', len(experimental_sequence), i, j, pdb_resid_i, pdb_resid_j, contiguous_noUP_residues

                        # If this residue also has no uniprot crossref, and is contiguous in terms of pdb resnum, then add it to the list, and move on to the next one
                        if (uniprot_res_index_j == None) and ((pdb_resid_j - pdb_resid_i) == (j-i)):
                            #print 'adding to list:', j, resname_j
                            contiguous_noUP_residues.append( resname_j )
                            pass

                        # If this residue does have a uniprot crossref, and if it is contiguous in terms of pdb resnum, then we add the list of residues without uniprot crossrefs at this position
                        elif (uniprot_res_index_j != None) and ((pdb_resid_j - pdb_resid_i) == (j-i)):
                            #print 'adding to sequence_aln:', j
                            experimental_sequence_aln[ (uniprot_res_index_j - j) : uniprot_res_index_j ] = contiguous_noUP_residues
                            experimental_sequence_aln_conflicts[ (uniprot_res_index_j - j) : uniprot_res_index_j ] = 'x' * len(contiguous_noUP_residues)
                            i = j
                            break

                        # If this residue is not contiguous in terms of pdb resnum, go back and check if the first of contiguous_noUP_residues is pdb-contiguous with the previous residue - if so, add contiguous_noUP_residues
                        elif (pdb_resid_j - pdb_resid_i) != (j-i):
                            #print 'checking backwards:', j
                            if (pdb_resid_i - experimental_sequence_pdb_resids[i-1]) == 1:
                                last_uniprot_res_index = experimental_sequence_uniprot_res_indices[i-1]
                                experimental_sequence_aln[ last_uniprot_res_index + 1 : last_uniprot_res_index + 1 + (j-i)] = contiguous_noUP_residues
                                experimental_sequence_aln_conflicts[ last_uniprot_res_index + 1 : last_uniprot_res_index + 1 + (j-i)] = 'x' * len(contiguous_noUP_residues)
                            i = j - 1
                            break

                        # If we have reached the end of experimental_sequence, go back and check if the first of contiguous_noUP_residues is pdb-contiguous with the previous residue - if so, add contiguous_noUP_residues
                        if j == len(experimental_sequence) - 1:
                            #print 'THIS IS THE END', len(experimental_sequence), i, j, pdb_resid_i, experimental_sequence_pdb_resids[i], experimental_sequence_pdb_resids[i-1], contiguous_noUP_residues
                            #print experimental_sequence_pdb_resids
                            if (pdb_resid_i - experimental_sequence_pdb_resids[i-1]) == 1:
                                last_uniprot_res_index = experimental_sequence_uniprot_res_indices[i-1]
                                experimental_sequence_aln[ last_uniprot_res_index + 1 : last_uniprot_res_index + 2 + (j-i)] = contiguous_noUP_residues
                                experimental_sequence_aln_conflicts[ last_uniprot_res_index + 1 : last_uniprot_res_index + 2 + (j-i)] = 'x' * len(contiguous_noUP_residues)
                            i = j
                            break
                        j += 1

                i += 1

            # In cases such as 3LAU and 1O6L, additional sequence at end makes experimental_sequence_aln longer than uniprot_sequence by 1
            if len(experimental_sequence_aln) != len(uniprot_sequence):
                experimental_sequence_aln = experimental_sequence_aln[0:len(uniprot_sequence)]
                experimental_sequence_aln_conflicts = experimental_sequence_aln_conflicts[0:len(uniprot_sequence)]
                
            #print ''.join(experimental_sequence_aln_conflicts)

            # Now add the various sequence data to kinDB
            experimental_sequence_aln = ''.join(experimental_sequence_aln)
            experimental_sequence_aln_conflicts = ''.join(experimental_sequence_aln_conflicts)
            #print k, pdbid, chainid, len(experimental_sequence), len(observed_sequence)
            #print k, pdbid, chainid, experimental_sequence_aln
            exp = etree.SubElement(chain_node,'experimental')
            etree.SubElement(exp,'sequence').text = '\n' + seqwrap(experimental_sequence)
            exp.set('length', str(len(experimental_sequence)))
            etree.SubElement(exp,'sequence_aln').text = '\n' + seqwrap(experimental_sequence_aln)
            etree.SubElement(exp,'sequence_aln_conflicts').text = '\n' + seqwrap(experimental_sequence_aln_conflicts)
            obs = etree.SubElement(chain_node,'observed')
            etree.SubElement(obs,'sequence').text = '\n' + seqwrap(observed_sequence)

            #if pdbid == '2W1C':
            #if pdbid == '3LAU':
            #if pdbid == '1O6L':
            if pdbid == '3O50':
                #sys.exit()
                pass

        # Only add if the chain matches that in the kinDB

# =======================
# Delete pk_pdb/chain entries with @DELETE_ME attrib. These were cases where the sifts_uniprotAC did not match the uniprotAC in kinDB (derived from the UniProt entry by gather-uniprot.py)
# =======================
pk_pdbs_with_chains_to_be_deleted = set( kinDB.findall('kinase/pk_pdb/chain[@DELETE_ME=""]/..') )
for p in pk_pdbs_with_chains_to_be_deleted:
    chains_to_be_deleted = p.findall('chain[@DELETE_ME=""]')
    for d in chains_to_be_deleted:
        p.remove(d)
    # If the pk_pdb node now has no children (almost certainly the case), delete that too
    if len(p.getchildren()) == 0:
        kinase = p.getparent()
        kinase.remove(p)

# write the XML DB
#ofile = open(okinDB_path , 'w')
#ofile.write( etree.tostring(kinDB, pretty_print=True) )
#ofile.close()

print 'Done.'

