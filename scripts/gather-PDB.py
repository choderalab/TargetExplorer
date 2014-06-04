#!/usr/bin/env python
#
# For each PDB in database-stage.xml, download the SIFTS residue-mapping .xml file
# Add experimental and resolved sequences to the DB (numbered according to the UniProt sequence)
# Also add alignments of these sequences against the UniProt sequence
#
# Daniel L. Parton <partond@mskcc.org> - 10 Apr 2013
#
# Perhaps get DSSP info from here: http://www.rcsb.org/pdb/rest/das/pdbchainfeatures/features?segment=5pti.A
#

#==============================================================================
# Imports
#==============================================================================

import sys,os,gzip,re,datetime,copy,urllib2,traceback
from lxml import etree
import TargetExplorer as clab
import Bio.PDB
import Bio.Data.SCOPData
from multiprocessing import Pool

#==============================================================================
# Parameters
#==============================================================================

if '-stage' in sys.argv:
    run_mode = 'stage'
elif '-dev' in sys.argv:
    run_mode = 'dev'
else:
    run_mode = 'nowrite'

print 'Running in mode: %s' % run_mode

database_dir = 'database'
external_data_dir = 'external-data'
local_pdb_dir = os.path.join(external_data_dir, 'PDB')
local_sifts_dir = os.path.join(external_data_dir, 'SIFTS')
if not os.path.exists(local_pdb_dir):
    os.mkdir(local_pdb_dir)
if not os.path.exists(local_sifts_dir):
    os.mkdir(local_sifts_dir)

DBstage_filepath = os.path.join(database_dir, 'database-stage.xml')
if not os.path.exists(DBstage_filepath):
    raise Exception, '%s not found.' % DBstage_filepath

if run_mode != 'nowrite':
    DB_out_filename = 'database-%(run_mode)s.xml' % vars()
    DB_out_filepath = os.path.join(database_dir, DB_out_filename)

verbose = False

now = datetime.datetime.utcnow()
datestamp = now.strftime(clab.DB.datestamp_format_string)

parser = etree.XMLParser(remove_blank_text=True)

#==============================================================================
# Definitions
#==============================================================================

class Gather_PDB_Results(object):
    '''
    Container class returned by the method gather_pdb
    '''
    def __init__(self, UniProt_entry_name=None):
        self.entry_name = UniProt_entry_name
        self.structures = []
    def add_structure_results(self, structures):
        self.structures.append(structures)

class Structure_Data(object):
    '''
    Class for communicating PDB structure data
    '''
    def __init__(self, structureID=None, exception_message=None):
        self.structureID = structureID
        self.exception_message = exception_message
        self.chains = []
        self.expression_data = None
        self.pdbcompoundID = None
    def add_chain_results(self, chain_results):
        self.chains.append(chain_results)
    def add_expression_data(self, expression_data):
        self.expression_data = expression_data
    def add_pdbcompoundID(self, pdbcompoundID):
        self.pdbcompoundID = pdbcompoundID

class Chain_Data(object):
    '''
    Class for communicating PDB chain data
    '''
    def __init__(self, chainID=None, experimental_sequence=None, experimental_sequence_aln=None, experimental_sequence_aln_conflicts=None, observed_sequence_aln_exp=None, observed_sequence_aln=None, ss_aln=None, exception_message=None):
        self.chainID = chainID
        self.experimental_sequence=experimental_sequence
        self.experimental_sequence_aln=experimental_sequence_aln
        self.experimental_sequence_aln_conflicts=experimental_sequence_aln_conflicts
        self.observed_sequence_aln_exp=observed_sequence_aln_exp
        self.observed_sequence_aln=observed_sequence_aln
        self.ss_aln=ss_aln
        self.exception_message=exception_message

def gather_pdb(e):

    # ======
    # Get some initial data for this DB entry
    # ======

    pdb_node = DB_root[e].find('PDB')
    if pdb_node == None:
        return None

    structure_nodes = pdb_node.findall('structure')
    uniprot_node = DB_root[e].find('UniProt')
    uniprot_sequence = uniprot_node.findtext('isoforms/canonical_isoform/sequence').strip()
    uniprot_sequence = ''.join(uniprot_sequence.split('\n'))
    uniprotAC = uniprot_node.get('AC')
    entry_name = uniprot_node.get('entry_name')

    #if entry_name != 'MLKL_HUMAN':
    #    return None

    # results will be added to this container object and returned at the end of the function
    #entry_results = []
    results_obj = Gather_PDB_Results(UniProt_entry_name=entry_name)

    # ======
    # Iterate through the PDB structures listed in the DB entry
    # ======

    for structure_node in structure_nodes:
        pdbid = structure_node.get('ID')

        #if pdbid != '2ITN':
        #    continue

        # ========
        # Get PDB and SIFTS files
        # ========

        # TODO define this via project metadata .yaml file.
        structure_paths = ['/Users/partond/kinome-MSMSeeder/structures/pdb', '/Users/partond/kinome-MSMSeeder/structures/sifts']

        local_pdb_filepath = os.path.join('external-data', 'PDB', pdbid + '.pdb')
        local_sifts_filepath = os.path.join('external-data', 'SIFTS', pdbid + '.xml.gz')

        # Check if PDB file/symlink already exists and is not empty
        search_for_pdb = True
        if os.path.exists(local_pdb_filepath):
            if os.path.getsize(local_pdb_filepath) > 0:
                search_for_pdb = False

        # If not, search any user-defined paths and create a symlink if found
        if search_for_pdb:
            for structure_dir in structure_paths:
                pdb_filepath = os.path.join(structure_dir, pdbid + '.pdb')
                if os.path.exists(pdb_filepath):
                    if os.path.getsize(pdb_filepath) > 0:
                        if os.path.exists(local_pdb_filepath):
                            os.remove(local_pdb_filepath)
                        os.symlink(pdb_filepath, local_pdb_filepath)
                        break

            # If still not found, download the PDB file
            if not os.path.exists(local_pdb_filepath):
                print 'Downloading PDB file and saving as:', local_pdb_filepath
                page = clab.PDB.retrieve_pdb(pdbid, compressed='yes')
                # download and write compressed file, read in and decompress with gzip package, write decompressed file, delete compressed file.
                with open(local_pdb_filepath + '.gz', 'wb') as local_pdb_file:
                    local_pdb_file.write(page)
                with gzip.open(local_pdb_filepath + '.gz', 'rb') as local_pdb_gz_file:
                    local_pdb_gz_contents = local_pdb_gz_file.read()
                    with open(local_pdb_filepath, 'w') as local_pdb_file:
                        local_pdb_file.write(local_pdb_gz_contents)
                os.remove(local_pdb_filepath + '.gz')

        # Check if SIFTS file already exists and is not empty
        search_for_sifts = True
        if os.path.exists(local_sifts_filepath):
            if os.path.getsize(local_sifts_filepath) > 0:
                search_for_sifts = False

        # If not, search any user-defined paths and create a symlink if found
        if search_for_sifts:
            for structure_dir in structure_paths:
                sifts_filepath = os.path.join(structure_dir, pdbid + '.xml.gz')
                if os.path.exists(sifts_filepath):
                    if os.path.getsize(sifts_filepath) > 0:
                        if os.path.exists(local_sifts_filepath):
                            os.remove(local_sifts_filepath)
                        os.symlink(sifts_filepath, local_sifts_filepath)
                        break

            # If still not found, download the SIFTS XML file
            if not os.path.exists(local_sifts_filepath):
                print 'Downloading SIFTS file (compressed) and saving as:', local_sifts_filepath
                try:
                    page = clab.PDB.retrieve_sifts(pdbid)
                except urllib2.URLError as urlerror:
                    if urlerror.reason == 'ftp error: [Errno ftp error] 550 Failed to change directory.':
                        # Check the PDB file has definitely been downloaded. If so, then the problem is probably that the SIFTS people have not yet created the file for this PDB entry, or they have not added it to their server yet.
                        if os.path.exists(local_pdb_filepath):
                            # In this case, just add a message telling the script to delete this PDB structure from the DB. The continue clause skips to the end of the function.
                            print '%s SIFTS file could not be downloaded - this PDB entry will be deleted from the DB' % pdbid
                            structure_results_obj = Structure_Data(structureID=pdbid)
                            structure_results_obj.exception_message='DELETE_ME - SIFTS file could not be downloaded'
                            results_obj.add_structure_results(structure_results_obj)
                            continue
                        else:
                            raise urlerror
                    else:
                        raise urlerror

                with gzip.open(local_sifts_filepath, 'wb') as local_sifts_file:
                    local_sifts_file.write(page)

        sifts = etree.fromstring( gzip.open(local_sifts_filepath, 'r').read() )

        # structure results will be added to this object
        structure_results_obj = Structure_Data(structureID=pdbid)

        # ======
        # From the PDB file, get the EXPRESSION_SYSTEM and related fields, using Bio.PDB.PDBParser
        # ======

        # First get the chains from DB_root
        DB_chain_nodes = structure_node.findall('chain')
        # and the chain IDs
        DB_chainIDs_lower = [chain_node.get('ID').lower() for chain_node in DB_chain_nodes]

        pdbparser = Bio.PDB.PDBParser(QUIET=True)
        #pdbparser = Bio.PDB.PDBParser(QUIET=False)
        pdbdata = pdbparser.get_structure(pdbid, local_pdb_filepath)
        pdbheader = pdbparser.get_header()
        # Bio PDB compound structure: {'compound': {'1': {'chain': 'a, b'}}}
        pdbcompounds = pdbheader['compound']
        #skip_structure = False
        matching_pdbcompoundID = None
        for pdbcompoundID in pdbcompounds.keys():
            try:
                for pdbchainID in pdbcompounds[pdbcompoundID]['chain'].split(', '):
                    if pdbchainID in DB_chainIDs_lower:
                        matching_pdbcompoundID = pdbcompoundID
                        break
            except Exception as e:
                print 'ERROR for entry %s PDB %s. PDB header dict as parsed by BioPython follows:' % (entry_name, pdbid)
                print pdbheader
                print traceback.format_exc()

                # XXX uncommenting this and the skip_structure commands means that the problematic PDB file will be deleted, and a new one downloaded. The script continues at this point without re-parsing the new PDB file, so the results will be wrong. But then can re-run the script with the proper PDB files.
                #os.remove(local_pdb_filepath)
                #page = clab.PDB.retrieve_pdb(pdbid, compressed='yes')
                ## download and write compressed file, read in and decompress with gzip package, write decompressed file, delete compressed file.
                #with open(local_pdb_filepath + '.gz', 'wb') as local_pdb_file:
                #    local_pdb_file.write(page)
                #with gzip.open(local_pdb_filepath + '.gz', 'rb') as local_pdb_gz_file:
                #    local_pdb_gz_contents = local_pdb_gz_file.read()
                #    with open(local_pdb_filepath, 'w') as local_pdb_file:
                #        local_pdb_file.write(local_pdb_gz_contents)
                #os.remove(local_pdb_filepath + '.gz')
                #skip_structure = True

                raise e

        #if skip_structure:
        #    continue

        expression_data = {}
        if matching_pdbcompoundID != None:
            # Bio PDB source structure: {'source': {'1': {'expression_system': 'escherichia coli'}}}
            pdbexpression_data = pdbheader['source'][matching_pdbcompoundID]
            for key in pdbexpression_data.keys():
                if key[0:10] == 'expression':
                    # Make expression data upper-case again. I think it looks better for single-case text.
                    expression_data[key.upper()] = pdbexpression_data[key].upper()

        # ======
        # Iterate through chains
        # ======

        for c, chain_node in enumerate(DB_chain_nodes):
            #DELETE_ME = False
            exception_message = None
            chainID = chain_node.get('ID')
            if verbose: print entry_name, uniprotAC, pdbid, chainID
            # First check whether the first residue with matching chainID and a UniProt crossref has the same UniProt AC as was picked up from UniProt (by gather-uniprot.py).
            # 3O50 and 3O51 are picked up by gather-uniprot.py from uniprot AC O14965. But these have uniprot AC B4DX16 in the sifts .xml files, which is a TrEMBL entry. Sequences are almost identical except for deletion of ~70 residues prior to PK domain of B4DX16. This means that experimental_sequence_aln and related sequences are not added by gather-pdb.py. Need to sort out a special case for these pdbs. Should check for similar cases in other kinases.
            # 3O50 and 3O51 can be ignored. (Plenty of other PDBs for that protein)
            # 3OG7 is picked up from uniprot AC P15056, but the PDB entry links to Q5IBP5 - this is the AKAP9-BRAF fusion protein.
            # XXX TODO XXX 3OG7 will be ignored for now, but at some point should make separate entries for fusion proteins, and add the PDB files accordingly.
            if verbose: print sifts
            first_matching_uniprot_resi = sifts.find('entity[@type="protein"]/segment/listResidue/residue/crossRefDb[@dbSource="PDB"][@dbChainId="%s"]/../crossRefDb[@dbSource="UniProt"]' % chainID)
            sifts_uniprotAC = first_matching_uniprot_resi.get('dbAccessionId')
            if uniprotAC != sifts_uniprotAC:
                print 'PDB %s chain %s picked up from UniProt entry %s %s. Non-matching UniProtAC in sifts: %s. This chain will be deleted before writing the database.' %  (pdbid, chainID, entry_name, uniprotAC, sifts_uniprotAC)
                #DELETE_ME = True
                exception_message = 'DELETE_ME'

            #
            #
            # TODO check if there are any PDBs where two proteins share the same chainID (I seem to remember that there are - check previous scripts)
            #
            #

            # ======
            # Extract sequence data from the SIFTS XML
            # ======

            # These are the sifts residues which include a PDB crossref with matching chainID
            chain_residues = sifts.findall('entity[@type="protein"]/segment/listResidue/residue/crossRefDb[@dbSource="PDB"][@dbChainId="%s"]/..' % chainID)
            experimental_sequence = ''
            experimental_sequence_pdb_resids = []
            experimental_sequence_uniprot_res_indices = []
            observed_sequence_aln_exp = ''
            experimental_sequence_aln = ['-'] * len(uniprot_sequence) # This will contain the alignment of the experimental sequence against the full UniProt sequence. Conflicting residues will be added if they are contiguous with non-conflicting segments. NOTE: this is no longer added to the database.
            experimental_sequence_aln_conflicts = ['-'] * len(uniprot_sequence) # Same, but conflicting residues are added as lower case
            observed_sequence_aln = ['-'] * len(uniprot_sequence) # This will contain the alignment of the observed sequence against the full UniProt sequence. Conflicting residues will be ignored.
            ss_aln = ['-'] * len(uniprot_sequence) # This will contain the alignment of the secondary structure codes against the full UniProt sequence. Conflicting residues will be ignored.
            n_crossref_uniprot_matches = 0
            for r in chain_residues:
                residue_details = r.findall('residueDetail')
                residue_detail_texts = [ detail.text.strip() for detail in residue_details ] # list of strings
                ss = r.findtext('residueDetail[@property="codeSecondaryStructure"]')
                resname = r.attrib['dbResName'] 
                if resname == None:
                    print 'ERROR: UniProt crossref not found for conflicting residue!', e, pdbid, chainID, r.attrib
                    raise Exception
                try:
                    # Note that this BioPython dict converts a modified aa to the single-letter code of its unmodified parent (e.g. "TPO":"T")
                    single_letter = Bio.Data.SCOPData.protein_letters_3to1[ resname ]
                except KeyError:
                    if resname == 'ACE': # Just ignore N-terminal ACE
                        continue
                    elif resname == 'CAS': # S-(dimethylarsenic)cysteine
                        single_letter = 'C'
                    elif resname == 'MHO': # S-oxymethionine
                        single_letter = 'M'
                    elif resname == 'LGY': # 3NX8. (E)-N-(4-oxobutylidene)lysine
                        single_letter = 'K'
                    elif resname == 'AME': # N-acetylmethionine
                        single_letter = 'M'
                    elif resname == 'NMM': # 3KB7
                        single_letter = 'R'
                    elif resname == 'OCY': # 2R9S
                        single_letter = 'C'
                    elif resname == 'CY0': # 2J5E
                        single_letter = 'C'
                    elif resname == 'CY7': # 2JIV
                        single_letter = 'C'
                    else:
                        print 'KeyError: Problem converting resname', resname, 'to single letter code.', e, pdbid, chainID, r.attrib
                        raise KeyError
                # Add residue to experimental_sequence
                experimental_sequence += single_letter

                # Also save the pdb resids, which we will use later
                pdb_resid = r.find('crossRefDb[@dbSource="PDB"]').attrib['dbResNum']
                # Some pdb resids are e.g. '464A'
                if pdb_resid.isdigit() == False:
                    if pdbid in ['1O6L','2JDO','2JDR','2UW9','2X39','2XH5']: # These pdbs include three residues with pdb resids 464A, 464B, 464C, (all with UniProt crossrefs) then continues from 465. We will change this so that the pdb resids continue to iterate
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
                    print 'Problem converting pdb_resid into int.' , uniprotAC, pdbid, chainID, pdb_resid
                    raise Exception

                # Also add residue to experimental_sequence_aln. Residues which do not match the uniprot sequence (and thus do not have a uniprot crossref) will be added later
                crossref_uniprot = r.find('crossRefDb[@dbSource="UniProt"][@dbAccessionId="%s"]' % uniprotAC)
                if crossref_uniprot != None:
                    n_crossref_uniprot_matches += 1
                    index = int(crossref_uniprot.attrib['dbResNum']) - 1
                    experimental_sequence_aln[index] = single_letter
                    if 'Conflict' in residue_detail_texts or 'Engineered mutation' in residue_detail_texts:
                        experimental_sequence_aln_conflicts[index] = single_letter.lower()
                    else:
                        experimental_sequence_aln_conflicts[index] = single_letter
                    experimental_sequence_uniprot_res_indices.append(index)
                    # Add residue to observed_sequence_aln if it is observed and is not a conflict
                    if 'Not_Observed' not in residue_detail_texts and ('Conflict' not in residue_detail_texts or 'Engineered mutation' in residue_detail_texts):
                        observed_sequence_aln[index] = single_letter
                        if ss != None:
                            ss_aln[index] = ss
                else:
                    experimental_sequence_uniprot_res_indices.append(None)
                    pass
                # Add residue to observed_sequence_aln_exp if it is observed, otherwise '-'
                if 'Not_Observed' in residue_detail_texts:
                    observed_sequence_aln_exp += '-'
                else:
                    observed_sequence_aln_exp += single_letter

            # Now check whether the number of non-observed residues is more than 90% of the experimental sequence length
            n_unobserved_residues = observed_sequence_aln_exp.count('-')
            if ( float(n_unobserved_residues) / float(len(experimental_sequence)) ) > 0.9:
                #DELETE_ME = True
                exception_message = 'DELETE_ME'

            # ======
            # Now we add the residues which do not have a UniProt crossref
            # ======

            #print e, uniprotAC, pdbid, chainID
            #print experimental_sequence
            #print ''.join(experimental_sequence_aln_conflicts)

            i = 0

            # But first we have to deal with cases where residues have been added at the N-terminus which extend before the start of the uniprot sequence. The excess residues will be ignored.
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
            try:
                if first_exp_seq_uniprot_res_index == 0:
                    # And if the value of the first pdb resid is lower than that of the pdb resid corresponding to the first uniprot residue
                    if first_exp_seq_pdb_resid < corresponding_pdb_resid:
                        # Then we will ignore the excess residues
                        ignore_excess_Nterm_residues_flag = True
            except:
                # XXX should do something better than this
                # exception occurs with P27791 (KAPCA_RAT)
                exception_message = 'DELETE_ME'

            # Now iterate through the residues in the experimental sequence and add residues which do not have a uniprot crossref, but are contiguous in terms of PDB numbering

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
                            experimental_sequence_aln_conflicts[ (uniprot_res_index_j - j) : uniprot_res_index_j ] = list(''.join(contiguous_noUP_residues).lower())
                            i = j
                            break

                        # If this residue is not contiguous in terms of pdb resnum, go back and check if the first of contiguous_noUP_residues is pdb-contiguous with the previous residue - if so, add contiguous_noUP_residues
                        elif (pdb_resid_j - pdb_resid_i) != (j-i):
                            #print 'checking backwards:', j
                            if (pdb_resid_i - experimental_sequence_pdb_resids[i-1]) == 1:
                                last_uniprot_res_index = experimental_sequence_uniprot_res_indices[i-1]
                                experimental_sequence_aln[ last_uniprot_res_index + 1 : last_uniprot_res_index + 1 + (j-i)] = contiguous_noUP_residues
                                experimental_sequence_aln_conflicts[ last_uniprot_res_index + 1 : last_uniprot_res_index + 1 + (j-i)] = list(''.join(contiguous_noUP_residues).lower())
                            i = j - 1
                            break

                        # If we have reached the end of experimental_sequence, go back and check if the first of contiguous_noUP_residues is pdb-contiguous with the previous residue - if so, add contiguous_noUP_residues
                        if j == len(experimental_sequence) - 1:
                            #print 'THIS IS THE END', len(experimental_sequence), i, j, pdb_resid_i, experimental_sequence_pdb_resids[i], experimental_sequence_pdb_resids[i-1], contiguous_noUP_residues
                            #print experimental_sequence_pdb_resids
                            if (pdb_resid_i - experimental_sequence_pdb_resids[i-1]) == 1:
                                last_uniprot_res_index = experimental_sequence_uniprot_res_indices[i-1]
                                experimental_sequence_aln[ last_uniprot_res_index + 1 : last_uniprot_res_index + 2 + (j-i)] = contiguous_noUP_residues
                                experimental_sequence_aln_conflicts[ last_uniprot_res_index + 1 : last_uniprot_res_index + 2 + (j-i)] = list(''.join(contiguous_noUP_residues).lower())
                            i = j
                            break
                        j += 1

                i += 1

            # ======
            # Some final processing
            # ======

            # In cases such as 3LAU and 1O6L, additional sequence at end makes experimental_sequence_aln longer than uniprot_sequence by 1
            # Handle this by removing the extraneous sequence
            if len(experimental_sequence_aln) != len(uniprot_sequence):
                experimental_sequence_aln = experimental_sequence_aln[0:len(uniprot_sequence)]
                experimental_sequence_aln_conflicts = experimental_sequence_aln_conflicts[0:len(uniprot_sequence)]
                
            experimental_sequence_aln = ''.join(experimental_sequence_aln)
            experimental_sequence_aln_conflicts = ''.join(experimental_sequence_aln_conflicts)
            observed_sequence_aln = ''.join(observed_sequence_aln)
            ss_aln = ''.join(ss_aln)

            # ======
            # Now finish constructing the results container object before returning
            # ======

            # XXX
            #chain_result = [experimental_sequence, experimental_sequence_aln, experimental_sequence_aln_conflicts, observed_sequence_aln_exp, observed_sequence_aln, ss_aln, DELETE_ME]
            chain_results_obj = Chain_Data(chainID=chainID, experimental_sequence=experimental_sequence, experimental_sequence_aln=experimental_sequence_aln, experimental_sequence_aln_conflicts=experimental_sequence_aln_conflicts, observed_sequence_aln_exp=observed_sequence_aln_exp, observed_sequence_aln=observed_sequence_aln, ss_aln=ss_aln, exception_message=exception_message)
            #structure_results[c] = chain_result
            structure_results_obj.add_chain_results(chain_results_obj)

        structure_results_obj.add_expression_data(expression_data)
        structure_results_obj.add_pdbcompoundID(matching_pdbcompoundID)
        results_obj.add_structure_results(structure_results_obj)

    #return entry_results
    return results_obj






# ====================
# Main
# ====================

if __name__ == '__main__':

    # Read in the exsiting DB
    print 'Reading', DBstage_filepath
    DB_root = etree.parse(DBstage_filepath, parser).getroot()
    nentries = len(DB_root)
    print 'Number of entries:', nentries

    entry_range = range(nentries)
    # Use multiprocessor pool to retrieve various data from the PDB, for each entry in the DB
    pool = Pool()
    results = pool.map(gather_pdb, entry_range)
    # results = map(gather_pdb, entry_range)   # serial version, for debugging

    for result in results:

        if result == None:
            continue

        entry_name = result.entry_name
        DB_PDB_node = DB_root.find('entry/UniProt[@entry_name="%s"]/../PDB' % entry_name)
        if DB_PDB_node == None:
            continue

        structures = result.structures

        for structure_result in structures:

            PDB_ID = structure_result.structureID
            DB_structure_node = DB_PDB_node.find('structure[@ID="%s"]' % PDB_ID)

            if structure_result.exception_message:
                if structure_result.exception_message == 'DELETE_ME - SIFTS file could not be downloaded':
                    DB_structure_node.set('DELETE_ME','')
                continue

            # Remove any existing data derived from gather-pdb.py before adding new data
            DB_expression_data_node = DB_structure_node.find('expression_data')
            if DB_expression_data_node != None:
                DB_structure_node.remove(DB_expression_data_node)

            chains = structure_result.chains

            for chain_result in chains:

                chainID = chain_result.chainID
                DB_chain_node = DB_structure_node.find('chain[@ID="%s"]' % chainID)

                if chain_result.exception_message:
                    if chain_result.exception_message == 'DELETE_ME':
                        DB_chain_node.set('DELETE_ME','')

                # Remove any existing data derived from gather-pdb.py before adding new data
                exp = DB_chain_node.find('experimental_sequence')
                if exp != None:
                    DB_chain_node.remove(exp)
                obs = DB_chain_node.find('observed_sequence')
                if obs != None:
                    DB_chain_node.remove(obs)

                # Add chain data to DB
                exp = etree.SubElement(DB_chain_node, 'experimental_sequence')
                etree.SubElement(exp, 'sequence').text = '\n' + clab.core.seqwrap(chain_result.experimental_sequence)
                exp.set('length', str(len(chain_result.experimental_sequence)))
                #etree.SubElement(exp, 'sequence_aln').text = '\n' + clab.core.seqwrap(chain_result.experimental_sequence_aln) # NOTE: this is no longer added to the database
                etree.SubElement(exp, 'sequence_aln_conflicts').text = '\n' + clab.core.seqwrap(chain_result.experimental_sequence_aln_conflicts)
                obs = etree.SubElement(DB_chain_node, 'observed_sequence')
                etree.SubElement(obs, 'sequence_aln_exp').text = '\n' + clab.core.seqwrap(chain_result.observed_sequence_aln_exp)
                etree.SubElement(obs, 'sequence_aln').text = '\n' + clab.core.seqwrap(chain_result.observed_sequence_aln)
                etree.SubElement(obs, 'ss_aln').text = '\n' + clab.core.seqwrap(chain_result.ss_aln)

            # Add expression data and the pdbcompoundID to DB
            expression_data = structure_result.expression_data
            if verbose: print expression_data
            DB_expression_data_node = etree.Element('expression_data')
            for expression_data_key in expression_data.keys():
                DB_expression_data_node.set(expression_data_key, expression_data[expression_data_key])
            DB_structure_node.insert(0, DB_expression_data_node)
            DB_structure_node.set('pdbcompoundID', str(structure_result.pdbcompoundID))

    # =======================
    # Delete PDB structure and chain entries with @DELETE_ME attrib. These were cases where the sifts_uniprotAC did not match the uniprotAC in DB_root (derived from the UniProt entry by gather-uniprot.py), or where more than 90% of the experimental sequence was unobserved
    # =======================

    structures_to_be_deleted = set( DB_root.findall('entry/PDB/structure[@DELETE_ME=""]') )
    for s in structures_to_be_deleted:
        PDB_node = s.getparent()
        PDB_node.remove(s)

    structures_with_chains_to_be_deleted = set( DB_root.findall('entry/PDB/structure/chain[@DELETE_ME=""]/..') )
    for s in structures_with_chains_to_be_deleted:
        chains_to_be_deleted = s.findall('chain[@DELETE_ME=""]')
        for c in chains_to_be_deleted:
            s.remove(c)
        # If the PDB node now has no children, delete that too
        if len(s.getchildren()) == 0:
            kinase = s.getparent()
            kinase.remove(s)

    # =======================
    # If staging, update date_run
    # =======================

    if run_mode == 'stage':
        DB_root.set('gather_pdb_last_run', datestamp)

    #==============================================================================
    # If staging, compare new and old DBs
    #==============================================================================

    data_modified = False

    if run_mode == 'stage':
        # Parse the old DB
        DBold_root = etree.parse(DB_out_filepath, parser).getroot()

        # First a quick check to see whether the numbers of PDB experimental sequence nodes match
        DB_seq_nodes = DB_root.findall('entry/PDB/structure/chain/experimental_sequence/sequence')
        DBold_seq_nodes = DBold_root.findall('entry/PDB/structure/chain/experimental_sequence/sequence')
        if len(DB_seq_nodes) != len(DBold_seq_nodes):
            print 'Comparison of latest PDB data with data in %s indicates changes. DB will be re-written with new data. Number of experimental sequence nodes differs by : %d' % (DB_out_filepath, len(DB_seq_nodes) - len(DBold_seq_nodes))
            data_modified = True

        # If the numbers of PDB chain nodes match, then proceeed to a full comparison of PDB-derived data using diff
        else:
            DB_comparison_root = etree.Element('database')
            DBold_comparison_root = etree.Element('database')

            for entry in DB_root:
                pdb_node = entry.find('PDB')
                if pdb_node != None:
                    DB_comparison_entry_node = etree.SubElement(DB_comparison_root, 'entry')
                    DB_comparison_entry_node.append(copy.deepcopy(pdb_node))

            DB_comparison_string = etree.tostring(DB_comparison_root, pretty_print=True)
            DB_comparison_string = '\n'.join( DB_comparison_string.splitlines()[1:] ) # Ignore the first line, which contains datestamps (don't want to include these in the comparison)

            for entry in DBold_root:
                pdb_node = entry.find('PDB')
                if pdb_node != None:
                    DBold_comparison_entry_node = etree.SubElement(DBold_comparison_root, 'entry')
                    DBold_comparison_entry_node.append(copy.deepcopy(pdb_node))

            DBold_comparison_string = etree.tostring(DBold_comparison_root, pretty_print=True)
            DBold_comparison_string = '\n'.join( DBold_comparison_string.splitlines()[1:] ) # Ignore the first line, which contains datestamps (don't want to include these in the comparison)

            # Now compare the two comparison strings using GNU diff
            diff_output = clab.DB.diff_DB_comparison_strings(DBold_comparison_string, DB_comparison_string)

            if len(diff_output) > 0:
                print 'Comparison of latest PDB data with data in %s indicates changes. File will be rewritten with new data. Lines in diff comparison: %s' % (DB_out_filepath, len(diff_output))
                data_modified = True
                if verbose:
                    print DB_comparison_string.split('\n')[0:20]
                    print DBold_comparison_string.split('\n')[0:20]

            else:
                print 'Comparison of latest PDB data with data in %s indicates no changes. File will be rewritten with updated gather_pdb_last_run attrib, but other data will not be modified.' % DB_out_filepath
                data_modified = False


    #==============================================================================
    # If staging and there have been modifications, update date_modified
    #==============================================================================

    if run_mode == 'stage' and data_modified:
        DB_root.set('gather_pdb_last_modif', datestamp)

    # =======================
    # write the XML DB
    # =======================
    if run_mode == 'stage' and data_modified:
        clab.DB.writeDB(DB_root, DB_out_filepath)

    elif run_mode == 'stage' and not data_modified:
        DBold_root = etree.parse(DB_out_filepath, parser).getroot()
        DBold_root.set('gather_pdb_last_run', datestamp)
        clab.DB.writeDB(DBold_root, DB_out_filepath)

    elif run_mode == 'dev':
        clab.DB.writeDB(DB_root, DB_out_filepath)

    print ''
    print 'Done.'

