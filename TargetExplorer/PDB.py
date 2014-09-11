# Methods for retrieving data from the PDB
#
# Daniel L. Parton <partond@mskcc.org> - 16 Feb 2013

#==============================================================================
# IMPORTS
#==============================================================================

import urllib2

#==============================================================================
# METHODS
#==============================================================================

def retrieve_sifts(pdb_id):
    '''Retrieves a SIFTS .xml file, given a PDB ID. Works by modifying the PDBe download URL.
    File is retrieved compressed and returned uncompressed.
    Also removes annoying namespace stuff.
    '''
    import re, gzip, StringIO
    sifts_download_base_url='ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/xml/'
    url = sifts_download_base_url + pdb_id.lower() + '.xml.gz'

    # NOTE: can do exception handling for this as shown in the commented lines. In general this is done in the parent class.
    #try:
    response = urllib2.urlopen(url)
    #except urllib2.URLError as urlerror:
    #    if urlerror.reason = 'ftp error: [Errno ftp error] 550 Failed to change directory.':
    #        raise 
    #    else:
    #        raise urlerror

    sifts_page = response.read(100000000) # Max 100MB
    # Decompress string
    sifts_page = gzip.GzipFile(fileobj=StringIO.StringIO(sifts_page)).read()

    # Remove all attribs from the entry tag, and the rdf tag and contents
    sifts_page_processed = ''
    skip_rdf_tag_flag = False
    for line in sifts_page.splitlines():
        if line[0:6] == '<entry':
            sifts_page_processed += '<entry>' + '\n'
        elif line[0:7] == '  <rdf:':
            skip_rdf_tag_flag = True
            pass
        elif line[0:8] == '  </rdf:':
            skip_rdf_tag_flag = False
            pass
        else:
            if skip_rdf_tag_flag:
                continue
            sifts_page_processed += line + '\n'
    return sifts_page_processed

def retrieve_pdb(pdb_id,compressed='no'):
    '''Retrieves a PDB file, given a PDB ID. Works by modifying the PDB download URL.
    '''
    pdb_download_base_url='http://www.rcsb.org/pdb/files/'
    url = pdb_download_base_url + pdb_id + '.pdb'
    if compressed == 'yes':
        url += '.gz'
    response = urllib2.urlopen(url)
    pdb_file = response.read(10000000) # Max 10MB
    return pdb_file

def extract_sifts_seq(sifts_filepath, uniprot_ac, uniprot_entry_name, pdbid, chain_id, uniprot_sequence):
    from lxml import etree
    import gzip
    import Bio.Data.SCOPData

    exception_message = None

    sifts = etree.fromstring( gzip.open(sifts_filepath, 'r').read() )

    # First check whether the first residue with matching chainID and a UniProt crossref has the same UniProt AC as was picked up from UniProt (by gather-uniprot.py).
    # 3O50 and 3O51 are picked up by gather-uniprot.py from uniprot AC O14965. But these have uniprot AC B4DX16 in the sifts .xml files, which is a TrEMBL entry. Sequences are almost identical except for deletion of ~70 residues prior to PK domain of B4DX16. This means that experimental_sequence_aln and related sequences are not added by gather-pdb.py. Need to sort out a special case for these pdbs. Should check for similar cases in other kinases.
    # 3O50 and 3O51 can be ignored. (Plenty of other PDBs for that protein)
    # 3OG7 is picked up from uniprot AC P15056, but the PDB entry links to Q5IBP5 - this is the AKAP9-BRAF fusion protein.
    # XXX TODO XXX 3OG7 will be ignored for now, but at some point should make separate entries for fusion proteins, and add the PDB files accordingly.

    first_matching_uniprot_resi = sifts.find('entity[@type="protein"]/segment/listResidue/residue/crossRefDb[@dbSource="PDB"][@dbChainId="%s"]/../crossRefDb[@dbSource="UniProt"]' % chain_id)
    sifts_uniprot_ac = first_matching_uniprot_resi.get('dbAccessionId')
    if uniprot_ac != sifts_uniprot_ac:
        print 'PDB %s chain %s picked up from UniProt entry %s %s. Non-matching UniProtAC in sifts: %s. This chain will be deleted when rewriting the database.' %  (pdbid, chain_id, uniprot_entry_name, uniprot_ac, sifts_uniprot_ac)
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
    chain_residues = sifts.findall('entity[@type="protein"]/segment/listResidue/residue/crossRefDb[@dbSource="PDB"][@dbChainId="%s"]/..' % chain_id)
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
            print 'ERROR: UniProt crossref not found for conflicting residue!', uniprot_ac, pdbid, chain_id, r.attrib
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
                print 'KeyError: Problem converting resname', resname, 'to single letter code.', chain_id, r.attrib
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
            print 'Problem converting pdb_resid into int.', uniprot_ac, pdbid, chain_id, pdb_resid
            raise Exception

        # Also add residue to experimental_sequence_aln. Residues which do not match the uniprot sequence (and thus do not have a uniprot crossref) will be added later
        crossref_uniprot = r.find('crossRefDb[@dbSource="UniProt"][@dbAccessionId="%s"]' % uniprot_ac)
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
        exception_message = 'DELETE_ME'

    # ======
    # Now we add the residues which do not have a UniProt crossref
    # ======

    #print e, uniprot_ac, pdbid, chain_id
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

        chain_results_dict = {
            'chain_id': chain_id,
            'experimental_seq': experimental_sequence,
            'experimental_seq_aln_conflicts': experimental_sequence_aln_conflicts,
            'observed_seq_aln_exp': observed_sequence_aln_exp,
            'observed_seq_aln': observed_sequence_aln,
            'observed_ss_aln': ss_aln,
            'exception_message': exception_message,
        }
        return chain_results_dict
