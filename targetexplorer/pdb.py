import os
import datetime
import gzip
import StringIO
import urllib2
import traceback
from lxml import etree
import Bio.PDB
import Bio.Data.SCOPData
# from multiprocessing import Pool
from targetexplorer.flaskapp import models, db
from targetexplorer.core import logger


class GatherPDB(object):
    def __init__(self, structure_dirs=None, run_main=True, commit_to_db=True):
        self.commit_to_db = commit_to_db
        if type(structure_dirs) == str:
            self.structure_dirs = [structure_dirs]
        else:
            self.structure_dirs = structure_dirs

        if run_main:
            self.main()

    def main(self):
        external_data_dir = os.path.join('external-data')
        local_pdb_dir = os.path.join(external_data_dir, 'PDB')
        local_sifts_dir = os.path.join(external_data_dir, 'SIFTS')
        if not os.path.exists(local_pdb_dir):
            os.mkdir(local_pdb_dir)
        if not os.path.exists(local_sifts_dir):
            os.mkdir(local_sifts_dir)

        now = datetime.datetime.utcnow()

        # get current crawl number
        crawldata_row = models.CrawlData.query.first()
        current_crawl_number = crawldata_row.current_crawl_number
        print 'Current crawl number: %d' % current_crawl_number

        verbose = False

        # ====================
        # Main
        # ====================

        # first get a list of PDB rows from the db
        db_pdb_rows = models.PDBEntry.query.filter_by(crawl_number=current_crawl_number).all()

        # prepare a dict for each PDB row to be passed to the function extract_pdb_data
        db_pdb_dicts = []
        for pdb_row in db_pdb_rows:
            db_entry = models.DBEntry.query.filter_by(id=pdb_row.db_entry_id).first()
            uniprot_row = db_entry.uniprot.first()
            canon_isoform_row = db_entry.uniprot_isoforms.filter_by(is_canonical=True).first()
            chain_data = [
                {
                    'chain_row_id': chain_row.id,
                    'chain_id': chain_row.chain_id
                } for chain_row in pdb_row.chains
            ]
            pdb_dict = {
                'pdb_row_id': pdb_row.id,
                'pdb_id': pdb_row.pdb_id,
                'ac': uniprot_row.ac,
                'entry_name': uniprot_row.entry_name,
                'seq': canon_isoform_row.sequence,
                'chain_data': chain_data,
                'structure_dirs': self.structure_dirs
            }
            db_pdb_dicts.append(pdb_dict)

        # TODO multiprocessing is causing problems with tests. Maybe refactor to use mpi4py instead?
        # Use multiprocessor pool to retrieve various data for each PDB
        # pool = Pool()
        results = map(extract_pdb_data, db_pdb_dicts)

        for pdb_results in results:
            pdb_row_id = pdb_results['pdb_row_id']
            pdb_row = models.PDBEntry.query.filter_by(id=pdb_row_id).first()

            if 'exception_message' in pdb_results and pdb_results.get('exception_message') == 'SIFTS file could not be downloaded':
                db.session.delete(pdb_row)
                continue

            for chain_row_id in pdb_results['chain_dicts']:
                chain_row = models.PDBChain.query.filter_by(id=chain_row_id).first()
                chain_dict = pdb_results['chain_dicts'][chain_row_id]
                for key in chain_dict:
                    if key in ['chain_id', 'exception_message']:
                        continue
                    setattr(chain_row, key, chain_dict[key])

                # Delete chain entries with exception message of 'DELETE_ME'.
                # These are cases where the sifts_uniprotAC does not match the
                # uniprotAC in DB_root (derived from the UniProt entry by
                # gather-uniprot.py), or where more than 90% of the
                # experimental sequence is unobserved
                if chain_dict['exception_message'] == 'DELETE_ME':
                    db.session.delete(chain_row)

            # If all chain entries have been deleted, delete the entire PDB row.
            if pdb_row.chains.count() == 0:
                db.session.delete(pdb_row)

        current_crawl_datestamp_row = models.DateStamps.query.filter_by(crawl_number=current_crawl_number).first()
        current_crawl_datestamp_row.pdb_datestamp = now
        if self.commit_to_db:
            db.session.commit()
        print 'Done.'










class Extract_PDB_Results(object):
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


def extract_pdb_data(pdb_dict):
    '''Extract data for a single PDB structure
    '''
    pdb_row_id = pdb_dict['pdb_row_id']
    pdb_id = pdb_dict['pdb_id']
    ac = pdb_dict['ac']
    entry_name = pdb_dict['entry_name']
    seq = pdb_dict['seq']
    chain_data = pdb_dict['chain_data']
    structure_dirs = pdb_dict['structure_dirs']

    # if entry_name != 'MLKL_HUMAN':
    #     return None
    #
    # if pdb_id != '2ITN':
    #     return None

    # ========
    # Get PDB and SIFTS files
    # PDB files are used to extract expression system metadata
    # SIFTS files are used to extract sequence data
    # ========

    # TODO define this via project metadata .yaml file.
    # structure_dirs = ['/Users/partond/tmp/kinome-MSMSeeder/structures/pdb', '/Users/partond/tmp/kinome-MSMSeeder/structures/sifts']

    local_pdb_filepath = os.path.join('external-data', 'PDB', pdb_id + '.pdb.gz')
    local_sifts_filepath = os.path.join('external-data', 'SIFTS', pdb_id + '.xml.gz')

    # Check if PDB file/symlink already exists and is not empty
    search_for_pdb = True
    if os.path.exists(local_pdb_filepath):
        if os.path.getsize(local_pdb_filepath) > 0:
            search_for_pdb = False

    # If not, search any user-defined paths and create a symlink if found
    if search_for_pdb:
        if structure_dirs:
            for structure_dir in structure_dirs:
                pdb_filepath = os.path.join(structure_dir, pdb_id + '.pdb.gz')
                if os.path.exists(pdb_filepath):
                    if os.path.getsize(pdb_filepath) > 0:
                        if os.path.exists(local_pdb_filepath):
                            os.remove(local_pdb_filepath)
                        os.symlink(pdb_filepath, local_pdb_filepath)
                        break

        # If still not found, download the PDB file
        if not os.path.exists(local_pdb_filepath):
            print 'Downloading PDB file and saving as:', local_pdb_filepath
            page = retrieve_pdb(pdb_id, compressed='yes')
            # download and write compressed file
            with open(local_pdb_filepath, 'wb') as local_pdb_file:
                local_pdb_file.write(page)

    # Check if SIFTS file already exists and is not empty
    search_for_sifts = True
    if os.path.exists(local_sifts_filepath):
        if os.path.getsize(local_sifts_filepath) > 0:
            search_for_sifts = False

    # If not, search any user-defined paths and create a symlink if found
    if search_for_sifts:
        if structure_dirs:
            for structure_dir in structure_dirs:
                sifts_filepath = os.path.join(structure_dir, pdb_id + '.xml.gz')
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
                page = retrieve_sifts(pdb_id)
            except urllib2.URLError as urlerror:
                if urlerror.reason == 'ftp error: [Errno ftp error] 550 Failed to change directory.':
                    # Check the PDB file has definitely been downloaded. If so, then the problem is probably that the SIFTS people have not yet created the file for this PDB entry, or they have not added it to their server yet.
                    if os.path.exists(local_pdb_filepath):
                        # In this case, just add a message telling the script to delete this PDB structure from the DB. The continue clause skips to the end of the function.
                        print '%s SIFTS file could not be downloaded - this PDB entry will be deleted from the DB' % pdb_id
                        return {'pdb_row_id': pdb_row_id, 'exception_message': 'SIFTS file could not be downloaded'}
                    else:
                        raise urlerror
                else:
                    raise urlerror

            with gzip.open(local_sifts_filepath, 'wb') as local_sifts_file:
                local_sifts_file.write(page)

    # ======
    # From PDB file, get EXPRESSION_SYSTEM and related fields, using Bio.PDB.PDBParser
    # ======

    db_chain_ids_lower = [chain_dict['chain_id'].lower() for chain_dict in chain_data]

    pdbparser = Bio.PDB.PDBParser(QUIET=True)
    with gzip.open(local_pdb_filepath) as local_pdb_file:
        pdbdata = pdbparser.get_structure(pdb_id, local_pdb_file)
    pdbheader = pdbparser.get_header()
    # Bio PDB compound structure: {'compound': {'1': {'chain': 'a, b'}}}
    pdb_compounds = pdbheader['compound']
    matching_pdb_compound_id = None
    try:
        for pdb_compound_id in pdb_compounds.keys():
            for pdb_chain_id in pdb_compounds[pdb_compound_id]['chain'].split(', '):
                if pdb_chain_id in db_chain_ids_lower:
                    matching_pdb_compound_id = pdb_compound_id
                    break
        assert matching_pdb_compound_id is not None
    except Exception as e:
        print 'ERROR for entry %s PDB %s. PDB header dict as parsed by BioPython follows:' % (entry_name, pdb_id)
        print pdbheader
        print traceback.format_exc()
        raise e

    expression_data = {}
    # Bio PDB source structure: {'source': {'1': {'expression_system': 'escherichia coli'}}}
    pdbexpression_data = pdbheader['source'][matching_pdb_compound_id]
    for key in pdbexpression_data.keys():
        if key[0:10] == 'expression':
            # Make expression data upper-case again. I think it looks better for single-case text.
            expression_data[key.upper()] = pdbexpression_data[key].upper()
            # expression_data_obj = models.PDBExpressionData(expression_data_type=key.upper(), expression_data_value=pdbexpression_data[key].upper(), pdb=pdb_row)
            # db.session.add(expression_data_obj)

    # ======
    # Iterate through chains in PDBRow and extract sequence data from SIFTS file, and add to database
    # ======

    results = {'pdb_row_id': pdb_row_id, 'expression_data': expression_data, 'chain_dicts': {}}

    for chain_dict in chain_data:
        chain_row_id = chain_dict['chain_row_id']
        chain_id = chain_dict['chain_id']
        logger.debug(entry_name, ac, pdb_id, chain_id)
        pdb_chain_dict = extract_sifts_seq(local_sifts_filepath, ac, entry_name, pdb_id, chain_id, seq)
        results['chain_dicts'][chain_row_id] = pdb_chain_dict

    return results


def retrieve_sifts(pdb_id):
    '''Retrieves a SIFTS .xml file, given a PDB ID. Works by modifying the PDBe download URL.
    File is retrieved compressed and returned uncompressed.
    Also removes annoying namespace stuff.
    '''
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


def extract_sifts_seq(sifts_filepath, uniprot_ac, uniprot_entry_name, pdb_id, chain_id, uniprot_sequence):
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
        logger.info('PDB %s chain %s picked up from UniProt entry %s %s. Non-matching UniProtAC in sifts: %s. This chain will be deleted.' % (pdb_id, chain_id, uniprot_entry_name, uniprot_ac, sifts_uniprot_ac))
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
            print 'ERROR: UniProt crossref not found for conflicting residue!', uniprot_ac, pdb_id, chain_id, r.attrib
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
        # TODO need to generalize this. Shift to manual_overrides.yaml or do something else? In the short-term, perhaps just skip these PDBs?
        # Some pdb resids are e.g. '464A'
        if pdb_resid.isdigit() == False:
            if pdb_id in ['1O6L','2JDO','2JDR','2UW9','2X39','2XH5']: # These pdbs include three residues with pdb resids 464A, 464B, 464C, (all with UniProt crossrefs) then continues from 465. We will change this so that the pdb resids continue to iterate
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
            print 'Problem converting pdb_resid into int.', uniprot_ac, pdb_id, chain_id, pdb_resid
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

    #print e, uniprot_ac, pdb_id, chain_id
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
