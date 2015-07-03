import urllib2
import subprocess
import os
import datetime
from targetexplorer.flaskapp import models, db
from multiprocessing import Pool


class GatherBindingDB(object):
    def __init__(self,
                 use_existing_bindingdb_data=False,
                 grep_path=False,
                 run_main=True
                 ):
        self.use_existing_bindingdb_data = use_existing_bindingdb_data
        self.grep_path = grep_path
        if run_main:
            self.main()

    def main(self):
        external_data_dir = os.path.join('external-data')
        bindingdb_data_dir = os.path.join(external_data_dir, 'BindingDB')

        if not os.path.exists(bindingdb_data_dir):
            os.mkdir(bindingdb_data_dir)

        bindingdb_all_data_filepath = os.path.join(bindingdb_data_dir, 'BindingDB_All.tab')

        now = datetime.datetime.utcnow()

        # get current crawl number
        crawldata_row = models.CrawlData.query.first()
        current_crawl_number = crawldata_row.current_crawl_number
        print 'Current crawl number: %d' % current_crawl_number

        # =================
        # Download new BindingDB data
        # =================

        # Unless use_existing_bindingdb_data is set to True, retrieve a new file from BindingDB
        if os.path.exists(bindingdb_all_data_filepath) and self.use_existing_bindingdb_data:
            print 'BindingDB data file found at:', bindingdb_all_data_filepath
        else:
            print 'Retrieving new BindingDB data file from BindingDB server...'
            retrieve_all_BindingDB_data(bindingdb_all_data_filepath, decompress=False)



        # clab.DB.retrieve_external_data(DB_script_ID, forcedl=force_bindingdb_all_download)

        # =================
        # Use grep to extract information for each entry in the DB
        # =================

        # Overview:
        #   * first use grep -E to extract all lines matching any kinase uniprot AC
        #   * using multiprocessing for each kinase:
        #       * grep through to create data file for each kinase
        #       * parse each kinase data file and add to kinDB
        #       * delete individual kinase data files

        db_uniprot_acs = [value_tuple[0] for value_tuple in models.UniProt.query.filter_by(crawl_number=current_crawl_number).values(models.UniProt.ac)]

        # db_uniprot_acs = ['P00519']




        regex_search = '|'.join(db_uniprot_acs)
        bindingdb_matches_filepath = os.path.join(bindingdb_data_dir, 'tmp-bindingdb-matches.tab')

        # TODO use proper logic here
        if True:
            print 'First pass through BindingDB data with grep... - extracting all lines which match UniProt ACs in the DB (performance can be variable - testing has indicated approx. 1 hr on a 2012 MBP, but only 2 mins or so on a Linux machine; probably due to different GNU grep versions)'
            with open(bindingdb_all_data_filepath, 'r') as bindingdb_all_data_file:
                #subprocess.call('%(self.grep_path)s -E "%(regex_search)s" %(bindingdb_all_data_filepath)s > %(bindingdb_matches_filepath)s' % vars(), shell=True)
                subprocess.call('%s -E "%s" %s > %s' % (self.grep_path, regex_search, bindingdb_all_data_filepath, bindingdb_matches_filepath), shell=True)

        def extract_bindingdb(input_data):
            AC, grep_path, bindingdb_matches_filepath = input_data
            print AC
            DB_entry_bindingdb_data_filepath = os.path.join(bindingdb_data_dir, AC + '.tab')
            subprocess.call('%s -E "%s" %s > %s' % (grep_path, AC, bindingdb_matches_filepath, DB_entry_bindingdb_data_filepath), shell=True)
            bioassays_data = []
            with open(DB_entry_bindingdb_data_filepath, 'r') as bindingdb_file:
                for line in bindingdb_file:
                    returnedACs = get_ACs(line)
                    for returnedAC in returnedACs:
                        if returnedAC == AC:
                            bioassay_data = get_bioassay_data(line) # returns a dict
                            bioassays_data.append( bioassay_data )
            os.remove(DB_entry_bindingdb_data_filepath)
            return (AC, bioassays_data)

        # pool = Pool()
        input_data = [(AC, self.grep_path, bindingdb_matches_filepath) for AC in db_uniprot_acs]
        # TODO get this working in parallel. Currently gives some pickling error when using multiprocessing
        results = map(extract_bindingdb, input_data)

        for result in results:
            ac, bioassays_data = result
            db_uniprot_row = models.UniProt.query.filter_by(crawl_number=current_crawl_number, ac=ac).first()
            dbentry_id = db_uniprot_row.dbentry_id
            dbentry_row = models.DBEntry.query.filter_by(id=dbentry_id).first()
            dbentry_row.nbioassays = len(bioassays_data)

            for bioassay_data in bioassays_data:
                bindingdb_bioassay_obj = models.BindingDBBioassay(
                    crawl_number=current_crawl_number,
                    bindingdb_source=bioassay_data['BindingDB_source'],
                    doi=bioassay_data['DOI'],
                    pmid=bioassay_data['PMID'],
                    temperature=bioassay_data['temperature'],
                    ph=bioassay_data['pH'],
                    target_name=bioassay_data['target_name'],
                    ligand_bindingdb_id=bioassay_data['ligand_BindingDB_ID'],
                    ligand_chembl_id=bioassay_data['ligand_ChEMBL_ID'],
                    ligand_smiles_string=bioassay_data['ligand_SMILES_string'],
                    ligand_zinc_id=bioassay_data['ligand_zinc_id'],

                    ki=bioassay_data['Ki'],
                    ic50=bioassay_data['IC50'],
                    kd=bioassay_data['Kd'],
                    ec50=bioassay_data['EC50'],
                    kon=bioassay_data['kon'],
                    koff=bioassay_data['koff'],

                    dbentry=dbentry_row,
                )
                db.session.add(bindingdb_bioassay_obj)

        # ==============================================================
        # Update db BindingDB datestamp
        # ==============================================================

        current_crawl_datestamp_row = models.DateStamps.query.filter_by(crawl_number=current_crawl_number).first()
        current_crawl_datestamp_row.bindingdb_datestamp = now
        db.session.commit()
        print 'Done.'





def retrieve_all_BindingDB_data(bindingdb_all_data_filepath, decompress=True):
    '''
    Retrieves all BindingDB data in compressed tab-separated format, then decompresses and writes to file.
    Uses gunzip for decompression.
    '''
    print 'Downloading BindingDB_All.tab.gz...'
    url = 'http://bindingdb.org/bind/downloads/BindingDB_All.tab.gz'
    response = urllib2.urlopen(url)
    page = response.read(1000000000)
    with open(bindingdb_all_data_filepath+'.gz', 'wb') as ofile:
        ofile.write(page)
    print 'Finished downloading BindingDB_All.tab.gz'

    # decompress
    if decompress:
        print 'Decompressing BindingDB_All.tab.gz...'
        subprocess.call('gunzip -f %s' % bindingdb_all_data_filepath+'.gz', shell=True)
        print 'Finished decompressing BindingDB_All.tab.gz'

def get_ACs(line):
    words = line.split('\t')
    UniProtACs = words[20].split(' ')
    return UniProtACs

def get_bioassay_data(line):
    words = line.split('\t')
    SMILES_string = words[0]
    BindingDB_monomerID = words[1]
    ChEMBL_ID = words[8]
    zinc_id = words[13]
    data_origin = words[15]
    target_biomolecule = words[16]
    #target_source_organism = words[17]
    #target_sequence = words[18]
    #pdbID = words[19]
    #UniProtAC = words[20]
    Ki = words[22]
    IC50 = words[23]
    Kd = words[24]
    EC50 = words[25]
    kon = words[26]
    koff = words[27]
    PMID = words[28]
    DOI = words[29]
    pH = words[32]
    temperature = words[33]

    bioassay_data = {}

    bioassay_data['ligand_SMILES_string'] = SMILES_string
    bioassay_data['ligand_BindingDB_ID'] = BindingDB_monomerID
    bioassay_data['ligand_ChEMBL_ID'] = ChEMBL_ID
    bioassay_data['ligand_zinc_id'] = zinc_id
    bioassay_data['BindingDB_source'] = data_origin
    bioassay_data['target_name'] = target_biomolecule
    #bioassay_data['target_sequence'] = target_sequence
    #bioassay_data['target_source_organism'] = target_source_organism

    bioassay_data['Ki'] = Ki # XXX remove
    bioassay_data['IC50'] = IC50
    bioassay_data['Kd'] = Kd
    bioassay_data['EC50'] = EC50
    bioassay_data['kon'] = kon
    bioassay_data['koff'] = koff
    bioassay_data['pH'] = pH
    bioassay_data['temperature'] = temperature
    # if Ki != 'n/a':
    #     bioassay_data['Ki'] = Ki
    # if IC50 != 'n/a':
    #     bioassay_data['IC50'] = IC50
    # if Kd != 'n/a':
    #     bioassay_data['Kd'] = Kd
    # if EC50 != 'n/a':
    #     bioassay_data['EC50'] = EC50
    # if kon != 'n/a':
    #     bioassay_data['kon'] = kon
    # if koff != 'n/a':
    #     bioassay_data['koff'] = koff
    # if pH != 'n/a':
    #     bioassay_data['pH'] = pH
    # if temperature != 'n/a':
    #     bioassay_data['temperature'] = temperature

    bioassay_data['PMID'] = PMID
    bioassay_data['DOI'] = DOI

    return bioassay_data

