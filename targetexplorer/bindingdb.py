import urllib2
import subprocess
import os
import datetime
from targetexplorer.flaskapp import models, db
from targetexplorer.core import external_data_dirpath, logger

bindingdb_data_dir = os.path.join(external_data_dirpath, 'BindingDB')
bindingdb_all_data_filepath = os.path.join(bindingdb_data_dir, 'BindingDB_All.tab')
bindingdb_matches_filepath = os.path.join(bindingdb_data_dir, 'bindingdb-matches.tab')


class GatherBindingDB(object):
    def __init__(self,
                 use_existing_bindingdb_data=False,
                 run_main=True
                 ):
        self.use_existing_bindingdb_data = use_existing_bindingdb_data
        self.now = datetime.datetime.utcnow()
        crawldata_row = models.CrawlData.query.first()
        self.current_crawl_number = crawldata_row.current_crawl_number
        if run_main:
            self.setup()
            self.get_bindingdb_data_file()
            self.get_uniprot_acs_from_db()
            extracted_bindingdb_data = extract_bindingdb_data(
                bindingdb_all_data_filepath,
                self.db_uniprot_acs
            )
            self.create_db_rows(extracted_bindingdb_data)
            self.commit_to_db()

    def setup(self):
        if not os.path.exists(bindingdb_data_dir):
            os.mkdir(bindingdb_data_dir)

    def get_bindingdb_data_file(self):
        # Unless use_existing_bindingdb_data is set to True, retrieve a new file from BindingDB
        if os.path.exists(bindingdb_all_data_filepath) and self.use_existing_bindingdb_data:
            logger.info('BindingDB data file found at: {0}'.format(bindingdb_all_data_filepath))
        else:
            logger.info('Retrieving new BindingDB data file from BindingDB server...')
            retrieve_all_BindingDB_data(bindingdb_all_data_filepath, decompress=False)

    def get_uniprot_acs_from_db(self):
        self.db_uniprot_acs = [
            value_tuple[0] for value_tuple in models.UniProt.query.filter_by(
                crawl_number=self.current_crawl_number
            ).values(models.UniProt.ac)
        ]

    def create_db_rows(self, extracted_bindingdb_data):
        # for bindingdb_data_tuple in extracted_bindingdb_data:
        for ac, bioassays_data in extracted_bindingdb_data.iteritems():
            # ac, bioassays_data = bindingdb_data_tuple
            db_uniprot_row = models.UniProt.query.filter_by(
                crawl_number=self.current_crawl_number, ac=ac
            ).first()
            dbentry_id = db_uniprot_row.dbentry_id
            dbentry_row = models.DBEntry.query.filter_by(id=dbentry_id).first()
            dbentry_row.nbioassays = len(bioassays_data)

            for bioassay_data in bioassays_data:
                bindingdb_bioassay_obj = models.BindingDBBioassay(
                    crawl_number=self.current_crawl_number,
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

    def commit_to_db(self):
        current_crawl_datestamp_row = models.DateStamps.query.filter_by(
            crawl_number=self.current_crawl_number
        ).first()
        current_crawl_datestamp_row.bindingdb_datestamp = self.now
        db.session.commit()
        print 'Done.'


def retrieve_all_BindingDB_data(bindingdb_all_data_filepath, decompress=True):
    """
    Retrieves all BindingDB data in compressed tab-separated format, then decompresses and writes to file.
    Uses gunzip for decompression.
    """
    print 'Downloading BindingDB_All.tab.gz...'
    url = 'http://bindingdb.org/bind/downloads/BindingDB_All_2015m5.tsv.zip'
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

def get_acs(line):
    words = line.split('\t')
    UniProtACs = words[20].split(' ')
    return UniProtACs

def get_bioassay_data(bioassay_data_line_tsplit):
    SMILES_string = bioassay_data_line_tsplit[0]
    BindingDB_monomerID = bioassay_data_line_tsplit[1]
    ChEMBL_ID = bioassay_data_line_tsplit[8]
    zinc_id = bioassay_data_line_tsplit[13]
    data_origin = bioassay_data_line_tsplit[15]
    target_biomolecule = bioassay_data_line_tsplit[16]
    # target_source_organism = bioassay_data_line_tsplit[17]
    # target_sequence = bioassay_data_line_tsplit[18]
    # pdbID = bioassay_data_line_tsplit[19]
    # UniProtAC = bioassay_data_line_tsplit[20]
    Ki = bioassay_data_line_tsplit[22]
    IC50 = bioassay_data_line_tsplit[23]
    Kd = bioassay_data_line_tsplit[24]
    EC50 = bioassay_data_line_tsplit[25]
    kon = bioassay_data_line_tsplit[26]
    koff = bioassay_data_line_tsplit[27]
    PMID = bioassay_data_line_tsplit[28]
    DOI = bioassay_data_line_tsplit[29]
    pH = bioassay_data_line_tsplit[32]
    temperature = bioassay_data_line_tsplit[33]

    bioassay_data = dict()

    bioassay_data['ligand_SMILES_string'] = SMILES_string
    bioassay_data['ligand_BindingDB_ID'] = BindingDB_monomerID
    bioassay_data['ligand_ChEMBL_ID'] = ChEMBL_ID
    bioassay_data['ligand_zinc_id'] = zinc_id
    bioassay_data['BindingDB_source'] = data_origin
    bioassay_data['target_name'] = target_biomolecule
    #bioassay_data['target_sequence'] = target_sequence
    #bioassay_data['target_source_organism'] = target_source_organism

    bioassay_data['Ki'] = Ki
    bioassay_data['IC50'] = IC50
    bioassay_data['Kd'] = Kd
    bioassay_data['EC50'] = EC50
    bioassay_data['kon'] = kon
    bioassay_data['koff'] = koff
    bioassay_data['pH'] = pH
    bioassay_data['temperature'] = temperature

    bioassay_data['PMID'] = PMID
    bioassay_data['DOI'] = DOI

    return bioassay_data


def extract_bindingdb_data(
        bindingdb_data_filepath,
        uniprot_acs
    ):
    bindingdb_data_dict = dict()
    with open(bindingdb_data_filepath) as bindingdb_data_file:
        first_line = True
        for line in bindingdb_data_file:
            if first_line:
                first_line = False
                continue
            words = line.split('\t')
            ac_field = words[20]
            ac = ac_field.split(' ')[0]
            if ac in uniprot_acs:
                if ac not in bindingdb_data_dict:
                    bindingdb_data_dict[ac] = []
                bindingdb_data_dict[ac].append(
                    get_bioassay_data(words)   # returns a dict
                )

    return bindingdb_data_dict
