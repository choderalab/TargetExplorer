__author__ = 'isikm'

import datetime
from targetexplorer.flaskapp import models, db
from bioservices import ChEMBL


class GatherChemblTarget(object):
    def __init__(self, run_main=True):

        if run_main:
            self.main()

    def main(self):

        # get current crawl number
        crawldata_row = models.CrawlData.query.first()
        current_crawl_number = crawldata_row.current_crawl_number
        print 'Current crawl number: %d' % current_crawl_number

        # Get list of uniprot accession numbers that are in the rows of current crawl number
        uniprot_acs = [str(uniprot_ac[0]) for uniprot_ac in models.UniProt.query.filter_by(crawl_number=current_crawl_number).values(models.UniProt.ac)]

        # Iterate through Uniprot ACs
        for uniprot_ac in uniprot_acs:

            # get ChEMBL target id using bioservices
            chembl=ChEMBL(verbose=False)
            target_info = chembl.get_target_by_uniprotId(uniprot_ac)
            target_chembl_id=target_info['chemblId']
            print "ChEMBL ID of target protein:", target_chembl_id

            target = models.ChemblTarget(crawl_number=current_crawl_number, target_chembl_id=target_chembl_id)
            db.session.add(target)

        # Updata datestamp for Chembl
        now = datetime.datetime.utcnow()
        current_crawl_datestamp_row = models.DateStamps.query.filter_by(crawl_number=current_crawl_number).first()
        current_crawl_datestamp_row.chembl_datestamp = now

        db.session.commit()
        print 'Done.'
