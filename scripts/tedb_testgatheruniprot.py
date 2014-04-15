#!/usr/bin/env python
from app import db, models

entry = models.DBEntry()
db.session.add(entry)
uniprot = models.UniProt(ac='P00519', entry_name='ABL1_HUMAN', taxonid=9606, dbentry=entry)
# uniprot = models.UniProt(ac='P00519', entry_name='ABL1_HUMAN', dbentry=entry)
db.session.add(uniprot)

entry = models.DBEntry()
db.session.add(entry)
uniprot = models.UniProt(ac='P12931', entry_name='SRC_HUMAN', taxonid=9606, dbentry=entry)
# uniprot = models.UniProt(ac='P12931', entry_name='SRC_HUMAN', dbentry=entry)
db.session.add(uniprot)

entry = models.DBEntry()
db.session.add(entry)
uniprot = models.UniProt(ac='P00533', entry_name='EGFR_HUMAN', taxonid=9606, dbentry=entry)
# uniprot = models.UniProt(ac='P00533', entry_name='EGFR_HUMAN', dbentry=entry)
db.session.add(uniprot)

db.session.commit()
