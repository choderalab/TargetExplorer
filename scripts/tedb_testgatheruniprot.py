#!/usr/bin/env python
from app import db, models

print models.DBEntry.query.all()

existing_entries = models.DBEntry.query.all()
for a in existing_entries:
    print 'Deleting', a
    db.session.delete(a)

existing_uniprot = models.UniProt.query.all()
for a in existing_uniprot:
    print 'Deleting', a
    db.session.delete(a)

existing_pdb = models.PDB.query.all()
for a in existing_pdb:
    print 'Deleting', a
    db.session.delete(a)

db.session.commit()

entry = models.DBEntry()
db.session.add(entry)
uniprot = models.UniProt(ac='P00519', entry_name='ABL1_HUMAN', taxonid=9606, dbentry=entry)
db.session.add(uniprot)

entry = models.DBEntry()
db.session.add(entry)
uniprot = models.UniProt(ac='P12931', entry_name='SRC_HUMAN', taxonid=9606, dbentry=entry)
db.session.add(uniprot)

entry = models.DBEntry()
db.session.add(entry)
uniprot = models.UniProt(ac='P00533', entry_name='EGFR_HUMAN', taxonid=9606, dbentry=entry)
db.session.add(uniprot)

db.session.commit()
