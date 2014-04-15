#!/usr/bin/env python
from app import db, models

entry = models.DBEntry()
db.session.add(entry)
uniprot = models.UniProt(ac='P00519', entry_name='ABL1_HUMAN', dbentry=entry)
db.session.add(uniprot)
pdbids = ['1AB2', '1ABL', '1AWO']
[db.session.add(models.PDB(pdbid=pdbid, dbentry=entry)) for pdbid in pdbids]

entry = models.DBEntry()
db.session.add(entry)
uniprot = models.UniProt(ac='P12931', entry_name='SRC_HUMAN', dbentry=entry)
db.session.add(uniprot)
pdbids = ['1A07', '1A08', '1A09']
[db.session.add(models.PDB(pdbid=pdbid, dbentry=entry)) for pdbid in pdbids]

entry = models.DBEntry()
db.session.add(entry)
uniprot = models.UniProt(ac='P00533', entry_name='EGFR_HUMAN', dbentry=entry)
db.session.add(uniprot)
pdbids = ['1DNQ', '1DNR', '1IVO']
[db.session.add(models.PDB(pdbid=pdbid, dbentry=entry)) for pdbid in pdbids]

db.session.commit()
