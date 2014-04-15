#!/usr/bin/env python
from app import db, models

dbentries = models.DBEntry.query.all()

for dbentry in dbentries:
    uniprot = dbentry.uniprot.all()[0]
    if uniprot.ac == 'P00519':
        pdbids = ['1AB2', '1ABL', '1AWO']
        [db.session.add(models.PDB(pdbid=pdbid, dbentry=dbentry)) for pdbid in pdbids]
    elif uniprot.ac == 'P12931':
        pdbids = ['1A07', '1A08', '1A09']
        [db.session.add(models.PDB(pdbid=pdbid, dbentry=dbentry)) for pdbid in pdbids]
    elif uniprot.ac == 'P00533':
        pdbids = ['1DNQ', '1DNR', '1IVO']
        [db.session.add(models.PDB(pdbid=pdbid, dbentry=dbentry)) for pdbid in pdbids]

db.session.commit()
