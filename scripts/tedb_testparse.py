#!/usr/bin/env python
from app import db, models

entries = models.DBEntry.query.all()
print 'Entries'
for entry in entries:
    print entry, entry.uniprot.all(), entry.pdbs.all()

print ''

uniprot = models.UniProt.query.all()
print 'UniProt'
for u in uniprot:
    print u, u.taxonid

print ''

pdbs = models.PDB.query.all()
print 'PDBs'
for pdb in pdbs:
    print pdb
