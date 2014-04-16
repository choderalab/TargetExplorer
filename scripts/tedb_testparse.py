#!/usr/bin/env python
from app import db, models

dbversion = models.Version.query.all()[0]
print dbversion
print 'gather_uniprot last run:', dbversion.uniprotdatestamp

entries = models.DBEntry.query.all()
print '= %d entries =' % len(entries)
for e in range(10):
    if e >= len(entries):
        continue
    entry = entries[e]
    uniprot = entry.uniprot.all()[0]
    print entry, uniprot,
    print 'ndomains:', entry.uniprotdomains.count(),
    print 'nPDBs:', entry.pdbs.count()
    # print db.session.query(models.UniProtGeneName).filter_by(dbentry=entry, gene_name_type='primary').first()
    # print db.session.query(models.UniProtDiseaseAssociation).filter_by(dbentry=entry).all()

    # isoforms = db.session.query(models.UniProtIsoform).filter_by(dbentry=entry).all()
    # for isoform in isoforms:
    #     if isoform.canonical:
    #         print isoform
    #     else:
    #         notes = db.session.query(models.UniProtIsoformNote).filter_by(uniprotisoform=isoform).all()
    #         print isoform, notes

print ''

# print db.session.query(models.UniProtIsoform).all()

# uniprot = models.UniProt.query.all()
# print 'UniProt'
# for u in uniprot:
#     print u, u.taxonid
#
# print ''
#
# pdbs = models.PDB.query.all()
# print 'PDBs'
# for pdb in pdbs:
#     print pdb
