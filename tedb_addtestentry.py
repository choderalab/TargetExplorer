from app import db, models

entry = models.DBEntry(ac='P00519', entry_name='ABL1_HUMAN')
db.session.add(entry)
entry = models.DBEntry(ac='P12931', entry_name='SRC_HUMAN')
db.session.add(entry)
entry = models.DBEntry(ac='P00533', entry_name='EGFR_HUMAN')
db.session.add(entry)

db.session.commit()
