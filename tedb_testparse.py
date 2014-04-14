from app import db, models

entries = models.DBEntry.query.all()
for entry in entries:
    print entry
