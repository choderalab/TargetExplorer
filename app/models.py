from app import db

class DBEntry(db.Model):
    id = db.Column(db.Integer, primary_key=True)
    ac = db.Column(db.String(120), index=True, unique=True)
    entry_name = db.Column(db.String(120), index=True, unique=True)

    def __repr__(self):
        return '<Entry AC %r entry_name %r>' % (self.ac, self.entry_name)
