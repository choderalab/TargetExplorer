from app import db

class DBEntry(db.Model):
    __tablename__ = 'dbentry'
    id = db.Column(db.Integer, primary_key=True)
    uniprot = db.relationship('UniProt', backref='dbentry', lazy='dynamic')
    pdbs = db.relationship('PDB', backref='dbentry', lazy='dynamic')
    def __repr__(self):
        return '<DBEntry %d>' % (self.id)

class UniProt(db.Model):
    __tablename__ = 'uniprot'
    id = db.Column(db.Integer, primary_key=True)
    ac = db.Column(db.String(64), unique=True)
    entry_name = db.Column(db.String(64), unique=True)
    taxonid = db.Column(db.String(64))
    dbentry_id = db.Column(db.Integer, db.ForeignKey('dbentry.id'))
    def __repr__(self):
        return '<UniProtData AC %r entry_name %r>' % (self.ac, self.entry_name)

class PDB(db.Model):
    __tablename__ = 'pdb'
    id = db.Column(db.Integer, primary_key=True)
    pdbid = db.Column(db.String(64), unique=True)
    dbentry_id = db.Column(db.Integer, db.ForeignKey('dbentry.id'))
    def __repr__(self):
        return '<PDB ID %r>' % (self.pdbid)

