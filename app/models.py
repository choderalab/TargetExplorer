from app import db

class Version(db.Model):
    __tablename__ = 'version'
    id = db.Column(db.Integer, primary_key=True)
    versionid = db.Column(db.Integer)
    uniprotdatestamp = db.Column(db.DateTime)
    def __repr__(self):
        if self.versionid:
            return '<DB Version %d>' % (self.versionid)
        else:
            return '<DB Version None>'

class DBEntry(db.Model):
    __tablename__ = 'dbentry'
    id = db.Column(db.Integer, primary_key=True)
    uniprot = db.relationship('UniProt', backref='dbentry', lazy='dynamic')
    uniprotdomains = db.relationship('UniProtDomain', backref='dbentry', lazy='dynamic')
    uniprotgenenames = db.relationship('UniProtGeneName', backref='dbentry', lazy='dynamic')
    pdbs = db.relationship('PDB', backref='dbentry', lazy='dynamic')
    def __repr__(self):
        return '<DBEntry %d>' % (self.id)

class UniProt(db.Model):
    __tablename__ = 'uniprot'
    id = db.Column(db.Integer, primary_key=True)
    ac = db.Column(db.String(64), unique=True)
    entry_name = db.Column(db.String(64), unique=True)
    family = db.Column(db.String(64))
    taxonid = db.Column(db.String(64))
    recommended_name = db.Column(db.Text)
    last_uniprot_update = db.Column(db.String(64))
    domains = db.relationship('UniProtDomain', backref='uniprot_entry', lazy='dynamic')
    genenames = db.relationship('UniProtGeneName', backref='uniprot_entry', lazy='dynamic')
    dbentry_id = db.Column(db.Integer, db.ForeignKey('dbentry.id'))
    def __repr__(self):
        return '<UniProt AC %r entry_name %r>' % (self.ac, self.entry_name)

class UniProtGeneName(db.Model):
    __tablename__ = 'uniprotgenename'
    id = db.Column(db.Integer, primary_key=True)
    gene_name = db.Column(db.String(64))
    gene_name_type = db.Column(db.String(64))
    dbentry_id = db.Column(db.Integer, db.ForeignKey('dbentry.id'))
    uniprotid = db.Column(db.Integer, db.ForeignKey('uniprot.id'))
    def __repr__(self):
        return '<UniProtGeneName %r>' % (self.gene_name)

class UniProtDomain(db.Model):
    __tablename__ = 'uniprotdomain'
    id = db.Column(db.Integer, primary_key=True)
    targetid = db.Column(db.String(64), unique=True)
    description = db.Column(db.Text())
    begin = db.Column(db.Integer)
    end = db.Column(db.Integer)
    length = db.Column(db.Integer)
    sequence = db.Column(db.Text())
    dbentry_id = db.Column(db.Integer, db.ForeignKey('dbentry.id'))
    uniprotid = db.Column(db.Integer, db.ForeignKey('uniprot.id'))
    def __repr__(self):
        return '<UniProtDomain targetid %r>' % (self.targetid)

class PDB(db.Model):
    __tablename__ = 'pdb'
    id = db.Column(db.Integer, primary_key=True)
    pdbid = db.Column(db.String(64))
    dbentry_id = db.Column(db.Integer, db.ForeignKey('dbentry.id'))
    def __repr__(self):
        return '<PDB ID %r>' % (self.pdbid)

