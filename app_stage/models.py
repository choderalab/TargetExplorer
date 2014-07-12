from app_stage import db

table_class_names = [
    'Version',
    'DBEntry',
    'UniProt',
    'UniProtGeneName',
    'UniProtIsoform',
    'UniProtIsoformNote',
    'UniProtDomain',
    'UniProtFunction',
    'UniProtDiseaseAssociation',
    'UniProtSubcellularLocation',
    'PDB',
    'NCBIGeneEntry',
    'EnsemblGeneEntry',
    'HGNCEntry']

frontend2backend_mappings = {
    'npdbs': ['DBEntry', 'npdbs'],
    'ndomains': ['DBEntry', 'ndomains'],
    'nisoforms': ['DBEntry', 'nisoforms'],
    'nfunctions': ['DBEntry', 'nfunctions'],
    'ndiseaseassociations': ['DBEntry', 'ndiseaseassociations'],
    'family': ['UniProt', 'family'],
    'species': ['UniProt', 'taxon_name_common'],
    'domain_length': ['UniProtDomain', 'length'],
    'subcellular_location': ['UniProtSubcellularLocation', 'subcellular_location'],
   #  'pseudodomain': ['UniProtDomain', 'pseudodomain'],
}


class Version(db.Model):
    __tablename__ = 'version'
    id = db.Column(db.Integer, primary_key=True)
    version_id = db.Column(db.Integer)
    uniprot_datestamp = db.Column(db.DateTime)
    pdb_datestamp = db.Column(db.DateTime)
    def __repr__(self):
        if self.version_id:
            return '<DB Version %d>' % (self.version_id)
        else:
            return '<DB Version None>'

class DBEntry(db.Model):
    __tablename__ = 'dbentry'
    id = db.Column(db.Integer, primary_key=True)
    npdbs = db.Column(db.Integer)
    ndomains = db.Column(db.Integer)
    nisoforms = db.Column(db.Integer)
    nfunctions = db.Column(db.Integer)
    ndiseaseassociations = db.Column(db.Integer)
    uniprot = db.relationship('UniProt', backref='dbentry', lazy='dynamic')
    uniprotdomains = db.relationship('UniProtDomain', backref='dbentry', lazy='dynamic')
    uniprotgenenames = db.relationship('UniProtGeneName', backref='dbentry', lazy='dynamic')
    uniprotisoforms = db.relationship('UniProtIsoform', backref='dbentry', lazy='dynamic')
    uniprotfunctions = db.relationship('UniProtFunction', backref='dbentry', lazy='dynamic')
    uniprotdiseaseassociations = db.relationship('UniProtDiseaseAssociation', backref='dbentry', lazy='dynamic')
    uniprotsubcellularlocations = db.relationship('UniProtSubcellularLocation', backref='dbentry', lazy='dynamic')
    pdbs = db.relationship('PDB', backref='dbentry', lazy='dynamic')
    ncbi_gene_entries = db.relationship('NCBIGeneEntry', backref='dbentry', lazy='dynamic')
    ensembl_gene_entries = db.relationship('EnsemblGeneEntry', backref='dbentry', lazy='dynamic')
    hgnc_entries = db.relationship('HGNCEntry', backref='dbentry', lazy='dynamic')
    def __repr__(self):
        return '<DBEntry %d>' % (self.id)

class UniProt(db.Model):
    __tablename__ = 'uniprot'
    id = db.Column(db.Integer, primary_key=True)
    ac = db.Column(db.String(64), unique=True)
    entry_name = db.Column(db.String(64), unique=True)
    family = db.Column(db.String(64))
    recommended_name = db.Column(db.Text)
    ncbi_taxonid = db.Column(db.String(64))
    taxon_name_scientific = db.Column(db.String(120))
    taxon_name_common = db.Column(db.String(120))
    lineage = db.Column(db.Text) # ascending comma-separated values
    last_uniprot_update = db.Column(db.String(64))
    isoforms = db.relationship('UniProtIsoform', backref='uniprot_entry', lazy='dynamic')
    domains = db.relationship('UniProtDomain', backref='uniprot_entry', lazy='dynamic')
    genenames = db.relationship('UniProtGeneName', backref='uniprot_entry', lazy='dynamic')
    functions = db.relationship('UniProtFunction', backref='uniprot_entry', lazy='dynamic')
    disease_associations = db.relationship('UniProtDiseaseAssociation', backref='uniprot_entry', lazy='dynamic')
    subcellular_locations = db.relationship('UniProtSubcellularLocation', backref='uniprot_entry', lazy='dynamic')
    dbentry_id = db.Column(db.Integer, db.ForeignKey('dbentry.id'))
    def __repr__(self):
        return '<UniProt AC %r entry_name %r>' % (self.ac, self.entry_name)

class UniProtGeneName(db.Model):
    __tablename__ = 'uniprotgenename'
    id = db.Column(db.Integer, primary_key=True)
    gene_name = db.Column(db.String(64))
    gene_name_type = db.Column(db.String(64))
    dbentry_id = db.Column(db.Integer, db.ForeignKey('dbentry.id'))
    uniprot_id = db.Column(db.Integer, db.ForeignKey('uniprot.id'))
    def __repr__(self):
        return '<UniProtGeneName %r>' % (self.gene_name)

class UniProtIsoform(db.Model):
    __tablename__ = 'uniprotisoform'
    id = db.Column(db.Integer, primary_key=True)
    ac = db.Column(db.String(64))
    canonical = db.Column(db.Boolean)
    length = db.Column(db.Integer)
    mass = db.Column(db.Integer)
    date_modified = db.Column(db.String(64))
    version = db.Column(db.Integer)
    sequence = db.Column(db.Text)
    notes = db.relationship('UniProtIsoformNote', backref='uniprotisoform', lazy='dynamic')
    dbentry_id = db.Column(db.Integer, db.ForeignKey('dbentry.id'))
    uniprot_id = db.Column(db.Integer, db.ForeignKey('uniprot.id'))
    def __repr__(self):
        return '<UniProtIsoform %r canonical: %r>' % (self.ac, self.canonical)

class UniProtIsoformNote(db.Model):
    __tablename__ = 'uniprotisoformnote'
    id = db.Column(db.Integer, primary_key=True)
    note = db.Column(db.Text)
    uniprotisoformid = db.Column(db.Integer, db.ForeignKey('uniprotisoform.id'))
    def __repr__(self):
        return '<UniProtIsoformNote %r>' % (self.note)

class UniProtDomain(db.Model):
    __tablename__ = 'uniprotdomain'
    id = db.Column(db.Integer, primary_key=True)
    targetid = db.Column(db.String(64), unique=True)
    description = db.Column(db.Text())
    begin = db.Column(db.Integer)
    end = db.Column(db.Integer)
    length = db.Column(db.Integer)
    sequence = db.Column(db.Text)
    # pseudodomain = db.Column(db.Boolean)
    dbentry_id = db.Column(db.Integer, db.ForeignKey('dbentry.id'))
    uniprot_id = db.Column(db.Integer, db.ForeignKey('uniprot.id'))
    def __repr__(self):
        return '<UniProtDomain targetid %r>' % (self.targetid)

class UniProtFunction(db.Model):
    __tablename__ = 'uniprotfunction'
    id = db.Column(db.Integer, primary_key=True)
    function = db.Column(db.Text)
    dbentry_id = db.Column(db.Integer, db.ForeignKey('dbentry.id'))
    uniprot_id = db.Column(db.Integer, db.ForeignKey('uniprot.id'))
    def __repr__(self):
        return '<UniProtFunction %r>' % (self.function)

class UniProtDiseaseAssociation(db.Model):
    __tablename__ = 'uniprotdiseaseassociation'
    id = db.Column(db.Integer, primary_key=True)
    disease_association = db.Column(db.Text)
    dbentry_id = db.Column(db.Integer, db.ForeignKey('dbentry.id'))
    uniprot_id = db.Column(db.Integer, db.ForeignKey('uniprot.id'))
    def __repr__(self):
        return '<UniProtDiseaseAssociation %r>' % (self.disease_association)

class UniProtSubcellularLocation(db.Model):
    __tablename__ = 'uniprotsubcellularlocation'
    id = db.Column(db.Integer, primary_key=True)
    subcellular_location = db.Column(db.Text)
    dbentry_id = db.Column(db.Integer, db.ForeignKey('dbentry.id'))
    uniprot_id = db.Column(db.Integer, db.ForeignKey('uniprot.id'))
    def __repr__(self):
        return '<UniProtSubcellularLocation %r>' % (self.subcellular_location)

class PDB(db.Model):
    __tablename__ = 'pdb'
    id = db.Column(db.Integer, primary_key=True)
    pdbid = db.Column(db.String(64))
    method = db.Column(db.Text)
    resolution = db.Column(db.Float)
    dbentry_id = db.Column(db.Integer, db.ForeignKey('dbentry.id'))
    def __repr__(self):
        return '<PDB ID %r>' % (self.pdbid)

# TODO
# class PDBChain(db.Model):

class NCBIGeneEntry(db.Model):
    __tablename__ = 'ncbi_gene_entry'
    id = db.Column(db.Integer, primary_key=True)
    gene_id = db.Column(db.Integer)
    dbentry_id = db.Column(db.Integer, db.ForeignKey('dbentry.id'))
    def __repr__(self):
        return '<NCBIGeneEntry ID %r>' % (self.gene_id)

class EnsemblGeneEntry(db.Model):
    __tablename__ = 'ensembl_gene_entry'
    id = db.Column(db.Integer, primary_key=True)
    gene_id = db.Column(db.String(64))
    dbentry_id = db.Column(db.Integer, db.ForeignKey('dbentry.id'))
    def __repr__(self):
        return '<EnsemblGeneEntry ID %r>' % (self.gene_id)

class HGNCEntry(db.Model):
    __tablename__ = 'hgnc_entry'
    id = db.Column(db.Integer, primary_key=True)
    gene_id = db.Column(db.String(64))
    approved_symbol = db.Column(db.String(64))
    dbentry_id = db.Column(db.Integer, db.ForeignKey('dbentry.id'))
    def __repr__(self):
        return '<HGNCEntry approved_symbol %r>' % (self.approved_symbol)

