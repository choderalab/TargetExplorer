from targetexplorer.flaskapp import db
from flask_sqlalchemy import _BoundDeclarativeMeta


frontend2backend_mappings = {
    'ac': ['UniProt', 'ac'],
    'name': ['UniProt', 'entry_name'],
    'target': ['UniProtDomain', 'targetid'],
    'npdbs': ['DBEntry', 'npdbs'],
    'ndomains': ['DBEntry', 'ndomains'],
    'nisoforms': ['DBEntry', 'nisoforms'],
    'nfunctions': ['DBEntry', 'nfunctions'],
    'ndiseaseassociations': ['DBEntry', 'ndiseaseassociations'],
    'npubs': ['DBEntry', 'npubs'],
    'nbioassays': ['DBEntry', 'nbioassays'],
    'family': ['UniProt', 'family'],
    'species': ['UniProt', 'taxon_name_common'],
    'domain_length': ['UniProtDomain', 'length'],
    'subcellular_location': ['UniProtSubcellularLocation', 'subcellular_location'],
    # TODO 'pseudodomain': ['UniProtDomain', 'pseudodomain'],
}


class CrawlData(db.Model):
    __tablename__ = 'crawldata'
    id = db.Column(db.Integer, primary_key=True)
    current_crawl_number = db.Column(db.Integer)
    safe_crawl_number = db.Column(db.Integer)
    safe_crawl_datestamp = db.Column(db.DateTime)
    def __repr__(self):
        return '<DB CrawlData current crawl number %d safe crawl number %d>' % (self.current_crawl_number, self.safe_crawl_number)

class DateStamps(db.Model):
    __tablename__ = 'datestamps'
    id = db.Column(db.Integer, primary_key=True)
    crawl_number = db.Column(db.Integer)
    uniprot_datestamp = db.Column(db.DateTime)
    pdb_datestamp = db.Column(db.DateTime)
    ncbi_gene_datestamp = db.Column(db.DateTime)
    bindingdb_datestamp = db.Column(db.DateTime)
    cbioportal_datestamp = db.Column(db.DateTime)
    commit_datestamp = db.Column(db.DateTime)
    def __repr__(self):
        return '<DB DateStamps crawl number %d>' % self.crawl_number

class DBEntry(db.Model):
    __tablename__ = 'dbentry'
    id = db.Column(db.Integer, primary_key=True)
    crawl_number = db.Column(db.Integer)
    npdbs = db.Column(db.Integer)
    ndomains = db.Column(db.Integer)
    nisoforms = db.Column(db.Integer)
    nfunctions = db.Column(db.Integer)
    ndiseaseassociations = db.Column(db.Integer)
    npubs = db.Column(db.Integer)
    nbioassays = db.Column(db.Integer)
    uniprot = db.relationship('UniProt', backref='dbentry', lazy='dynamic')
    uniprotdomains = db.relationship('UniProtDomain', backref='dbentry', lazy='dynamic')
    uniprotgenenames = db.relationship('UniProtGeneName', backref='dbentry', lazy='dynamic')
    uniprotisoforms = db.relationship('UniProtIsoform', backref='dbentry', lazy='dynamic')
    uniprotfunctions = db.relationship('UniProtFunction', backref='dbentry', lazy='dynamic')
    uniprotdiseaseassociations = db.relationship('UniProtDiseaseAssociation', backref='dbentry', lazy='dynamic')
    uniprotsubcellularlocations = db.relationship('UniProtSubcellularLocation', backref='dbentry', lazy='dynamic')
    pdbs = db.relationship('PDB', backref='dbentry', lazy='dynamic')
    ncbi_gene_entries = db.relationship('NCBIGeneEntry', backref='dbentry', lazy='dynamic')
    ensembl_genes = db.relationship('EnsemblGene', backref='dbentry', lazy='dynamic')
    hgnc_entries = db.relationship('HGNCEntry', backref='dbentry', lazy='dynamic')
    bindingdb_bioassays = db.relationship('BindingDBBioassay', backref='dbentry', lazy='dynamic')
    cbioportal_mutations = db.relationship('CbioportalMutation', backref='dbentry', lazy='dynamic')
    def __repr__(self):
        return '<DBEntry %d>' % self.id

class UniProt(db.Model):
    __tablename__ = 'uniprot'
    id = db.Column(db.Integer, primary_key=True)
    crawl_number = db.Column(db.Integer)
    ac = db.Column(db.String(64))
    entry_name = db.Column(db.String(64))
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
    crawl_number = db.Column(db.Integer)
    gene_name = db.Column(db.String(64))
    gene_name_type = db.Column(db.String(64))
    dbentry_id = db.Column(db.Integer, db.ForeignKey('dbentry.id'))
    uniprot_id = db.Column(db.Integer, db.ForeignKey('uniprot.id'))
    def __repr__(self):
        return '<UniProtGeneName %r>' % self.gene_name

class UniProtIsoform(db.Model):
    __tablename__ = 'uniprotisoform'
    id = db.Column(db.Integer, primary_key=True)
    crawl_number = db.Column(db.Integer)
    ac = db.Column(db.String(64))
    canonical = db.Column(db.Boolean)
    length = db.Column(db.Integer)
    mass = db.Column(db.Integer)
    date_modified = db.Column(db.String(64))
    version = db.Column(db.Integer)
    sequence = db.Column(db.Text)
    notes = db.relationship('UniProtIsoformNote', backref='uniprotisoform', lazy='dynamic')
    ensembl_transcripts = db.relationship('EnsemblTranscript', backref='uniprotisoform', lazy='dynamic')
    dbentry_id = db.Column(db.Integer, db.ForeignKey('dbentry.id'))
    uniprot_id = db.Column(db.Integer, db.ForeignKey('uniprot.id'))
    def __repr__(self):
        return '<UniProtIsoform %r canonical: %r>' % (self.ac, self.canonical)

class UniProtIsoformNote(db.Model):
    __tablename__ = 'uniprotisoformnote'
    id = db.Column(db.Integer, primary_key=True)
    crawl_number = db.Column(db.Integer)
    note = db.Column(db.Text)
    uniprotisoformid = db.Column(db.Integer, db.ForeignKey('uniprotisoform.id'))
    def __repr__(self):
        return '<UniProtIsoformNote %r>' % self.note

class UniProtDomain(db.Model):
    __tablename__ = 'uniprotdomain'
    id = db.Column(db.Integer, primary_key=True)
    crawl_number = db.Column(db.Integer)
    targetid = db.Column(db.String(64))
    description = db.Column(db.Text())
    begin = db.Column(db.Integer)
    end = db.Column(db.Integer)
    length = db.Column(db.Integer)
    sequence = db.Column(db.Text)
    # pseudodomain = db.Column(db.Boolean)
    cbioportal_mutations = db.relationship('CbioportalMutation', backref='uniprot_domain', lazy='dynamic')
    dbentry_id = db.Column(db.Integer, db.ForeignKey('dbentry.id'))
    uniprot_id = db.Column(db.Integer, db.ForeignKey('uniprot.id'))
    def __repr__(self):
        return '<UniProtDomain targetid %r>' % self.targetid

class UniProtFunction(db.Model):
    __tablename__ = 'uniprotfunction'
    id = db.Column(db.Integer, primary_key=True)
    crawl_number = db.Column(db.Integer)
    function = db.Column(db.Text)
    dbentry_id = db.Column(db.Integer, db.ForeignKey('dbentry.id'))
    uniprot_id = db.Column(db.Integer, db.ForeignKey('uniprot.id'))
    def __repr__(self):
        return '<UniProtFunction %r>' % self.function

class UniProtDiseaseAssociation(db.Model):
    __tablename__ = 'uniprotdiseaseassociation'
    id = db.Column(db.Integer, primary_key=True)
    crawl_number = db.Column(db.Integer)
    disease_association = db.Column(db.Text)
    dbentry_id = db.Column(db.Integer, db.ForeignKey('dbentry.id'))
    uniprot_id = db.Column(db.Integer, db.ForeignKey('uniprot.id'))
    def __repr__(self):
        return '<UniProtDiseaseAssociation %r>' % self.disease_association

class UniProtSubcellularLocation(db.Model):
    __tablename__ = 'uniprotsubcellularlocation'
    id = db.Column(db.Integer, primary_key=True)
    crawl_number = db.Column(db.Integer)
    subcellular_location = db.Column(db.Text)
    dbentry_id = db.Column(db.Integer, db.ForeignKey('dbentry.id'))
    uniprot_id = db.Column(db.Integer, db.ForeignKey('uniprot.id'))
    def __repr__(self):
        return '<UniProtSubcellularLocation %r>' % self.subcellular_location

class PDB(db.Model):
    __tablename__ = 'pdb'
    id = db.Column(db.Integer, primary_key=True)
    crawl_number = db.Column(db.Integer)
    pdbid = db.Column(db.String(64))
    method = db.Column(db.Text)
    resolution = db.Column(db.Float)
    chains = db.relationship('PDBChain', backref='pdb', lazy='dynamic')
    expression_data = db.relationship('PDBExpressionData', backref='pdb', lazy='dynamic')
    dbentry_id = db.Column(db.Integer, db.ForeignKey('dbentry.id'))
    def __repr__(self):
        return '<PDB ID %r>' % self.pdbid

class PDBChain(db.Model):
    __tablename__ = 'pdbchain'
    id = db.Column(db.Integer, primary_key=True)
    crawl_number = db.Column(db.Integer)
    chain_id = db.Column(db.String(64))
    domain_id = db.Column(db.Integer)
    begin = db.Column(db.Integer)
    end = db.Column(db.Integer)
    experimental_seq = db.Column(db.Text)
    experimental_seq_aln_conflicts = db.Column(db.Text)
    experimental_seq_len = db.Column(db.Integer)
    observed_seq_aln_exp = db.Column(db.Text)
    observed_seq_aln = db.Column(db.Text)
    observed_ss_aln = db.Column(db.Text)
    pdb_id = db.Column(db.Integer, db.ForeignKey('pdb.id'))
    def __repr__(self):
        return '<PDBChain ID %r>' % self.chain_id

class PDBExpressionData(db.Model):
    __tablename__ = 'pdbexpressiondata'
    id = db.Column(db.Integer, primary_key=True)
    crawl_number = db.Column(db.Integer)
    expression_data_type = db.Column(db.String(64))
    expression_data_value = db.Column(db.Text)
    pdb_id = db.Column(db.Integer, db.ForeignKey('pdb.id'))
    def __repr__(self):
        return '<PDBExpressionData id %r>' % self.id

class NCBIGeneEntry(db.Model):
    __tablename__ = 'ncbi_gene_entry'
    id = db.Column(db.Integer, primary_key=True)
    crawl_number = db.Column(db.Integer)
    gene_id = db.Column(db.Integer)
    publications = db.relationship('NCBIGenePublication', backref='ncbi_gene_entry', lazy='dynamic')
    dbentry_id = db.Column(db.Integer, db.ForeignKey('dbentry.id'))
    def __repr__(self):
        return '<NCBIGeneEntry ID %r>' % self.gene_id

class NCBIGenePublication(db.Model):
    __tablename__ = 'ncbi_gene_publication'
    id = db.Column(db.Integer, primary_key=True)
    crawl_number = db.Column(db.Integer)
    pmid = db.Column(db.Integer)
    ncbi_gene_entry_id = db.Column(db.Integer, db.ForeignKey('ncbi_gene_entry.id'))
    def __repr__(self):
        return '<NCBIGenePublication PMID %r>' % self.pmid

class EnsemblGene(db.Model):
    __tablename__ = 'ensembl_gene'
    id = db.Column(db.Integer, primary_key=True)
    crawl_number = db.Column(db.Integer)
    gene_id = db.Column(db.String(64))
    ensembl_transcripts = db.relationship('EnsemblTranscript', backref='ensembl_gene', lazy='dynamic')
    ensembl_proteins = db.relationship('EnsemblProtein', backref='ensembl_gene', lazy='dynamic')
    dbentry_id = db.Column(db.Integer, db.ForeignKey('dbentry.id'))
    def __repr__(self):
        return '<EnsemblGene ID %r>' % self.gene_id

class EnsemblTranscript(db.Model):
    __tablename__ = 'ensembl_transcript'
    id = db.Column(db.Integer, primary_key=True)
    crawl_number = db.Column(db.Integer)
    transcript_id = db.Column(db.String(64))
    ensembl_proteins = db.relationship('EnsemblProtein', backref='ensembl_transcript', lazy='dynamic')
    ensembl_gene_id = db.Column(db.Integer, db.ForeignKey('ensembl_gene.id'))
    uniprot_isoform_id = db.Column(db.Integer, db.ForeignKey('uniprotisoform.id'))
    def __repr__(self):
        return '<EnsemblTranscript ID %r>' % self.transcript_id

class EnsemblProtein(db.Model):
    __tablename__ = 'ensembl_protein'
    id = db.Column(db.Integer, primary_key=True)
    crawl_number = db.Column(db.Integer)
    protein_id = db.Column(db.String(64))
    ensembl_gene_id = db.Column(db.Integer, db.ForeignKey('ensembl_gene.id'))
    ensembl_transcript_id = db.Column(db.Integer, db.ForeignKey('ensembl_transcript.id'))
    def __repr__(self):
        return '<EnsemblProtein ID %r>' % self.protein_id

class HGNCEntry(db.Model):
    __tablename__ = 'hgnc_entry'
    id = db.Column(db.Integer, primary_key=True)
    crawl_number = db.Column(db.Integer)
    gene_id = db.Column(db.String(64))
    approved_symbol = db.Column(db.String(64))
    dbentry_id = db.Column(db.Integer, db.ForeignKey('dbentry.id'))
    def __repr__(self):
        return '<HGNCEntry approved_symbol %r>' % self.approved_symbol

class BindingDBBioassay(db.Model):
    __tablename__ = 'bindingdb_bioassay'
    id = db.Column(db.Integer, primary_key=True)
    crawl_number = db.Column(db.Integer)
    bindingdb_source = db.Column(db.Text)
    doi = db.Column(db.String(64))
    pmid = db.Column(db.Integer)
    temperature = db.Column(db.String(64))
    ph = db.Column(db.String(64))
    target_name = db.Column(db.Text)
    ligand_bindingdb_id = db.Column(db.Integer)
    ligand_chembl_id = db.Column(db.String(64))
    ligand_smiles_string = db.Column(db.Text)
    ligand_zinc_id = db.Column(db.String(64))

    ki = db.Column(db.String(64))
    ic50 = db.Column(db.String(64))
    kd = db.Column(db.String(64))
    ec50 = db.Column(db.String(64))
    kon = db.Column(db.String(64))
    koff = db.Column(db.String(64))

    # measurements = db.relationship('BindingDBMeasurement', backref='bindingdb_bioassay', lazy='dynamic')
    dbentry_id = db.Column(db.Integer, db.ForeignKey('dbentry.id'))
    def __repr__(self):
        return '<BindingDB bioassay>'

# class BindingDBMeasurement(db.Model):
#     __tablename__ = 'bindingdb_measurement'
#     id = db.Column(db.Integer, primary_key=True)
#     crawl_number = db.Column(db.Integer)
#     measurement_type = db.Column(db.String(64))  # {Kd, Ki, IC50, EC50, koff, kon}
#     measurement_value = db.Column(db.String(64))
#     bindingdb_bioassay_id = db.Column(db.Integer, db.ForeignKey('bindingdb_bioassay.id'))
#     def __repr__(self):
#         return '<BindingDB measurement type: %r value: %r>' % (self.measurement_type, self.measurement_value)

class CbioportalCase(db.Model):
    __tablename__ = 'cbioportal_case'
    id = db.Column(db.Integer, primary_key=True)
    crawl_number = db.Column(db.Integer)
    study = db.Column(db.Text)
    case_id = db.Column(db.Text)
    mutations = db.relationship('CbioportalMutation', backref='cbioportal_case', lazy='dynamic')
    def __repr__(self):
        return '<CbioportalCase ID {0}>'.format(self.case_id)

class CbioportalMutation(db.Model):
    __tablename__ = 'cbioportal_mutation'
    id = db.Column(db.Integer, primary_key=True)
    crawl_number = db.Column(db.Integer)
    type = db.Column(db.Text)
    cbioportal_aa_change_string = db.Column(db.Text)
    mutation_origin = db.Column(db.Text)
    validation_status = db.Column(db.Text)
    functional_impact_score = db.Column(db.Text)
    chromosome_index = db.Column(db.Integer)
    chromosome_startpos = db.Column(db.Integer)
    chromosome_endpos = db.Column(db.Integer)
    reference_dna_allele = db.Column(db.Text)
    variant_dna_allele = db.Column(db.Text)
    oncotator_aa_pos = db.Column(db.Integer)
    oncotator_reference_aa = db.Column(db.Text)
    oncotator_variant_aa = db.Column(db.Text)
    oncotator_ensembl_transcript_id = db.Column(db.Text)
    in_uniprot_domain = db.Column(db.Boolean)
    dbentry_id = db.Column(db.Integer, db.ForeignKey('dbentry.id'))
    uniprot_domain_id = db.Column(db.Integer, db.ForeignKey('uniprotdomain.id'))
    cbioportal_case_id = db.Column(db.Integer, db.ForeignKey('cbioportal_case.id'))
    def __repr__(self):
        return '<CbioportalMutation ID {0} type {1} aa_change {2} in_uniprot_domain {3}>'.format(self.id, self.type, self.oncotator_aa_pos, self.in_uniprot_domain)


_module_local_names = [key for key in locals().keys()]
table_class_names = [
    key for key in _module_local_names if isinstance(locals()[key], _BoundDeclarativeMeta)
]
