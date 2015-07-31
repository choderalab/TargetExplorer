from targetexplorer.flaskapp import db
from flask_sqlalchemy import _BoundDeclarativeMeta


frontend2backend_mappings = {
    'ac': ['UniProtEntry', 'ac'],
    'name': ['UniProtEntry', 'entry_name'],
    'domain': ['UniProtDomain', 'domain_id'],
    'npdbs': ['DBEntry', 'npdbs'],
    'ndomains': ['DBEntry', 'ndomains'],
    'nisoforms': ['DBEntry', 'nisoforms'],
    'nfunctions': ['DBEntry', 'nfunctions'],
    'ndiseaseassociations': ['DBEntry', 'ndisease_associations'],
    'npubs': ['DBEntry', 'npubs'],
    'nbioassays': ['DBEntry', 'nbioassays'],
    'family': ['UniProtEntry', 'family'],
    'species': ['UniProtEntry', 'taxon_name_common'],
    'domain_length': ['UniProtDomain', 'length'],
    'subcellular_location': ['UniProtSubcellularLocation', 'subcellular_location'],
    'pseudodomain': ['UniProtDomain', 'is_pseudodomain'],
}


class CrawlData(db.Model):
    __tablename__ = 'crawldata'
    id = db.Column(db.Integer, primary_key=True)
    current_crawl_number = db.Column(db.Integer)
    safe_crawl_number = db.Column(db.Integer)
    safe_crawl_datestamp = db.Column(db.DateTime)
    def __repr__(self):
        return '<DB CrawlData current crawl number {} safe crawl number {}>'.format(
            self.current_crawl_number, self.safe_crawl_number
        )


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
        return '<DB DateStamps crawl number {}>'.format(self.crawl_number)


class DBEntry(db.Model):
    """
    Essentially a "top-level" table for uniting all data types.
    Represents a "biological entity" - incorporating gene, transcript, and protein data.
    """
    __tablename__ = 'db_entries'
    id = db.Column(db.Integer, primary_key=True)
    crawl_number = db.Column(db.Integer)
    npdbs = db.Column(db.Integer)
    ndomains = db.Column(db.Integer)
    nisoforms = db.Column(db.Integer)
    nfunctions = db.Column(db.Integer)
    ndisease_associations = db.Column(db.Integer)
    npubs = db.Column(db.Integer)
    nbioassays = db.Column(db.Integer)
    uniprot = db.relationship('UniProtEntry', backref='db_entry', lazy='dynamic')
    uniprot_domains = db.relationship('UniProtDomain', backref='db_entry', lazy='dynamic')
    uniprot_gene_names = db.relationship('UniProtGeneName', backref='db_entry', lazy='dynamic')
    uniprot_isoforms = db.relationship('UniProtIsoform', backref='db_entry', lazy='dynamic')
    uniprot_functions = db.relationship('UniProtFunction', backref='db_entry', lazy='dynamic')
    uniprot_disease_associations = db.relationship('UniProtDiseaseAssociation', backref='db_entry', lazy='dynamic')
    uniprot_subcellular_locations = db.relationship('UniProtSubcellularLocation', backref='db_entry', lazy='dynamic')
    pdbs = db.relationship('PDBEntry', backref='db_entry', lazy='dynamic')
    ncbi_gene_entries = db.relationship('NCBIGeneEntry', backref='db_entry', lazy='dynamic')
    ensembl_genes = db.relationship('EnsemblGene', backref='db_entry', lazy='dynamic')
    hgnc_entries = db.relationship('HGNCEntry', backref='db_entry', lazy='dynamic')
    bindingdb_bioassays = db.relationship('BindingDBBioassay', backref='db_entry', lazy='dynamic')
    cbioportal_mutations = db.relationship('CbioportalMutation', backref='db_entry', lazy='dynamic')
    def __repr__(self):
        return '<DBEntry {}>'.format(self.id)


class UniProtEntry(db.Model):
    __tablename__ = 'uniprot_entries'
    id = db.Column(db.Integer, primary_key=True)
    crawl_number = db.Column(db.Integer)
    ac = db.Column(db.String(64))
    entry_name = db.Column(db.String(64))
    family = db.Column(db.String(64))
    recommended_name = db.Column(db.Text)
    ncbi_taxon_id = db.Column(db.String(64))
    taxon_name_scientific = db.Column(db.String(120))
    taxon_name_common = db.Column(db.String(120))
    lineage = db.Column(db.Text)   # ascending comma-separated values
    last_uniprot_update = db.Column(db.String(64))
    isoforms = db.relationship('UniProtIsoform', backref='uniprot_entry', lazy='dynamic')
    domains = db.relationship('UniProtDomain', backref='uniprot_entry', lazy='dynamic')
    gene_names = db.relationship('UniProtGeneName', backref='uniprot_entry', lazy='dynamic')
    functions = db.relationship('UniProtFunction', backref='uniprot_entry', lazy='dynamic')
    disease_associations = db.relationship('UniProtDiseaseAssociation', backref='uniprot_entry', lazy='dynamic')
    subcellular_locations = db.relationship('UniProtSubcellularLocation', backref='uniprot_entry', lazy='dynamic')
    db_entry_id = db.Column(db.Integer, db.ForeignKey('db_entries.id'))
    def __repr__(self):
        return '<UniProtEntry AC {} entry_name {}>'.format(self.ac, self.entry_name)


class UniProtGeneName(db.Model):
    __tablename__ = 'uniprot_gene_names'
    id = db.Column(db.Integer, primary_key=True)
    crawl_number = db.Column(db.Integer)
    gene_name = db.Column(db.String(64))
    gene_name_type = db.Column(db.String(64))
    db_entry_id = db.Column(db.Integer, db.ForeignKey('db_entries.id'))
    uniprot_id = db.Column(db.Integer, db.ForeignKey('uniprot_entries.id'))
    def __repr__(self):
        return '<UniProtGeneName {}>'.format(self.gene_name)


class UniProtIsoform(db.Model):
    __tablename__ = 'uniprot_isoforms'
    id = db.Column(db.Integer, primary_key=True)
    crawl_number = db.Column(db.Integer)
    ac = db.Column(db.String(64))
    is_canonical = db.Column(db.Boolean)
    length = db.Column(db.Integer)
    mass = db.Column(db.Integer)
    date_modified = db.Column(db.String(64))
    version = db.Column(db.Integer)
    sequence = db.Column(db.Text)
    notes = db.relationship('UniProtIsoformNote', backref='uniprot_isoform', lazy='dynamic')
    ensembl_transcripts = db.relationship('EnsemblTranscript', backref='uniprot_isoform', lazy='dynamic')
    db_entry_id = db.Column(db.Integer, db.ForeignKey('db_entries.id'))
    uniprot_id = db.Column(db.Integer, db.ForeignKey('uniprot_entries.id'))
    def __repr__(self):
        return '<UniProtIsoform {} canonical: {}>'.format(self.ac, self.is_canonical)


class UniProtIsoformNote(db.Model):
    __tablename__ = 'uniprot_isoform_notes'
    id = db.Column(db.Integer, primary_key=True)
    crawl_number = db.Column(db.Integer)
    note = db.Column(db.Text)
    uniprot_isoform_id = db.Column(db.Integer, db.ForeignKey('uniprot_isoforms.id'))
    def __repr__(self):
        return '<UniProtIsoformNote {}>'.format(self.note)


class UniProtDomain(db.Model):
    __tablename__ = 'uniprot_domains'
    id = db.Column(db.Integer, primary_key=True)
    crawl_number = db.Column(db.Integer)
    domain_id = db.Column(db.Integer)
    target_id = db.Column(db.String(64))   # ABL1_HUMAN_D0 (Protein kinase)
    is_target_domain = db.Column(db.Boolean)
    description = db.Column(db.Text())
    begin = db.Column(db.Integer)
    end = db.Column(db.Integer)
    length = db.Column(db.Integer)
    sequence = db.Column(db.Text)
    is_pseudodomain = db.Column(db.Boolean)
    pseudodomain_notes = db.Column(db.Text)
    pdb_chains = db.relationship('PDBChain', backref='uniprot_domain', lazy='dynamic')
    cbioportal_mutations = db.relationship('CbioportalMutation', backref='uniprot_domain', lazy='dynamic')
    db_entry_id = db.Column(db.Integer, db.ForeignKey('db_entries.id'))
    uniprot_id = db.Column(db.Integer, db.ForeignKey('uniprot_entries.id'))
    def __repr__(self):
        return '<UniProtDomain description {}>'.format(self.description)


class UniProtFunction(db.Model):
    __tablename__ = 'uniprot_functions'
    id = db.Column(db.Integer, primary_key=True)
    crawl_number = db.Column(db.Integer)
    function = db.Column(db.Text)
    db_entry_id = db.Column(db.Integer, db.ForeignKey('db_entries.id'))
    uniprot_id = db.Column(db.Integer, db.ForeignKey('uniprot_entries.id'))
    def __repr__(self):
        return '<UniProtFunction {}>'.format(self.function)


class UniProtDiseaseAssociation(db.Model):
    __tablename__ = 'uniprot_disease_associations'
    id = db.Column(db.Integer, primary_key=True)
    crawl_number = db.Column(db.Integer)
    disease_association = db.Column(db.Text)
    db_entry_id = db.Column(db.Integer, db.ForeignKey('db_entries.id'))
    uniprot_id = db.Column(db.Integer, db.ForeignKey('uniprot_entries.id'))
    def __repr__(self):
        return '<UniProtDiseaseAssociation {}>'.format(self.disease_association)


class UniProtSubcellularLocation(db.Model):
    __tablename__ = 'uniprot_subcellular_locations'
    id = db.Column(db.Integer, primary_key=True)
    crawl_number = db.Column(db.Integer)
    subcellular_location = db.Column(db.Text)
    db_entry_id = db.Column(db.Integer, db.ForeignKey('db_entries.id'))
    uniprot_id = db.Column(db.Integer, db.ForeignKey('uniprot_entries.id'))
    def __repr__(self):
        return '<UniProtSubcellularLocation {}>'.format(self.subcellular_location)


class PDBEntry(db.Model):
    __tablename__ = 'pdb_entries'
    id = db.Column(db.Integer, primary_key=True)
    crawl_number = db.Column(db.Integer)
    pdb_id = db.Column(db.String(64))
    method = db.Column(db.Text)
    resolution = db.Column(db.Float)
    chains = db.relationship('PDBChain', backref='pdb_entry', lazy='dynamic')
    expression_data = db.relationship('PDBExpressionData', backref='pdb_entry', lazy='dynamic')
    db_entry_id = db.Column(db.Integer, db.ForeignKey('db_entries.id'))
    def __repr__(self):
        return '<PDBEntry ID {}>'.format(self.pdb_id)


class PDBChain(db.Model):
    __tablename__ = 'pdb_chains'
    id = db.Column(db.Integer, primary_key=True)
    crawl_number = db.Column(db.Integer)
    chain_id = db.Column(db.String(64))
    begin = db.Column(db.Integer)
    end = db.Column(db.Integer)
    experimental_seq = db.Column(db.Text)
    experimental_seq_aln_conflicts = db.Column(db.Text)
    experimental_seq_len = db.Column(db.Integer)
    observed_seq_aln_exp = db.Column(db.Text)
    observed_seq_aln = db.Column(db.Text)
    observed_ss_aln = db.Column(db.Text)
    pdb_entry_id = db.Column(db.Integer, db.ForeignKey('pdb_entries.id'))
    uniprot_domain_id = db.Column(db.Integer, db.ForeignKey('uniprot_domains.id'))
    def __repr__(self):
        return '<PDBChain ID {}>'.format(self.chain_id)


class PDBExpressionData(db.Model):
    __tablename__ = 'pdb_expression_data'
    id = db.Column(db.Integer, primary_key=True)
    crawl_number = db.Column(db.Integer)
    expression_data_type = db.Column(db.String(64))
    expression_data_value = db.Column(db.Text)
    pdb_entry_id = db.Column(db.Integer, db.ForeignKey('pdb_entries.id'))
    def __repr__(self):
        return '<PDBExpressionData type {}>'.format(self.expression_data_type)


class NCBIGeneEntry(db.Model):
    __tablename__ = 'ncbi_gene_entries'
    id = db.Column(db.Integer, primary_key=True)
    crawl_number = db.Column(db.Integer)
    gene_id = db.Column(db.Integer)
    publications = db.relationship('NCBIGenePublication', backref='ncbi_gene_entry', lazy='dynamic')
    db_entry_id = db.Column(db.Integer, db.ForeignKey('db_entries.id'))
    def __repr__(self):
        return '<NCBIGeneEntry ID {}>'.format(self.gene_id)


class NCBIGenePublication(db.Model):
    __tablename__ = 'ncbi_gene_publication'
    id = db.Column(db.Integer, primary_key=True)
    crawl_number = db.Column(db.Integer)
    pmid = db.Column(db.Integer)
    ncbi_gene_entry_id = db.Column(db.Integer, db.ForeignKey('ncbi_gene_entries.id'))
    def __repr__(self):
        return '<NCBIGenePublication PMID {}>'.format(self.pmid)


class EnsemblGene(db.Model):
    __tablename__ = 'ensembl_genes'
    id = db.Column(db.Integer, primary_key=True)
    crawl_number = db.Column(db.Integer)
    gene_id = db.Column(db.String(64))
    ensembl_transcripts = db.relationship('EnsemblTranscript', backref='ensembl_gene', lazy='dynamic')
    ensembl_proteins = db.relationship('EnsemblProtein', backref='ensembl_gene', lazy='dynamic')
    db_entry_id = db.Column(db.Integer, db.ForeignKey('db_entries.id'))
    def __repr__(self):
        return '<EnsemblGene ID {}>'.format(self.gene_id)


class EnsemblTranscript(db.Model):
    __tablename__ = 'ensembl_transcripts'
    id = db.Column(db.Integer, primary_key=True)
    crawl_number = db.Column(db.Integer)
    transcript_id = db.Column(db.String(64))
    ensembl_proteins = db.relationship('EnsemblProtein', backref='ensembl_transcript', lazy='dynamic')
    ensembl_gene_id = db.Column(db.Integer, db.ForeignKey('ensembl_genes.id'))
    uniprot_isoform_id = db.Column(db.Integer, db.ForeignKey('uniprot_isoforms.id'))
    def __repr__(self):
        return '<EnsemblTranscript ID {}>'.format(self.transcript_id)


class EnsemblProtein(db.Model):
    __tablename__ = 'ensembl_proteins'
    id = db.Column(db.Integer, primary_key=True)
    crawl_number = db.Column(db.Integer)
    protein_id = db.Column(db.String(64))
    ensembl_gene_id = db.Column(db.Integer, db.ForeignKey('ensembl_genes.id'))
    ensembl_transcript_id = db.Column(db.Integer, db.ForeignKey('ensembl_transcripts.id'))
    def __repr__(self):
        return '<EnsemblProtein ID {}>'.format(self.protein_id)


class HGNCEntry(db.Model):
    __tablename__ = 'hgnc_entries'
    id = db.Column(db.Integer, primary_key=True)
    crawl_number = db.Column(db.Integer)
    gene_id = db.Column(db.String(64))
    approved_symbol = db.Column(db.String(64))
    db_entry_id = db.Column(db.Integer, db.ForeignKey('db_entries.id'))
    def __repr__(self):
        return '<HGNCEntry approved_symbol {}>'.format(self.approved_symbol)


class BindingDBBioassay(db.Model):
    __tablename__ = 'bindingdb_bioassays'
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
    db_entry_id = db.Column(db.Integer, db.ForeignKey('db_entries.id'))
    def __repr__(self):
        return '<BindingDBBioassay source {}>'.format(self.bindingdb_source)


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
    __tablename__ = 'cbioportal_cases'
    id = db.Column(db.Integer, primary_key=True)
    crawl_number = db.Column(db.Integer)
    study = db.Column(db.Text)
    case_id = db.Column(db.Text)
    mutations = db.relationship('CbioportalMutation', backref='cbioportal_case', lazy='dynamic')
    def __repr__(self):
        return '<CbioportalCase ID {0}>'.format(self.case_id)


class CbioportalMutation(db.Model):
    __tablename__ = 'cbioportal_mutations'
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
    db_entry_id = db.Column(db.Integer, db.ForeignKey('db_entries.id'))
    uniprot_domain_id = db.Column(db.Integer, db.ForeignKey('uniprot_domains.id'))
    cbioportal_case_id = db.Column(db.Integer, db.ForeignKey('cbioportal_cases.id'))
    def __repr__(self):
        return '<CbioportalMutation ID type {} aa_change {} in_uniprot_domain {}>'.format(
            self.type, self.oncotator_aa_pos, self.in_uniprot_domain
        )


_module_local_names = [key for key in locals().keys()]
table_class_names = [
    key for key in _module_local_names if isinstance(locals()[key], _BoundDeclarativeMeta)
]
