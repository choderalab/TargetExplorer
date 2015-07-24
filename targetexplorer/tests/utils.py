import os
import shutil
import gzip
import functools
from nose import SkipTest
from contextlib import contextmanager
from targetexplorer.utils import installation_testdir_filepath, get_installed_resource_filepath
from targetexplorer.core import external_data_dirpath
from targetexplorer.flaskapp import app, db
from targetexplorer.initproject import InitProject
from targetexplorer.uniprot import GatherUniProt
from targetexplorer.pdb import GatherPDB
from targetexplorer.ncbi_gene import GatherNCBIGene
from targetexplorer.bindingdb import GatherBindingDB
from targetexplorer.cbioportal import GatherCbioportalData
from targetexplorer.commit import Commit


@contextmanager
def projecttest_context(set_up_project_stage='init'):
    """
    Context manager used within TargetExplorer testing framework.
    Uses the SetUpSampleProject class to set up a TargetExplorer database at various points of
    completion, within a temporary directory set up by the nose plugin `setup_tmp_db_plugin`.

    The following command will set up the database at the point after retrieval of UniProt data:

    >>> with projecttest_context(set_up_project_stage='uniprot'):
    >>>     ...

    See docs for SetUpSampleProject for other stages which can be set up.
    """
    # Set up
    with open(installation_testdir_filepath) as installation_testdir_file:
        temp_dir = installation_testdir_file.read()
    cwd = os.getcwd()
    os.chdir(temp_dir)

    set_up_sample_project(stage=set_up_project_stage)

    # Test is run at this point
    yield

    # Tear down
    db.drop_all()
    db.create_all()
    os.chdir(cwd)


def set_up_sample_project(stage='init'):
    sample_project = SetUpSampleProject()
    setup_method = getattr(sample_project, stage)
    setup_method()


class SetUpSampleProject(object):
    """
    Uses reference data files to set up a TargetExplorer database at various points of completion.
    The following commands will set up the database at the point after retrieval of UniProt data:

    >>> sample_project = SetUpSampleProject()
    >>> sample_project.uniprot()

    Other stages which can be set up are:
    * 'blank' - just a blank database file; no contained data
    * 'init' - as if user had run DoraInit.py
    * 'uniprot' - as if user had run DoraInit.py, then DoraGatherUniProt.py
    * 'pdb' - as if user had run above scripts, then DoraGatherPDB.py
    * 'ncbi_gene' - etc.
    * 'bindingdb' - etc.
    * 'cbioportal' - etc.
    * 'committed' - as if all gather scripts have been run, and DoraCommit.py
    """
    def __init__(self):
        self.structure_dirs = get_installed_resource_filepath(
            os.path.join('resources', 'structures')
        )
        self.uniprot_query = 'mnemonic:ABL1_HUMAN'
        self.uniprot_domain_regex = '^Protein kinase(?!; truncated)(?!; inactive)'

    def blank(self):
        pass

    def init(self):
        InitProject(
            db_name='test',
            project_path=os.getcwd(),
            uniprot_query=self.uniprot_query,
            uniprot_domain_regex=self.uniprot_domain_regex
        )
        self.copy_source_data_files()

    def copy_source_data_files(self):
        uniprot_ref_filepath = get_installed_resource_filepath(
            os.path.join('resources', 'uniprot-search-abl1.xml.gz')
        )
        if not os.path.exists(os.path.join(external_data_dirpath, 'UniProt')):
            os.mkdir(os.path.join(external_data_dirpath, 'UniProt'))
        with gzip.open(uniprot_ref_filepath) as uniprot_ref_file:
            with open(
                os.path.join(external_data_dirpath, 'UniProt', 'uniprot-search.xml'), 'w'
            ) as uniprot_test_file:
                uniprot_test_file.write(uniprot_ref_file.read())

        gene2pubmed_ref_filepath = get_installed_resource_filepath(
            os.path.join('resources', 'gene2pubmed-abl1.gz')
        )
        if not os.path.exists(os.path.join(external_data_dirpath, 'NCBI_Gene')):
            os.mkdir(os.path.join(external_data_dirpath, 'NCBI_Gene'))
        shutil.copy(
            gene2pubmed_ref_filepath,
            os.path.join(external_data_dirpath, 'NCBI_Gene', 'gene2pubmed.gz')
        )

        bindingdb_ref_filepath = get_installed_resource_filepath(
            os.path.join('resources', 'BindingDB-abl1.tab.gz')
        )
        if not os.path.exists(os.path.join(external_data_dirpath, 'BindingDB')):
            os.mkdir(os.path.join(external_data_dirpath, 'BindingDB'))
        with gzip.open(bindingdb_ref_filepath) as bindingdb_ref_file:
            with open(
                os.path.join(external_data_dirpath, 'BindingDB', 'BindingDB_All.tab'), 'w'
            ) as bindingdb_test_file:
                bindingdb_test_file.write(bindingdb_ref_file.read())

        cbioportal_ref_filepath = get_installed_resource_filepath(
            os.path.join('resources', 'cbioportal-mutations-abl1.xml.gz')
        )
        if not os.path.exists(os.path.join(external_data_dirpath, 'cBioPortal')):
            os.mkdir(os.path.join(external_data_dirpath, 'cBioPortal'))
        with gzip.open(cbioportal_ref_filepath) as cbioportal_ref_file:
            with open(
                os.path.join(external_data_dirpath, 'cBioPortal', 'cbioportal-mutations.xml'), 'w'
            ) as cbioportal_test_file:
                cbioportal_test_file.write(cbioportal_ref_file.read())

        oncotator_ref_filepath = get_installed_resource_filepath(
            os.path.join('resources', 'oncotator-data-abl1.json.gz')
        )
        shutil.copy(
            oncotator_ref_filepath,
            os.path.join(external_data_dirpath, 'cBioPortal', 'oncotator-data.json.gz')
        )

    def uniprot(self):
        self.init()
        GatherUniProt(
            uniprot_query=self.uniprot_query,
            uniprot_domain_regex=self.uniprot_domain_regex,
            use_existing_data=True,
        )

    def pdb(self):
        self.uniprot()
        GatherPDB(structure_dirs=self.structure_dirs)

    def ncbi_gene(self):
        self.pdb()
        GatherNCBIGene(use_existing_gene2pubmed=True)

    def bindingdb(self):
        self.ncbi_gene()
        GatherBindingDB(use_existing_bindingdb_data=True)

    def cbioportal(self):
        self.bindingdb()
        GatherCbioportalData(use_existing_cbioportal_data=True, use_existing_oncotator_data=True)

    def committed(self):
        self.cbioportal()
        Commit()


def expected_failure(test):
    @functools.wraps(test)
    def inner(*args, **kwargs):
        try:
            test(*args, **kwargs)
        except Exception:
            raise SkipTest
        else:
            raise AssertionError(
                'A failure was expected, but this test appeared to pass. You may want to remove the expected_failure decorator.'
            )
    return inner
