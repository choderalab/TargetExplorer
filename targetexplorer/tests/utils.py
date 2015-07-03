import os
import shutil
import gzip
import functools
from nose import SkipTest
from contextlib import contextmanager
from targetexplorer.utils import installation_testdir_filepath, get_installed_resource_filepath
from targetexplorer.core import external_data_dirpath
from targetexplorer.flaskapp import db
from targetexplorer.initproject import initialize_crawldata_and_datestamps
from targetexplorer.uniprot_tmp import GatherUniProt
from targetexplorer.pdb_tmp import GatherPDB
from targetexplorer.ncbi_gene_tmp import GatherNCBIGene
from targetexplorer.bindingdb_tmp import GatherBindingDB
from targetexplorer.cbioportal_tmp import GatherCbioportalData
from targetexplorer.commit import Commit


@contextmanager
def projecttest_context(set_up_project_stage='blankdb'):
    with open(installation_testdir_filepath) as installation_testdir_file:
        temp_dir = installation_testdir_file.read()
    cwd = os.getcwd()
    os.chdir(temp_dir)

    sample_project = SetUpSampleProject()
    setup_method = getattr(sample_project, set_up_project_stage)
    setup_method()
    yield

    db.drop_all()
    db.create_all()
    os.chdir(cwd)


class SetUpSampleProject(object):
    def __init__(self):
        self.uniprot_query ='mnemonic:ABL1_HUMAN'
        self.uniprot_domain_regex ='^Protein kinase(?!; truncated)(?!; inactive)'

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

    def blankdb(self):
        pass

    def init(self):
        initialize_crawldata_and_datestamps()

    def uniprot(self):
        self.init()
        GatherUniProt(
            uniprot_query=self.uniprot_query,
            uniprot_domain_regex=self.uniprot_domain_regex
        )

    def pdb(self):
        self.uniprot()
        GatherPDB()

    def ncbi_gene(self):
        self.pdb()
        GatherNCBIGene(use_existing_gene2pubmed=True)

    def bindingdb(self):
        self.ncbi_gene()
        import platform
        if platform.system() == 'Darwin':
            grep_path = '/usr/local/bin/grep'
            if not os.path.exists(grep_path):
                raise Exception('Please use Homebrew version of grep')
        else:
            grep_path = False
        GatherBindingDB(use_existing_bindingdb_data=True, grep_path=grep_path)

    def cbioportal(self):
        self.bindingdb()
        GatherCbioportalData()

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
