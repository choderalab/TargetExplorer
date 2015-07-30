import os
from lxml import etree
import gzip
from targetexplorer.flaskapp import models
from targetexplorer.uniprot import retrieve_uniprot, GatherUniProt
from targetexplorer.core import xml_parser
from targetexplorer.utils import get_installed_resource_filepath
from targetexplorer.tests.utils import projecttest_context
from nose.plugins.attrib import attr


@attr('network')
def test_retrieve_uniprot():
    xml_text = retrieve_uniprot('mnemonic:ABL1_HUMAN')
    xml_root = etree.fromstring(xml_text, xml_parser)


@attr('network')
def test_retrieve_uniprot_matches_reference():
    xml_text = retrieve_uniprot('mnemonic:ABL1_HUMAN')
    xml_root = etree.fromstring(xml_text, xml_parser)
    ref_xml_filepath = get_installed_resource_filepath(os.path.join('resources', 'uniprot-search-abl1.xml.gz'))
    with gzip.open(ref_xml_filepath) as ref_xml_file:
        ref_xml_root = etree.parse(ref_xml_file, xml_parser)

    entry_name = xml_root.find('entry/name').text
    ref_entry_name = ref_xml_root.find('entry/name').text
    assert entry_name == ref_entry_name


@attr('unit')
def test_gather_uniprot():
    with projecttest_context(set_up_project_stage='init'):
        GatherUniProt(
            uniprot_query='mnemonic:ABL1_HUMAN',
            uniprot_domain_regex='^Protein kinase(?!; truncated)(?!; inactive)',
            use_existing_data=True
        )
        first_uniprot_entry = models.UniProt.query.first()
        first_uniprot_domain = models.UniProtDomain.query.first()
        first_target_domain = models.UniProtDomain.query.filter_by(is_target_domain=True).first()
        first_pdb_chain = models.PDBChain.query.first()
        assert first_uniprot_entry.entry_name == 'ABL1_HUMAN'
        assert first_uniprot_domain.domain_id == 0
        assert first_target_domain.target_id == 'ABL1_HUMAN_D0'
        assert first_target_domain.domain_id == 2
        assert first_pdb_chain.uniprotdomain.domain_id == 1


@attr('unit')
def test_gather_uniprot_no_domain_regex():
    with projecttest_context(set_up_project_stage='init'):
        GatherUniProt(
            uniprot_query='mnemonic:ABL1_HUMAN',
            use_existing_data=True
        )
        first_uniprot_entry = models.UniProt.query.first()
        first_uniprot_domain = models.UniProtDomain.query.first()
        assert first_uniprot_entry.entry_name == 'ABL1_HUMAN'
        assert first_uniprot_domain.domain_id == 0
        assert len(models.UniProtDomain.query.filter_by(is_target_domain=True).all()) == 0


@attr('network')
def test_gather_uniprot_using_network():
    with projecttest_context(set_up_project_stage='init'):
        GatherUniProt(
            uniprot_query='mnemonic:ABL1_HUMAN',
            uniprot_domain_regex='^Protein kinase(?!; truncated)(?!; inactive)',
        )
        first_uniprot_entry = models.UniProt.query.first()
        first_uniprot_domain = models.UniProtDomain.query.first()
        first_target_domain = models.UniProtDomain.query.filter_by(is_target_domain=True).first()
        assert first_uniprot_entry.entry_name == 'ABL1_HUMAN'
        assert first_uniprot_domain.domain_id == 0
        assert first_target_domain.target_id == 'ABL1_HUMAN_D0'
        assert first_target_domain.domain_id == 2
