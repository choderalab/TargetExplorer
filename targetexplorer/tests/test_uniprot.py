import os

from lxml import etree

from targetexplorer.flaskapp import models
from targetexplorer.uniprot import retrieve_uniprot, GatherUniProt
from targetexplorer.core import xml_parser
from targetexplorer.utils import get_installed_resource_filepath
from targetexplorer.tests.utils import projecttest_context


def test_retrieve_uniprot():
    xml_text = retrieve_uniprot('mnemonic:ABL1_HUMAN')
    xml_root = etree.fromstring(xml_text, xml_parser)


def test_retrieve_uniprot_matches_reference():
    xml_text = retrieve_uniprot('mnemonic:ABL1_HUMAN')
    xml_root = etree.fromstring(xml_text, xml_parser)
    ref_xml_filepath = get_installed_resource_filepath(os.path.join('resources', 'uniprot-search.xml'))
    with open(ref_xml_filepath) as ref_xml_file:
        ref_xml_root = etree.parse(ref_xml_file, xml_parser)

    entry_name = xml_root.find('entry/name').text
    ref_entry_name = ref_xml_root.find('entry/name').text
    assert entry_name == ref_entry_name


def test_gather_uniprot():
    with projecttest_context(set_up_project_stage='init'):
        GatherUniProt(
            uniprot_query='mnemonic:ABL1_HUMAN',
            uniprot_domain_regex='^Protein kinase(?!; truncated)(?!; inactive)',
        )
        first_uniprot_entry = models.UniProt.query.first()
        first_uniprot_domain = models.UniProtDomain.query.first()
        assert first_uniprot_entry.entry_name == 'ABL1_HUMAN'
        assert first_uniprot_domain.targetid == 'ABL1_HUMAN_D0'
