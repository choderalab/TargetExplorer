import urllib2
import datetime
import TargetExplorer
from lxml import etree

parser = etree.XMLParser(remove_blank_text=True)


def retrieve_gene2pubmed(gene2pubmed_gzfilepath):
    '''
    Download compressed gene2pubmed from NCBI FTP site and write to file.
    '''
    url = 'ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2pubmed.gz'
    response = urllib2.urlopen(url)
    page = response.read(100000000)
    with open(gene2pubmed_gzfilepath, 'w') as gene2pubmed_gzfile:
        gene2pubmed_gzfile.write(page)
    # # Update metadata
    # metadata_root = etree.parse(TargetExplorer.DB.external_data_metadata_filepath, parser).getroot()
    # gene2pubmed_node = metadata_root.find('NCBI_Gene/gene2pubmed')
    # if gene2pubmed_node is None:
    #     NCBI_Gene_node = etree.SubElement(metadata_root, 'NCBI_Gene')
    #     gene2pubmed_node = etree.SubElement(NCBI_Gene_node, 'gene2pubmed')
    # gene2pubmed_node.set('filepath', gene2pubmed_gzfilepath)
    # now = datetime.datetime.utcnow()
    # datestamp = now.strftime(TargetExplorer.DB.datestamp_format_string)
    # gene2pubmed_node.set('datestamp', datestamp)
    # with open(TargetExplorer.DB.external_data_metadata_filepath, 'w') as xml_out_file:
    #     xml_out_file.write(etree.tostring(metadata_root, pretty_print=True))


def retrieve_xml(GeneID):
    '''
    Don't use this to retrieve publication data.
    '''
    url = 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id=%(GeneID)s&retmode=xml' % vars()
    response = urllib2.urlopen(url)
    page = response.read(100000000000)
    lines = page.splitlines()

    ofilepath = 'tmp-gene.xml'
    with open(ofilepath, 'w') as ofile:
        ofile.write(page)
