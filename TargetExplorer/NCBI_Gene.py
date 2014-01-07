import urllib2, re, sys, os, gzip, datetime
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
    # Update metadata
    metadata_root = etree.parse(TargetExplorer.DB.external_data_metadata_filepath, parser).getroot()
    gene2pubmed_node = metadata_root.find('NCBI_Gene/gene2pubmed')
    gene2pubmed_node.set('filepath', gene2pubmed_gzfilepath)
    now = datetime.datetime.utcnow()
    datestamp = now.strftime(TargetExplorer.DB.datestamp_format_string)
    gene2pubmed_node.set('datestamp', datestamp)
    with open(TargetExplorer.DB.external_data_metadata_filepath, 'w') as xml_out_file:
        xml_out_file.write(etree.tostring(metadata_root, pretty_print=True))

def extract_PMID_from_Gene2PubMed_line(line, query_GeneIDs):
    words = line.split('\t', 3)
    #tax_id = words[0]
    #GeneID = words[1]
    #PMID = words[2]
    if words[1] in query_GeneIDs:
        return (words[1], words[2].strip())

def get_publications(gene2pubmed_filepath, query_GeneIDs):
    '''
    gene2pubmed_filepath: path to compressed gene2pubmed txt file
    query_GeneIDs: list of GeneID strings
    returns XML object with structure: gene[@GeneID=""]/publications/publication[@PMID=""]
    '''
    from multiprocessing import Pool
    from functools import partial

    with gzip.open(gene2pubmed_filepath, 'r') as gene2pubmed_file:
        gene2pubmed_lines = gene2pubmed_file.readlines()[1:] # ignore first line ("#format: ...")
    print 'Uncompressed Gene2PubMed file contains %d lines' % len(gene2pubmed_lines)

    #xml_root = etree.Element('root')
    partial_extract_PMID_from_Gene2PubMed_line = partial(extract_PMID_from_Gene2PubMed_line, query_GeneIDs=query_GeneIDs)

    pool = Pool()
    results = pool.map(partial_extract_PMID_from_Gene2PubMed_line, gene2pubmed_lines)
    # the mapped method will return None for each non-matching line, so strip these out
    results = [ r for r in results if r != None ]

    PMIDs_by_GeneID = {GeneID : [] for GeneID in query_GeneIDs}
    for result in results:
        PMIDs_by_GeneID[result[0]].append(result[1])

    return PMIDs_by_GeneID

    #for result in results:
    #    GeneID = result[0]
    #    PMID = result[1]

    #for GeneID in query_GeneIDs:
    #    gene_node = etree.SubElement(xml_root, 'gene')
    #    gene_node.set('GeneID', GeneID)
    #    etree.SubElement(gene_node, 'publications')
    #    for line in gene2pubmed_file:
    #        if line[0:8] == '#Format:':
    #            continue
    #        words = line.split('\t')
    #        tax_id = words[0]
    #        GeneID = words[1]
    #        if GeneID in query_GeneIDs:
    #            PMID = words[2]
    #            gene_pubs_node = xml_root.find('gene[@GeneID="%(GeneID)s"]/publications' % vars())
    #            pub_node = etree.SubElement(gene_pubs_node, 'publication')
    #            pub_node.set('PMID', PMID)
    #return xml_root

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


