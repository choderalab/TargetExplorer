import os, urllib, urllib2
from lxml import etree
from StringIO import StringIO
import Bio.Seq, Bio.Alphabet
import TargetExplorer as clab


webpage_index = 'http://plasmid.med.harvard.edu/PLASMID/GetCollection.do?collectionName=HIP%20FLEXGene%20human%20kinase%20collection%20(pDNR-Dual;%20C-terminal%20tag%20only)'

with open('cloneIDs.txt', 'r') as cloneIDs_file:
    # cloneIDs copied manually from the index page
    cloneIDs = [line.split()[0] for line in cloneIDs_file.readlines()]

cloneIDs_ints = [int(cloneID.replace('HsCD', '')) for cloneID in cloneIDs]

baseurl = 'http://plasmid.med.harvard.edu/PLASMID/GetCloneDetail.do?cloneid='

maxreadlength=10000000 # 1 MB

# read in DB

parser = etree.XMLParser(remove_blank_text=True)
DB_filepath = os.path.join('..', '..', '..', '..', 'database', 'database.xml')
DB_root = etree.parse(DB_filepath, parser).getroot()

# iterate through clones from plasmid library

with open('plasmid_insert_data.csv', 'w') as ofile:

    output_header_string = ','.join(['cloneID_int', 'cloneID', 'gene_symbol', 'UniProtAC', 'UniProt_entry_name', 'mutation?', 'discrepancy?', 'insert_dna_seq', 'insert_aa_seq', 'UniProt_canonseq']) + '\n'
    ofile.write(output_header_string)

    for cloneID_int in cloneIDs_ints:
        url = baseurl + str(cloneID_int)
        response = urllib2.urlopen(url)
        page = response.read(maxreadlength)

        parser = etree.HTMLParser()
        html_tree = etree.parse(StringIO(page), parser).getroot()
        insert_info_tablerow = html_tree.find('body/table/tr[@class="tableinfo"]')
        insert_info_data = [data.text for data in insert_info_tablerow.findall('td')]
        mutation = insert_info_data[3] # 'Yes' or 'No'
        discrepancy = insert_info_data[4] # 'Yes' or 'No'
        gene_symbol = insert_info_data[8] # ?HGNC

        cloneID = html_tree.find('body/table[1]/tr[1]/td[2]').text # e.g. 'HsCD00001029'

        insert_dna_seq = ''.join( html_tree.find('body/pre').text.splitlines() )

        # convert dna to aa
        insert_aa_seq = Bio.Seq.Seq(insert_dna_seq, Bio.Alphabet.generic_dna).translate(to_stop=True)


        # match to DB entry
        # first try HGNC gene symbol
        DB_HGNC = DB_root.find('entry/HGNC/entry[@Approved_Symbol="%s"]' % gene_symbol)
        if DB_HGNC != None:
            UniProt_node = DB_HGNC.getparent().getparent().find('UniProt')
            UniProtAC = UniProt_node.get('AC')
            print 'Gene symbol %s matched to UniProt AC %s (via HGNC gene symbol)' % (gene_symbol, UniProtAC)
        else:
            # if not found, try the UniProt gene symbols
            DB_UniProt_gene_name = DB_root.xpath('entry/UniProt/gene_names/gene_name[text()="%s"]' % gene_symbol)
            if DB_UniProt_gene_name == []:
                print 'Gene symbol %s could not be matched to the DB\n' % gene_symbol
                continue
            UniProt_node = DB_UniProt_gene_name[0].getparent().getparent()
            UniProtAC = UniProt_node.get('AC')
            print 'Gene symbol %s matched to UniProt AC %s (via UniProt gene name)' % (gene_symbol, UniProtAC)

        UniProt_entry_name = UniProt_node.get('entry_name')
        UniProt_canonseq = clab.core.sequnwrap( UniProt_node.find('isoforms/canonical_isoform/sequence').text )

        print ''

        output_string = ','.join( [str(cloneID_int), cloneID, gene_symbol, UniProtAC, UniProt_entry_name, mutation, discrepancy, insert_dna_seq, str(insert_aa_seq), UniProt_canonseq] ) + '\n'
        ofile.write(output_string)

