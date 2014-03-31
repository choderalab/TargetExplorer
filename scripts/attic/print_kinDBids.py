# Simple script to print kinDB_ids for kinases returned by a given xpath query
# xpath query should point to kinome/kinase/uniprot/pk_domain elements
#
# Daniel L. Parton <partond@mskcc.org> - 3 June 2013
#

import sys
from lxml import etree

xpath_query = sys.argv[1]
print xpath_query

kinDB = etree.parse('kinDB-complete.xml').getroot()

report = kinDB.xpath(xpath_query)

for r in report:
    print r.get('kinDB_id')

