# Retrieve HGNC flat txt file and parse gene nomenclature info
#
# Daniel L. Parton <partond@mskcc.org> - 11 Oct 2013
#

import os
import choderalab as clab
from lxml import etree

external_data_dir = os.path.join('external-data', 'HGNC')
HGNC_dir = os.path.join('external-data', 'HGNC')
if not os.path.exists(HGNC_dir):
    os.mkdir(HGNC_dir)

HGNC_data_filepath = os.path.join(HGNC_dir, 'HGNC-data.txt.gz')

clab.HGNC.retrieve_data_file(HGNC_data_filepath)

HGNC_data = clab.HGNC.HGNC_data(HGNC_data_filepath)

print HGNC_data.IDkeyed_data['HGNC:9376']

