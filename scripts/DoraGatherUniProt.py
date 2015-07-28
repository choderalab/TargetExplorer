#!/usr/bin/env python
#==============================================================================
# Imports
#==============================================================================

import argparse
from targetexplorer.core import read_project_config
from targetexplorer.uniprot import GatherUniProt

def parse_arguments():
        argparser = argparse.ArgumentParser(description='Gather UniProt')
        argparser.add_argument('--use_existing_data', help='Do not download a new UniProt document. Only works if an existing document is present.', action='store_true', default=False)
        argparser.add_argument('--count_nonselected_domain_names', help='Count the number of occurrences of domain names which are not selected by the regex and write to the "selected_domain_names.txt" file (may take a little while)', action='store_true', default=False)
        return argparser.parse_args()

args = parse_arguments()

project_config = read_project_config()

GatherUniProt(
    use_existing_data=args.use_existing_data,
    count_nonselected_domain_names=args.count_nonselected_domain_names,
    uniprot_query=project_config['uniprot_query'],
    uniprot_domain_regex=project_config['uniprot_domain_regex'],
)

print('Done.')
