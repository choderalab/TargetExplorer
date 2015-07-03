#!/usr/bin/env python
import os
import argparse
from targetexplorer.core import project_config_filename
from targetexplorer.initproject import InitProject


def parse_arguments():
        argparser = argparse.ArgumentParser(description='Initialize TargetExplorer database')
        argparser.add_argument(
            '--db_name', type=str, required=False, help='Database name, without extension'
        )
        return argparser.parse_args()

args = parse_arguments()

InitProject(
    db_name=args.db_name,
    project_path=os.path.abspath('.')
)
print(
    'Please now edit the UniProt search options in {0} before running the '
    'database generation scripts.'.format(project_config_filename)
)
