#!/usr/bin/env python
#
# Initialize a SQLite database with the TargetExplorer schema
#
# Daniel L. Parton <daniel.parton@choderalab.org> - 9 Apr 2014

import os, argparse, yaml
import sqlite3

# ========
# Parse command-line arguments
# ========

# argparser = argparse.ArgumentParser(description='Initialize a MySQL database with the TargetExplorer schema.')
# argparser.add_argument('--db_name', type=str, help='Database name (e.g. "kinome")')
# args = argparser.parse_args()

if not os.path.exists('config.yaml'):
    raise Exception, 'config.yaml not found'
with open('config.yaml') as config_file:
    config = yaml.load(config_file)
    print config

if not os.path.exists('database'):
    os.mkdir('database')

if not os.path.exists('external-data'):
    os.mkdir('external-data')

db_path = os.path.join('database', config.get('db_name') + '.db')

with sqlite3.connect(db_path) as connection:
    cursor = connection.cursor()
    connection.commit()

