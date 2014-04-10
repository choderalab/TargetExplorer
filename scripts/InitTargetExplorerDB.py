#!/usr/bin/env python
#
# Initialize a MySQL database with the TargetExplorer schema
#
# Daniel L. Parton <daniel.parton@choderalab.org> - 9 Apr 2014

import os, argparse, getpass
import TargetExplorer
import MySQLdb as mdb

# ========
# Parse command-line arguments
# ========

argparser = argparse.ArgumentParser(description='Initialize a MySQL database with the TargetExplorer schema.')
argparser.add_argument('--db_name', type=str, help='Database name (e.g. kinome)')
# argparser.add_argument('--host', type=str, help='MySQL hostname')
argparser.add_argument('--user', type=str, help='MySQL username')
args = argparser.parse_args()

user_password = getpass.getpass('Please enter your MySQL password:')

connection = mdb.connect(user=args.user, passwd=user_password)
with connection:
    cursor = connection.cursor()
    cursor.execute("CREATE DATABASE IF NOT EXISTS %s" % args.db_name)

connection = mdb.connect(host='localhost', user=args.user, passwd=user_password, db=args.db_name)

with connection:
    cursor = connection.cursor(mdb.cursors.DictCursor)
    cursor.execute("CREATE TABLE IF NOT EXISTS db_entries(_id INT PRIMARY KEY AUTO_INCREMENT)")

