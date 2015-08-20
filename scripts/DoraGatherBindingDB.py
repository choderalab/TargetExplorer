import argparse
from targetexplorer.bindingdb import GatherBindingDB

argparser = argparse.ArgumentParser(description='Gather BindingDB data')
argparser.add_argument(
    '--use_existing_bindingdb_data',
    help='Do not download a new BindingDB_All.tab file. Only works if an existing file is present.',
    action='store_true',
    default=False
)
argparser.add_argument(
    '--nocommit',
    help='Run script, but do not commit to database.',
    action='store_true',
    default=False
)
args = argparser.parse_args()

GatherBindingDB(
    use_existing_bindingdb_data=args.use_existing_bindingdb_data,
    commit_to_db=not args.nocommit
)
