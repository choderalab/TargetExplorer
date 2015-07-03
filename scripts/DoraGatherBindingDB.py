import argparse
from targetexplorer.BindingDB import GatherBindingDB

argparser = argparse.ArgumentParser(description='Gather BindingDB data')
argparser.add_argument(
    '--use_existing_bindingdb_data',
    help='Do not download a new BindingDB_All.tab file. Only works if an existing file is present.',
    action='store_true',
    default=False
)
argparser.add_argument(
    '--grep_path',
    help='Provide explicit path for grep. Otherwise, searches in the usual places.',
    type=str,
    default='grep'
)
args = argparser.parse_args()

GatherBindingDB(
    use_existing_bindingdb_data=args.use_existing_bindingdb_data,
    grep_path=args.grep_path
)
