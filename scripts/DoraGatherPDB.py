import argparse
from targetexplorer.pdb import GatherPDB

argparser = argparse.ArgumentParser(description='Gather PDB data')
argparser.add_argument(
    '--nocommit',
    help='Run script, but do not commit to database.',
    action='store_true',
    default=False
)

args = argparser.parse_args()
GatherPDB(commit_to_db=not args.nocommit)
