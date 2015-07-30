import argparse
from targetexplorer.pdb import GatherPDB


argparser = argparse.ArgumentParser(description='Gather PDB data')
args = argparser.parse_args()
GatherPDB()
