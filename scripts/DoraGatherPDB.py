import argparse
from targetexplorer.pdb_tmp import GatherPDB


argparser = argparse.ArgumentParser(description='Gather cBioPortal data')
args = argparser.parse_args()
GatherPDB()
