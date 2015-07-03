import argparse
from targetexplorer.PDB import GatherPDB


argparser = argparse.ArgumentParser(description='Gather cBioPortal data')
args = argparser.parse_args()
GatherPDB()
