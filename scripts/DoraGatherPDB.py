import argparse
from targetexplorer.protein_databank import GatherPDB


argparser = argparse.ArgumentParser(description='Gather cBioPortal data')
args = argparser.parse_args()
GatherPDB()
