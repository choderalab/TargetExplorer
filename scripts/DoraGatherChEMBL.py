__author__ = 'isikm'

import argparse
from targetexplorer.chembl import GatherChemblTarget

argparser = argparse.ArgumentParser(description='Gather ChEMBL Target data')
args = argparser.parse_args()

GatherChemblTarget()