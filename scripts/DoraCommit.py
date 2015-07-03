import argparse
from targetexplorer.commit import Commit

argparser = argparse.ArgumentParser(description='Commit database')
args = argparser.parse_args()

Commit()
