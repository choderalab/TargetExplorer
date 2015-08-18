import argparse
from targetexplorer.cbioportal import GatherCbioportalData, external_cbioportal_data_filepath


argparser = argparse.ArgumentParser(description='Gather cBioPortal data')
argparser.add_argument(
    '--use_existing_data',
    help='Only works if an existing file {0} is present.'.format(external_cbioportal_data_filepath),
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

GatherCbioportalData(
    use_existing_cbioportal_data=args.use_existing_data,
    use_existing_oncotator_data=args.use_existing_data,
    commit_to_db=args.nocommit
)
