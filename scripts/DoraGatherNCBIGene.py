import argparse
from targetexplorer.ncbi_gene_tmp import GatherNCBIGene

argparser = argparse.ArgumentParser(description='Gather NCBI Gene data')
argparser.add_argument(
    '--use_existing_gene2pubmed',
    help='Do not download a new gene2pubmed.gz file. Only works if an existing file is present.',
    action='store_true',
    default=False
)
args = argparser.parse_args()

GatherNCBIGene(args.use_existing_gene2pubmed)
