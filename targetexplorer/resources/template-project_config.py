import os
# UniProt search options
uniprot_query_string = 'EXAMPLE... domain:"protein kinase" AND reviewed:yes'
uniprot_domain_regex = 'EXAMPLE... ^Protein kinase(?!; truncated)(?!; inactive)'

# General database options
ncrawls_to_save = 5

# Don't edit the code below here
db_name = 'DB_NAME'
project_basedir = os.path.abspath(os.path.dirname(__file__))
targetexplorer_install_dir = 'TARGETEXPLORER_INSTALL_DIR'
