TargetExplorer
==============

Authors
-------

* Daniel L. Parton | daniel.parton@choderalab.org
* John D. Chodera | john.chodera@choderalab.org

Overview
--------

_Note:_ The code is in the process of being migrated from an XML-based database
to an SQLite implementation. The old, stable XML-based code can be found in the
branch _XML_.

The database is generated using a series of scripts which gather in data from
various public web resources. The first script to run is tedb\_init.py, which
initializes a new database. This should be followed by
tedb\_gather\_uniprot.py, which retrieves a set of UniProt entries defined by a
given search term. These UniProt entries define the top-level nodes of this
database. Subsequent scripts add in data from various other databases such as
the PDB, NCBI Gene, cBioPortal, and BindingDB.

A Flask HTTP server is also provided, which serves database requests in JSON
format. An example frontend client can be seen here:

http://ec2-54-227-62-182.compute-1.amazonaws.com/kinomeDB/

Manifest
--------

scripts/ - scripts for generating and analyzing the database

TargetExplorer/ - main code library

flaskapp/ - Flask HTTP server

resources/ - other miscellaneous files

Dependencies
------------

* BioPython, v1.62 or higher
* Various other Python packages commonly used in scientific computing. Recommended aproach is to install either Continuum Anaconda (https://store.continuum.io/) or Enthought Canopy (https://www.enthought.com/products/canopy/)
