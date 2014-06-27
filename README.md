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
_XML_ branch.

The database is generated using a series of scripts which gather in data from
various public web resources. The first script to run is DoraCommit.py, which
initializes the necessary files and directory structure for a new database.
This should be followed by DoraGatherUniProt.py, which retrieves a set of
UniProt entries defined by a given search term. Subsequent scripts (once they
have been migrated over from the XML branch) add in data from various other
databases such as the PDB, NCBI Gene, cBioPortal, and BindingDB.

The above-mentioned scripts add data to the staging database. To commit this to
the public-facing, master database, run DoraCommit.py. A Flask HTTP server is
provided for the master database. It serves requests in JSON format.

A [frontend web client](https://github.com/choderalab/kinomeDB-webclient) is
currently in development, and an early working version can be seen here:

http://plfah2.mskcc.org/kinomeDB/

Manifest
--------

scripts/ - scripts for generating and analyzing the database

TargetExplorer/ - main code library

app\_master/ - data model and Flask HTTP server (master)

app\_stage/ - data model (staging)

resources/ - other miscellaneous files

Dependencies
------------

* BioPython, v1.62 or higher
* Various other Python packages commonly used in scientific computing. Recommended aproach is to install either Continuum Anaconda (https://store.continuum.io/) or Enthought Canopy (https://www.enthought.com/products/canopy/)
