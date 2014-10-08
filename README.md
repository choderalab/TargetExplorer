TargetExplorer
==============

Authors
-------

* Daniel L. Parton | daniel.parton@choderalab.org
* John D. Chodera | john.chodera@choderalab.org

Overview
--------

Database framework for storing genomic, structural and functional data for a
given protein family, with RESTful API.

The code is built on [Flask](http://flask.pocoo.org/) (a Python web framework)
and [SQLAlchemy](http://www.sqlalchemy.org/) (an object-relational mapper which
maps between SQL databases and Python objects).

The database is generated using a series of scripts which gather in data from
various public web resources. The first script to run is DoraInit.py, which
initializes the necessary files and directory structure for a new database.
This should be followed by DoraGatherUniProt.py, which retrieves a set of
[UniProt](http://www.uniprot.org/) entries defined by a given search term.
Subsequent scripts (once they have been migrated over from the old _XML_
branch) add in data from various other databases such as the
[PDB](http://www.rcsb.org), [NCBI Gene](http://www.ncbi.nlm.nih.gov/gene),
[cBioPortal](http://www.cbioportal.org), and
[BindingDB](http://www.bindingdb.org/bind/index.jsp). Finally, DoraCommit.py is
used to carry out a sanity check on the new data; if this passes, the new data
will be exposed through the API.

A [frontend web client](https://github.com/choderalab/kinomeDB-webclient) is
currently in development, and an early working version can be seen here:

http://plfah2.mskcc.org/kinomeDB/

Manifest
--------

scripts/ - scripts for generating and analyzing the database

TargetExplorer/ - main code library

flask\_app/ - data model and Flask HTTP server (master)

resources/ - other miscellaneous files

Dependencies
------------

* BioPython, v1.62 or higher
* Various other Python packages commonly used in scientific computing. Recommended aproach is to install either Continuum Anaconda (https://store.continuum.io/) or Enthought Canopy (https://www.enthought.com/products/canopy/)

Notes on database generation process
------------------------------------

A "crawl number" is iteratively assigned for each pass through the database
generation process, from DoraGatherUniProt.py to DoraCommit.py. If the process
is completed successfully, DoraCommit.py will update the "safe crawl number",
which tells the API to work with only the data corresponding to that crawl
number. The number of crawls to store in the database can be defined by the
user, and is set by default to 5. Older crawls are deleted by DoraCommit.py.
