TargetExplorer
==============

[![Build Status](https://travis-ci.org/choderalab/TargetExplorer.svg)](https://travis-ci.org/choderalab/TargetExplorer)

Database framework with RESTful API for aggregating genomic, structural, and functional data for target protein families.

Authors
-------

* Daniel L. Parton | daniel.parton@choderalab.org
* Mehtap Isik | mehtap.isik@choderalab.org
* John D. Chodera | john.chodera@choderalab.org

Overview
--------

The code is built using the [Flask](http://flask.pocoo.org/) a Python web framework
and [SQLAlchemy](http://www.sqlalchemy.org/) - an object-relational mapper which
maps between Python objects and SQL databases.

The database is generated using a series of scripts which gather in data from
various public web resources. The first script to run is DoraInit.py, which
initializes the necessary files and directory structure for a new database.
This should be followed by DoraGatherUniProt.py, which retrieves a set of
[UniProt](http://www.uniprot.org/) entries defined by a given search term.
Subsequent scripts add in data from various other databases such as the
[PDB](http://www.rcsb.org), [NCBI Gene](http://www.ncbi.nlm.nih.gov/gene),
[cBioPortal](http://www.cbioportal.org), and
[BindingDB](http://www.bindingdb.org/bind/index.jsp). Finally, DoraCommit.py is
used to ensure that all gather scripts have been run since the last commit;
if this passes, the new data is *committed*, meaning it can then be exposed through the API.

A [frontend web client](https://github.com/choderalab/kinomeDB-webclient) is
currently in development.

Installation using Anaconda (recommended approach)
--------------------------------------------------

First install [Anaconda](https://store.continuum.io/cshop/anaconda/) - a free and awesome Python
distribution for scientific and data-intensive applications.

```.sh
conda config --add channels http://anaconda.org/choderalab
conda install targetexplorer
```

Notes on database generation process
------------------------------------

A "crawl number" is iteratively assigned for each pass through the database
generation process, from DoraGatherUniProt.py to DoraCommit.py. If the process
is completed successfully, DoraCommit.py will update the "safe crawl number",
which tells the API to work with only the data corresponding to that crawl
number. The number of crawls to store in the database can be defined by the
user, and is set by default to 5. Older crawls are deleted by DoraCommit.py.
