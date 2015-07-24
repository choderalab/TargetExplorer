# Testing framework

Testing a Flask-SQLAlchemy database can be tricky.
The following is a general overview of the custom testing framework used in TargetExplorer.

## General description

"Unit tests" should not have any external database dependencies (e.g. UniProt).
They should use reference files bundled with TargeteExplorer as a source of data.

Tests which do require network access to external database dependencies should be given the
attribute "network".

Tests of the database are run using a sample TargetExplorer project which includes only the human
protein Abl1.

* Custom nose plugin (`targetexplorer.tests.noseplugins.setup_tmp_db_plugin`):
    * This allows code to be executed by nosetests before and after the actual tests are run. 
    * At beginning:
        * Creates a temporary directory
        * Writes the temporary directory path to a file installed in the TargetExplorer resources
        directory: 'resources/testdir'
        * chdir to temp dir
        * Create a blank TargetExplorer database in temporary directory
    * All nosetests are then run in this temporary directory, and can access the database
    (using the path in 'resources/testdir')
    * At end:
        * Finally the plugin deletes the temporary directory
        * (the user remains in the original directory - the temporary directory is used only
        internally by nosetests)

* Custom context manager for tests which require a database
(`targetexplorer.tests.utils.projecttest_context`):
    * Uses reference data files to set up the database at various points of completion.
    * This allows tests to run without retrieving data from external databases over a network
    connection
    * Following command will set up the database at the point after retrieval of UniProt data:
    ```.py
    with projecttest_context(set_up_project_stage='uniprot'):
        ...
    ```
    * Other stages which can be set up are:
        * 'blank' - just a blank database file; no contained data
        * 'init' - as if user had run DoraInit.py
        * 'uniprot' - as if user had run DoraInit.py, then DoraGatherUniProt.py
        * 'pdb' - as if user had run above scripts, then DoraGatherPDB.py
        * 'ncbi_gene' - etc.
        * 'bindingdb' - etc.
        * 'cbioportal' - etc.
        * 'committed' - as if all gather scripts have been run, and DoraCommit.py
    * Order of operation:
        * First reads the temporary test directory from 'resources/testdir'
        * chdir to temp dir
        * Sets up the sample project at the given stage
        * `yield` (test is run at this point)
        * After the test completes (whether successful or not), the context manager finishes by
        dropping all data in the database, allowing the next test to start from a fresh blank database
