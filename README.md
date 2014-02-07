TargetExplorerDB
================

Suggested procedure for integrating these web pages and associated scripts with the main TargetExplorerDB branch
----------------------------------------------------------------------------------------------------------------

Make a separate clone of the TargetExplorerDB repo.

cd into the cloned repo.

    git checkout --orphan gh-pages

Go back to your original development clone of TargetExplorerDB.

Copy the new repo directory to the original repo directory, renaming it to 'gh-pages'.

'gh-pages' should be present in .gitignore in the top directory of the original repo.

The gh-pages branch should access files from the main branch by using the appropriate relative paths.

