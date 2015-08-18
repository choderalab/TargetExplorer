#!/bin/bash
# This installs the program and runs unit tests
set -e
conda build devtools/conda-recipe
conda install --yes --use-local targetexplorer-dev
conda install --yes nose
pushd .; cd /
nosetests targetexplorer -v --exe -a unit
popd