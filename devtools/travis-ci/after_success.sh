#!/bin/bash

echo $TRAVIS_PULL_REQUEST
echo $TRAVIS_BRANCH

if [[ "$TRAVIS_PULL_REQUEST" == "true" ]]; then
    echo "This is a pull request. No deployment will be done."; exit 0
fi


if [[ "$TRAVIS_BRANCH" != "master" ]]; then
    echo "No deployment on BRANCH='$TRAVIS_BRANCH'"; exit 0
fi

# Deploy to binstar.
conda install --yes anaconda-client jinja2
binstar -t $BINSTAR_TOKEN upload --force -u choderalab -p targetexplorer-dev $HOME/miniconda/conda-bld/*/targetexplorer-dev-*.tar.bz2
