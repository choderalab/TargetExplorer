#!/bin/bash
echo $TRAVIS_PULL_REQUEST $TRAVIS_BRANCH

deploy="true"

if [[ "$TRAVIS_PULL_REQUEST" != "false" ]]; then
    echo "This is a pull request. No deployment will be done."
    deploy="false"
fi

if [[ "$TRAVIS_BRANCH" != "master" ]]; then
    echo "No deployment on BRANCH='$TRAVIS_BRANCH'"
    deploy="false"
fi

echo $deploy

if [[ "$deploy" == "true" ]]; then
    binstar -t $BINSTAR_TOKEN upload --force -u choderalab -p targetexplorer-dev $HOME/miniconda/conda-bld/linux-64/targetexplorer-dev-*
fi
