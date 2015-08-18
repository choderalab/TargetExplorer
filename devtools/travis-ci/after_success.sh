#!/bin/bash
echo $TRAVIS_PULL_REQUEST $TRAVIS_BRANCH

deploy="true"

if [[ "$TRAVIS_PULL_REQUEST" == "true" ]]; then
    echo "This is a pull request. No deployment will be done."
    deploy="false"
fi
if [[ "$TRAVIS_BRANCH" != "master" ]]; then
    echo "No deployment on BRANCH='$TRAVIS_BRANCH'"
    deploy="false"
fi

if [[ "$deploy" == "true" ]]; then
    binstar -t $BINSTAR_TOKEN upload --force -u choderalab -p targetexplorer-dev $HOME/miniconda/conda-bld/linux-64/targetexplorer-dev-*
fi

# Run tests with external dependencies (Travis build should not depend upon this since, for example, external sites could be temporarily down)
nosetests targetexplorer -v --exe -A 'network and not slow'
