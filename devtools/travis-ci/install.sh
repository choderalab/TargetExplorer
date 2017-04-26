#!/bin/bash
MINICONDA=Miniconda-latest-Linux-x86_64.sh
MINICONDA_HOME=$HOME/miniconda
MINICONDA_MD5=$(curl -s http://repo.continuum.io/miniconda/ | grep -A3 $MINICONDA | sed -n '4p' | sed -n 's/ *<td>\(.*\)<\/td> */\1/p')
wget http://repo.continuum.io/miniconda/$MINICONDA
if [[ $MINICONDA_MD5 != $(md5sum $MINICONDA | cut -d ' ' -f 1) ]]; then
    echo "Miniconda MD5 mismatch"
    exit 1
fi
bash $MINICONDA -b -p $MINICONDA_HOME

export PATH=$MINICONDA_HOME/bin:$PATH

sudo apt-get update

conda update --yes conda
source activate $python
conda install --yes conda-build anaconda-client jinja2
conda config --add channels http://conda.anaconda.org/choderalab
