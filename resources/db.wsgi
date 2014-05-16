import os, sys

this_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, this_dir)
from config import targetexplorer_install_dir
sys.path.insert(0, targetexplorer_install_dir)

from flaskapp import app as app

if __name__ == '__main__':
    app.debug = True
    app.run()
