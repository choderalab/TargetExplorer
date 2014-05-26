import os, sys

this_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, this_dir)
from config import targetexplorer_install_dir
sys.path.insert(0, targetexplorer_install_dir)

import app_master

if __name__ == '__main__':
    app_master.app.debug = True
    app_master.app.run()
