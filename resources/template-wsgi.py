import os, sys

this_dir = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0, this_dir)
from config import targetexplorer_install_dir
sys.path.insert(0, targetexplorer_install_dir)

from app_master import app

if __name__ == '__main__':
    import app_master
    app_master.app.debug = True
    app_master.app.run()
