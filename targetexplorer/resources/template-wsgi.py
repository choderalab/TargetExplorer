# import os, sys

# this_dir = os.path.dirname(os.path.realpath(__file__))
# sys.path.insert(0, this_dir)
# from project_config import targetexplorer_install_dir
# sys.path.insert(0, targetexplorer_install_dir)

from targetexplorer.flaskapp import app

if __name__ == '__main__':
    from targetexplorer.flaskapp import app, views
    app.debug = True
    app.run()
