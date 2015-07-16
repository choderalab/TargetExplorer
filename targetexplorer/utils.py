import os
import json
import targetexplorer
from pkg_resources import resource_filename


def get_installed_resource_filepath(relative_path):
    """Get the full path to one of the resource files shipped with TargetExplorer.
    In the source distribution, these files are in ``targetexplorer/resources``,
    but on installation, they are moved to somewhere in the user's python
    site-packages directory.
    This function uses the pkg_resources package to find the file within the installation directory
    structure.

    Parameters
    ----------
    relative_path : str
        Path of the file to load (relative to the ``targetexplorer`` source folder).

    Examples
    --------
    get_installed_resource_filename('resources/template-wsgi.py')
    """
    fn = resource_filename(targetexplorer.__name__, relative_path)
    if not os.path.exists(fn):
        raise ValueError(
            "Sorry! {0} does not exist."
            "If you just added it, you'll have to re-install".format(relative_path)
        )
    return fn


installation_testdir_filepath = get_installed_resource_filepath(
    os.path.join('resources', 'testdir')
)


def json_dump_pretty(obj, filepath):
    json.dump(obj, filepath, indent=4, separators=(',', ': '))


class DatabaseException(Exception):
    pass
