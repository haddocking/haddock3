"""General utilities."""
import logging
import shutil


logger = logging.getLogger(__name__)


def get_result_or_same_in_list(function, value):
    """
    Return the result if True or the value within a list.

    Applies `function` to `value` and returns its result if it evaluates
    to True. Otherwise, return the value within a list.

    `function` should receive a single argument, the `value`.
    """
    result = function(value)
    return result if result else [value]


def make_list_if_string(item):
    if isinstance(item, str):
        return [item]
    return item


def copy_files_to_dir(paths, directory):
    """
    Copy files to directory.

    Parameters
    ----------
    paths : iterable of paths
        Source files.

    directory : path
        Where to copy files to.
    """
    for path in paths:
        shutil.copy(path, directory)


def remove_folder(folder):
    """Removes a folder if it exists."""
    if folder.exists():
        logger.warning(f'{folder} exists and it will be REMOVED!')
        shutil.rmtree(folder)
