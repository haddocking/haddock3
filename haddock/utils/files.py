from os.path import dirname, join


MAIN_DIRECTORY = dirname(join(dirname(__file__), '..', '..', '..'))


def get_full_path(*path):
    return join(MAIN_DIRECTORY, *path)
