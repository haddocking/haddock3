"""I/O helper."""

def open_files_to_lines(*files):
    """
    Open files to lines.

    New-lines are stripped.

    Returns
    -------
    list of lists of strings
        The lines of the files.
        Input order is maintained.
    """
    f_paths = map(Path, files)
    return [f.read_text().split(os.linesep) for f in f_paths]
