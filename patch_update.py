#!/usr/bin/env python
# -*- encoding: utf-8 -*-

"""Update HADDOCK patch version using git logs information."""

import json
import os
import re
from importlib import reload
from subprocess import PIPE, Popen, check_output

import haddock


class ProcessFailure(Exception):
    """Exception class devoted to unfound merge PR in `git log`.

    Parameters
    ----------
    nb_last_merges : int
        Total number of merges that have been parsed.
    """

    def __init__(self, nb_last_merges: int):
        self.nb = nb_last_merges

    def __str__(self):
        return f"ProcessFailure: {self.__repr__()}"

    def __repr__(self):
        return (f'After {self.nb} log merges, still no Merged Pull Request'
                ' found using git log.')


def get_N_recent_merges(n: int = 10):
    """Obtain the `N` most recent git merges.

    Parameters
    ----------
    n : int
        Number of most recent merges to return.

    Return
    ------
    loaded_logs : list
        List of `N` most recent merges.
    """
    # initiate git cmd line
    git_cmd_ = [
        'git', 'log',  # to get the logs
        '--merges',  # to get only merges
        '-n', str(n),  # to get only last N merges
        '--pretty=format:{"commit": "%H", "subject": "%f"}',  # noqa : E501 # gather commit hash and subject as dict
        ]

    # run the command
    git_logs = Popen(git_cmd_, stdout=PIPE)

    # gather outputs
    logs = git_logs.stdout.read().decode('utf-8')

    # cast string to dict
    loaded_logs = [json.loads(logstr) for logstr in logs.split(os.linesep)]

    # return list of last N logs
    return loaded_logs


def find_closest_pr_merge(logs: list) -> dict:
    """Find most recent merged Pull Request.

    Parameters
    ----------
    logs : list
        List of git logs.
        NOTE: here we consider they are sorted from most recent to ancient.

    Return
    ------
    logdt : dict
        First git log data containing a merged pull request in the list.
    """
    search_flag = re.compile(r'Merge-pull-request-\d{1,6}')
    for logdt in logs:
        if search_flag.search(logdt["subject"]):
            return logdt


def find_last_pr_merge(init_n: int = 10) -> dict:
    """Find most recent merged Pull Request.

    Parameters
    ----------
    init_n : int
        Initiale number of merges logs to parse.

    Return
    ------
    closest_merge_log : dict
        Dictionary containing the most recent merged PR log info.
    """
    previous_logs = []
    while True:
        # gather last n merges
        logs = get_N_recent_merges(init_n)

        # make sure we are still getting new log entries
        if logs == previous_logs:
            raise ProcessFailure(init_n)

        # search last pull request merge
        closest_merge_log = find_closest_pr_merge(logs)

        # check if something was returned
        if closest_merge_log:
            # return clostest merge log
            return closest_merge_log

        # increase number of merges commit to parse
        init_n += 10
        # modify controle variable
        previous_logs = logs
        # retry


def update_haddock_version(commit_id: str, shorten: int = 6):
    """Update HADDOCK version with `commit_id` as new patch.

    Parameters
    ----------
    commit_id : str
        Version of the current patch to apply.
    shorten : int
        Number of character to shorten commit id.

    Return
    ------
    current_version : str
        The current HADDOCK version.
    """
    # generate shorten version of the commit
    short_commit = commit_id[:shorten]
    # generate number version of commit
    number_commit = commit_to_number(short_commit)
    # modify patch version with current commit
    modify_patch_version(number_commit)
    # get current (just modified) version
    current_version = get_current_version()
    # check it was well applied
    assert current_version.split('.')[-1] == number_commit
    assert get_current_patch_version() == number_commit
    # return current version
    return current_version


def commit_to_number(commit_id: str):
    """Convert a commit id to a number.
    
    Parameters
    ----------
    commit_id : str
        A string containing letters / numbers.

    Return
    number_seq : str
        A string containing only numbers.
    """
    # initiate holding list
    number_seq_ = []
    # loop over letters
    for letter in commit_id:
        try:
            # check if an integer
            number = int(letter)
        except ValueError:
            # convert an ascii to letter
            number = ord(letter)
        finally:
            # add string casted version of the number
            number_seq_.append(str(number))
    # join the sequence
    number_seq = ''.join(number_seq_)
    # return it
    return number_seq


def modify_patch_version(patch_v: str):
    """Optain current haddock3 patch version.

    Parameters
    ----------
    patch_v : str
        Version of the current patch to apply.
    """
    # path where it is stored
    fpath = haddock.__file__
    # build sed cmd line
    sed_cmd_ = f"""sed -i 's|v_patch =.*|v_patch = "{patch_v}"|g' {fpath}"""
    # use sed to modify patch version
    check_output(sed_cmd_, shell=True)


def load_current_haddock():
    """Reload HADDOCK."""
    reload(haddock)


def get_current_version():
    """Optain current haddock version."""
    load_current_haddock()
    return haddock.version


def get_current_patch_version():
    """Optain current haddock3 patch version."""
    load_current_haddock()
    return haddock.v_patch


def main() -> str:
    """Process the finding of last merge of pull request.

    Return
    ------
    version : str
        Current HADDOCK version
    """
    # get current patch version
    init_patch_version = get_current_patch_version()

    # search last pull request merge
    try:
        last_merge_log = find_last_pr_merge(10)
    except ProcessFailure:
        return get_current_version()
    else:
        # get last merge PR commit hash
        commit_id = last_merge_log["commit"]

    # update haddock3 version
    try:
        version = update_haddock_version(commit_id)
    except AssertionError:
        version = update_haddock_version(init_patch_version)
    
    # return current version
    return version


def set_version() -> str:
    """Explicit name for main wrapper."""
    return main()


############################
# Command line entry point #
############################
if __name__ == "__main__":
    main()
