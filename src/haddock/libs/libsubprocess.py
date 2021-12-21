"""Run subprocess jobs."""
import os
import shlex
import subprocess
from pathlib import Path

from haddock.core.defaults import cns_exec
from haddock.core.exceptions import CNSRunningError, JobRunningError


class Job:
    """A job to be executed by the engine."""

    def __init__(self, input, output, executable, *args):
        self.input = input
        self.output = output
        self.executable = executable
        self.args = args

    def run(self):
        """Execute subprocess job."""
        cmd = " ".join([
            os.fspath(self.executable),
            ''.join(map(str, self.args)),  # empty string if no args
            os.fspath(self.input),
            ])

        with open(self.output, 'w') as outf:
            p = subprocess.Popen(shlex.split(cmd),
                                 stdout=outf,
                                 close_fds=True)
            out, error = p.communicate()

        p.kill()

        if error:
            raise JobRunningError(error)

        return out


class CNSJob:
    """A CNS job script."""

    def __init__(
            self,
            input_file,
            output_file,
            envvars=None,
            ):
        """
        CNS subprocess.

        To execute the job, call the `.run()` method.

        Parameters
        ----------
        input_file : str or pathlib.Path
            The path to the .inp CNS file.

        output_file : str or pathlib.Path
            The path to the .out CNS file, where the standard output
            will be saved.

        envvars : dict
            A dictionary containing the environment variables needed for
            the CNSJob. These will be passed to subprocess.Popen.env
            argument.
        """
        self.input_file = input_file
        self.output_file = output_file
        self.envvars = envvars
        self.cns_exec = cns_exec

    @property
    def envvars(self):
        """CNS environment vars."""
        return self._envvars

    @envvars.setter
    def envvars(self, envvars):
        """CNS environment vars."""
        self._envvars = envvars or {}
        if not isinstance(self._envvars, dict):
            raise ValueError('`envvars` must be a dictionary.')

        for k, v in self._envvars.items():
            if isinstance(v, Path):
                self._envvars[k] = str(v)

    @property
    def cns_exec(self):
        """CNS executable path."""
        return self._cns_exec

    @cns_exec.setter
    def cns_exec(self, cns_exec_path):
        if cns_exec_path is None:
            cns_exec_path = cns_exec  # global cns_exec

        if not os.access(cns_exec_path, mode=os.X_OK):
            raise ValueError(
                f'{str(cns_exec_path)!r} binary file not found, '
                'or is not executable.'
                )

        self._cns_exec = cns_exec_path

    def run(self):
        """Run this CNS job script."""
        with open(self.input_file) as inp, \
                open(self.output_file, 'w+') as outf:

            p = subprocess.Popen(
                self.cns_exec,
                stdin=inp,
                stdout=outf,
                stderr=subprocess.PIPE,
                close_fds=True,
                env=self.envvars,
                )

            out, error = p.communicate()
            p.kill()

        if error:
            raise CNSRunningError(error)

        return out
