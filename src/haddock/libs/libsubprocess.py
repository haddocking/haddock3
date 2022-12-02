"""Run subprocess jobs."""
import os
import shlex
import subprocess
from contextlib import suppress
from pathlib import Path

from haddock.core.defaults import cns_exec as global_cns_exec
from haddock.core.exceptions import CNSRunningError, JobRunningError
from haddock.libs.libio import gzip_files


class BaseJob:
    """Base class for a subprocess job."""

    def __init__(self, input_, output, executable, *args):
        self.input = input_
        self.output = output
        self.executable = executable
        self.args = args

    def run(self):
        """Execute job in subprocess."""
        self.make_cmd()

        with open(self.output, 'w') as outf:
            p = subprocess.Popen(
                shlex.split(self.cmd),
                stdout=outf,
                close_fds=True,
                )
            out, error = p.communicate()

        p.kill()

        if error:
            raise JobRunningError(error)

        return out


def Job(BaseJob):
    """
    Instantiate a standard job.

    Runs with the following scheme:

        $ cmd ARGS INPUT
    """
    def make_cmd(self):
        """Execute subprocess job."""
        self.cmd = " ".join([
            os.fspath(self.executable),
            ''.join(map(str, self.args)),  # empty string if no args
            os.fspath(self.input),
            ])
        return


class JobInputFirst(BaseJob):
    """
    Instantiate a subprocess job with inverted args and input.

    Runs with the following scheme, INPUT comes first:

        $ cmd INPUT ARGS
    """

    def make_cmd(self):
        """Execute job in subprocess."""
        self.cmd = " ".join([
            os.fspath(self.executable),
            os.fspath(self.input),
            ''.join(map(str, self.args)),  # empty string if no args
            ])
        return


class CNSJob:
    """A CNS job script."""

    def __init__(
            self,
            input_file,
            output_file,
            envvars=None,
            cns_exec=None,
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

    def __repr__(self):
        return (
            f"CNSJob({self.input_file}, {self.output_file}, "
            f"envvars={self.envvars}, cns_exec={self.cns_exec})"
            )

    def __str__(self):
        return repr(self)

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
        if not cns_exec_path:
            cns_exec_path = global_cns_exec  # global cns_exec

        if not os.access(cns_exec_path, mode=os.X_OK):
            raise ValueError(
                f'{str(cns_exec_path)!r} binary file not found, '
                'or is not executable.'
                )

        self._cns_exec = cns_exec_path

    def run(
            self,
            compress_inp=False,
            compress_out=True,
            compress_seed=False,
            ):
        """
        Run this CNS job script.

        Parameters
        ----------
        compress_inp : bool
            Compress the *.inp file to '.gz' after the run. Defaults to
            ``False``.

        compress_out : bool
            Compress the *.out file to '.gz' after the run. Defaults to
            ``True``.

        compress_seed : bool
            Compress the *.seed file to '.gz' after the run. Defaults to
            ``False``.
        """
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

        if compress_inp:
            gzip_files(self.input_file, remove_original=True)

        if compress_out:
            gzip_files(self.output_file, remove_original=True)

        if compress_seed:
            with suppress(FileNotFoundError):
                gzip_files(
                    Path(Path(self.output_file).stem).with_suffix('.seed'),
                    remove_original=True)

        if error:
            raise CNSRunningError(error)

        return out
