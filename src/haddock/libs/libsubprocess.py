"""Run subprocess jobs."""

import os
import shlex
import subprocess
from contextlib import suppress
from pathlib import Path

from haddock.core.defaults import cns_exec as global_cns_exec
from haddock.core.exceptions import (
    CNSRunningError,
    JobRunningError,
    KnownCNSError,
    )
from haddock.core.typing import Any, FilePath, Optional, ParamDict
from haddock.gear.known_cns_errors import KNOWN_ERRORS as KNOWN_CNS_ERRORS
from haddock.libs.libio import gzip_files


class BaseJob:
    """Base class for a subprocess job."""

    def __init__(
        self, input_: Path, output: Path, executable: Path, *args: Any
    ) -> None:
        self.input = input_
        self.output = output
        self.executable = executable
        self.args = args
        self.cmd: str

    def run(self) -> bytes:
        """Execute job in subprocess."""
        self.make_cmd()

        with open(self.output, "w") as outf:
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

    def make_cmd(self) -> None:
        """Execute subprocess job."""
        raise NotImplementedError()


class Job(BaseJob):
    """
    Instantiate a standard job.

    Runs with the following scheme:

        $ cmd ARGS INPUT
    """

    def make_cmd(self) -> None:
        """Execute subprocess job."""
        self.cmd = " ".join(
            [
                os.fspath(self.executable),
                "".join(map(str, self.args)),  # empty string if no args
                os.fspath(self.input),
            ]
        )
        return


class JobInputFirst(BaseJob):
    """
    Instantiate a subprocess job with inverted args and input.

    Runs with the following scheme, INPUT comes first:

        $ cmd INPUT ARGS
    """

    def make_cmd(self) -> None:
        """Execute job in subprocess."""
        self.cmd = " ".join(
            [
                os.fspath(self.executable),
                os.fspath(self.input),
                "".join(map(str, self.args)),  # empty string if no args
            ]
        )
        return


class CNSJob:
    """A CNS job script."""

    def __init__(
        self,
        input_file: FilePath,
        output_file: Optional[FilePath] = None,
        error_file: Optional[FilePath] = None,
        envvars: Optional[ParamDict] = None,
        cns_exec: Optional[FilePath] = None,
    ) -> None:
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
        self.error_file = error_file
        self.envvars = envvars
        self.cns_exec = cns_exec

    def __repr__(self) -> str:
        _input_file = self.input_file
        if isinstance(self.input_file, str):
            _input_file = "IO Stream"
        return (
            f"CNSJob({_input_file}, {self.output_file}, "
            f"envvars={self.envvars}, cns_exec={self.cns_exec})"
        )

    def __str__(self) -> str:
        return repr(self)

    @property
    def envvars(self) -> ParamDict:
        """CNS environment vars."""
        return self._envvars

    @envvars.setter
    def envvars(self, envvars: Optional[ParamDict]) -> None:
        """CNS environment vars."""
        self._envvars = envvars or {}
        if not isinstance(self._envvars, dict):
            raise ValueError("`envvars` must be a dictionary.")

        for k, v in self._envvars.items():
            if isinstance(v, Path):
                self._envvars[k] = str(v)

    @property
    def cns_exec(self) -> FilePath:
        """CNS executable path."""
        return self._cns_exec

    @cns_exec.setter
    def cns_exec(self, cns_exec_path: Optional[FilePath]) -> None:
        if not cns_exec_path:
            cns_exec_path = global_cns_exec  # global cns_exec

        if not os.access(cns_exec_path, mode=os.X_OK):
            raise ValueError(
                f"{str(cns_exec_path)!r} binary file not found, "
                "or is not executable."
            )

        self._cns_exec = cns_exec_path

    def run(
        self,
        compress_inp: bool = False,
        compress_out: bool = True,
        compress_seed: bool = False,
        compress_err: bool = True,
    ) -> bytes:
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

        if isinstance(self.input_file, str):
            p = subprocess.Popen(
                self.cns_exec,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                close_fds=True,
                env=self.envvars,
            )
            out, error = p.communicate(input=self.input_file.encode())
            p.kill()

        elif isinstance(self.input_file, Path) and self.output_file is not None:
            with open(self.input_file) as inp:
                p = subprocess.Popen(
                    self.cns_exec,
                    stdin=inp,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    close_fds=True,
                    env=self.envvars,
                )
                out, error = p.communicate()
                p.kill()
                # Write out file
                with open(self.output_file, "wb+") as outf:
                    outf.write(out)

            if compress_inp:
                gzip_files(self.input_file, remove_original=True)

            if compress_out:
                gzip_files(self.output_file, remove_original=True)

            if compress_seed:
                with suppress(FileNotFoundError):
                    gzip_files(
                        Path(Path(self.output_file).stem).with_suffix(".seed"),
                        remove_original=True,
                    )

        # If undetected error or detect an error in the STDOUT
        if error or self.contains_cns_stdout_error(out):
            # Write .err file
            with open(self.error_file, "wb+") as errf:
                errf.write(out)
            # Compress it
            if compress_err:
                gzip_files(self.error_file, remove_original=True)
            if error:
                raise CNSRunningError(error)

        # Return STDOUT
        return out

    @staticmethod
    def contains_cns_stdout_error(out: bytes) -> bool:
        # Decode end of STDOUT
        # Search in last 24000 characters (300 lines * 80 characters)
        sout = out[-24000:].split(bytes(os.linesep, "utf-8"))
        # Reverse loop on lines (read backward)
        for bytes_line in reversed(sout):
            line = bytes_line.decode("utf-8")
            # This checks for an unknown CNS error
            # triggered when CNS is about to crash due to internal error
            if "^^^^^" in line:
                return True
            # Check if a known error is found
            elif any([error in line for error in KNOWN_CNS_ERRORS.keys()]):
                return True
        return False
