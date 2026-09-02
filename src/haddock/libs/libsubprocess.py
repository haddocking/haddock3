"""Run subprocess jobs."""

import os
import re
import shlex
import subprocess
from contextlib import suppress
from pathlib import Path

from haddock import log
from haddock.core.defaults import cns_exec as global_cns_exec
from haddock.core.exceptions import (
    CachedCNSFailure,
    CNSRunningError,
    JobRunningError,
)
from haddock.libs.libcache import (
    append_cache_record,
    append_debug_command,
    lookup_cache_record,
    source_cache_has_pdb_path,
    verify_and_restore,
)
from haddock.core.typing import Any, FilePath, Iterable, Optional, ParamDict
from haddock.gear.known_cns_errors import KNOWN_ERRORS as KNOWN_CNS_ERRORS
from haddock.libs.libcnsoutput import normalize_cns_pdb, normalize_cns_psf
from haddock.libs.libio import gzip_files
from haddock.libs.libseamless import (
    canonical_mapping_for_job,
    job_checksum,
    result_checksum_for_paths,
    stage_debug_synthesis,
)


CNS_DENORMAL_STDERR = (
    b"Note: The following floating-point exceptions are signalling: IEEE_DENORMAL\n"
)


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
        output_files: Optional[Iterable[FilePath]] = None,
        output_pdb_files: Optional[Iterable[FilePath]] = None,
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
        output_pdb_files : iterable
            PDB files expected from this CNS job. When provided, run-volatile
            CNS header lines are removed after successful execution.
        output_files : iterable
            Files expected from this CNS job. Files with ``.pdb`` or ``.psf``
            suffixes are normalized after successful execution.
        """
        self.input_file = input_file
        # Resolve in-memory CNS input consistently in worker processes.
        self.work_dir = Path.cwd().resolve()
        self.cache_context = None
        self.cache_debug = False
        self.cache_hit = False
        self.output_file = output_file
        self.error_file = error_file
        self.envvars = envvars
        self.cns_exec = cns_exec
        self.output_files = [Path(output_file) for output_file in output_files or []]
        self.output_pdb_files = list(
            dict.fromkeys(
                [
                    *(
                        Path(output_pdb_file)
                        for output_pdb_file in output_pdb_files or []
                    ),
                    *(
                        output_file
                        for output_file in self.output_files
                        if output_file.suffix.lower() == ".pdb"
                    ),
                ]
            )
        )
        for output_pdb_file in self.output_pdb_files:
            if output_pdb_file not in self.output_files:
                self.output_files.append(output_pdb_file)
        self._pending_partial_outputs: dict[Path, Path] = {}
        self._assert_declared_output_bindings()
        self._declared_input_script = self._input_script()

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
                f"{str(cns_exec_path)!r} binary file not found, or is not executable."
            )

        self._cns_exec = cns_exec_path

    def run(
        self,
        compress_inp: bool = False,
        compress_out: bool = True,
        compress_seed: bool = False,
        compress_err: bool = True,
    ) -> bytes:
        """Run this CNS job, serving it from a cache when one can supply it."""
        if self.cache_context is None:
            out = self._run_direct(
                compress_inp=compress_inp,
                compress_out=compress_out,
                compress_seed=compress_seed,
                compress_err=compress_err,
            )
            return out

        mapping, checksum = self._cache_mapping_and_checksum()
        cached_failure = False
        for source_index in self.cache_context.source_indexes:
            source_record = lookup_cache_record(source_index, checksum)
            if source_record is None:
                continue
            if source_record.result_checksum == "FAILED":
                cached_failure = True
                continue
            destinations = tuple(self._absolute_output(path) for path in mapping.output_paths)
            miss_reason = verify_and_restore(
                source_index,
                source_record,
                destinations,
                lambda paths: result_checksum_for_paths(mapping.canonical_output_names, paths),
            )
            if miss_reason is None:
                self.cache_hit = True
                return b""
            log.warning(
                "Cache miss for job %s from %s: %s",
                checksum,
                source_index.source_run,
                miss_reason,
            )
        if cached_failure:
            raise CachedCNSFailure(checksum)
        return self._run_direct(
            compress_inp=compress_inp,
            compress_out=compress_out,
            compress_seed=compress_seed,
            compress_err=compress_err,
        )


    def _run_direct(
        self,
        compress_inp: bool = False,
        compress_out: bool = True,
        compress_seed: bool = False,
        compress_err: bool = True,
    ) -> bytes:
        """
        Execute this CNS job script, without consulting any cache.

        Parameters
        ----------
        compress_inp : bool
            Compress the ``*.inp`` file to '.gz' after the run. Defaults to
            ``False``.

        compress_out : bool
            Compress the ``*.out`` file to '.gz' after the run. Defaults to
            ``True``.

        compress_seed : bool
            Compress the ``*.seed`` file to '.gz' after the run. Defaults to
            ``False``.
        """

        script = self.prepare_execution_input()

        p = subprocess.Popen(
            self.cns_exec,
            stdin=subprocess.PIPE,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            close_fds=True,
            env=self.envvars,
            cwd=self.work_dir,
        )
        out, error = p.communicate(input=script.encode())
        p.kill()

        # GNU Fortran emits this note for otherwise complete CNS calculations.
        error = (error or b"").replace(CNS_DENORMAL_STDERR, b"").strip()
        failed = bool(error) or self.contains_cns_stdout_error(out)
        if isinstance(p.returncode, int) and p.returncode != 0:
            failed = True
            error = error or f"CNS exited with status {p.returncode}".encode()
        if failed:
            # Write .err file (only if an error file was provided, otherwise
            # the diagnostic raise below would be masked by a TypeError)
            if self.error_file is not None:
                error_file = self._output_path(Path(self.error_file))
                with open(error_file, "wb+") as errf:
                    errf.write(out)
                # Compress it
                if compress_err:
                    gzip_files(
                        error_file,
                        remove_original=True,
                    )
            self._discard_partial_outputs()
            raise CNSRunningError(error or out)

        self.publish_outputs()

        if isinstance(self.input_file, Path) and compress_inp:
            gzip_files(self._output_path(self.input_file), remove_original=True)
        if isinstance(self.input_file, Path) and self.output_file is not None:
            output_file = self._output_path(Path(self.output_file))
            output_file.write_bytes(out)
            if compress_out:
                gzip_files(output_file, remove_original=True)
            if compress_seed:
                with suppress(FileNotFoundError):
                    seed_file = output_file.with_suffix(".seed")
                    if not seed_file.exists():
                        seed_file = self.work_dir / seed_file.name
                    gzip_files(
                        seed_file,
                        remove_original=True,
                    )

        # Return STDOUT
        return out

    def normalize_outputs(self) -> None:
        """Normalize known CNS-generated outputs."""
        self._normalize_paths(
            self._output_path(output_file) for output_file in self.output_files
        )

    @staticmethod
    def _normalize_paths(output_files: Iterable[Path]) -> None:
        """Normalize the supplied public or hidden CNS artifacts."""
        for output_file in output_files:
            suffix = output_file.suffix.lower()
            if suffix == ".pdb":
                normalize_cns_pdb(output_file)
            elif suffix == ".psf":
                normalize_cns_psf(output_file)

    def _input_script(self) -> Optional[str]:
        """Return the materialized CNS input, if it is available now."""
        if isinstance(self.input_file, str):
            return self.input_file
        path = self._output_path(self.input_file)
        if path.exists():
            return path.read_text(encoding="utf-8")
        return None

    def _assert_declared_output_bindings(self) -> None:
        """Reject a job whose metadata names outputs other than its CNS script."""
        script = self._input_script()
        if script is None:
            return
        expected = {
            ".pdb": self._output_path(next(iter(self.output_pdb_files))).resolve()
            if self.output_pdb_files
            else None,
            ".psf": next(
                (
                    self._output_path(output).resolve()
                    for output in self.output_files
                    if output.suffix.lower() == ".psf"
                ),
                None,
            ),
        }
        matches = {}
        for variable, suffix in (
            ("output_pdb_filename", ".pdb"),
            ("output_psf_filename", ".psf"),
        ):
            match = re.search(rf'\${variable}\s*=\s*"(?P<path>[^"]+)"', script)
            if match is None:
                continue
            matches[suffix] = match
            declared = self._output_path(Path(match.group("path"))).resolve()
            if expected[suffix] is None or declared != expected[suffix]:
                raise ValueError(
                    f"CNS ${variable}={match.group('path')!r} does not match "
                    f"the declared output for this job"
                )
        if matches:
            missing = [
                suffix
                for suffix, output in expected.items()
                if output and suffix not in matches
            ]
            if missing:
                raise ValueError(
                    "CNS input does not bind every declared output: "
                    + ", ".join(missing)
                )

    def _partial_output_script(self) -> tuple[str, dict[Path, Path]]:
        """Direct CNS to same-directory staging outputs until they are complete."""
        script = self._declared_input_script or self._input_script()
        if script is None:
            raise ValueError(f"CNS input {self.input_file!r} does not exist")
        self._declared_input_script = script
        partial_outputs: dict[Path, Path] = {}
        for variable in ("output_pdb_filename", "output_psf_filename"):
            match = re.search(rf'\${variable}\s*=\s*"(?P<path>[^"]+)"', script)
            if match is None:
                continue
            public = self._output_path(Path(match.group("path"))).resolve()
            partial = public.with_name(f"{public.stem}.partial{public.suffix}")
            partial_outputs[public] = partial
            partial_name = str(Path(match.group("path")).with_name(partial.name))
            script = (
                script[: match.start("path")]
                + partial_name
                + script[match.end("path") :]
            )
        return script, partial_outputs

    def prepare_execution_input(self) -> str:
        """Return and retain the exact CNS script that will be executed.

        Batch backends call this before writing their temporary input file; the
        matching :meth:`publish_outputs` call then validates and publishes the
        complete output set after the process has succeeded. Path-backed inputs
        are updated so the retained ``.inp`` and CNS stdout describe the same
        execution.
        """
        script, partial_outputs = self._partial_output_script()
        for partial in partial_outputs.values():
            partial.unlink(missing_ok=True)
        self._pending_partial_outputs = partial_outputs
        if isinstance(self.input_file, Path):
            self._output_path(self.input_file).write_text(script, encoding="utf-8")
        return script

    def publish_outputs(self, check_output_log: bool = False) -> None:
        """Publish a prepared complete output set, or normalize legacy outputs."""
        if check_output_log and self.output_file is not None:
            output_file = self._output_path(Path(self.output_file))
            if output_file.exists() and self.contains_cns_stdout_error(
                output_file.read_bytes()
            ):
                self._discard_partial_outputs()
                raise CNSRunningError(f"CNS reported an error in {output_file}")
        if self._pending_partial_outputs:
            self._publish_partial_outputs(self._pending_partial_outputs)
            self._pending_partial_outputs = {}
        else:
            self.normalize_outputs()

    def _publish_partial_outputs(self, partial_outputs: dict[Path, Path]) -> None:
        """Validate and atomically publish a complete normalized output set."""
        missing = [
            path
            for path in partial_outputs.values()
            if not path.is_file() or path.stat().st_size == 0
        ]
        if missing:
            self._discard_partial_outputs(partial_outputs)
            raise CNSRunningError(
                "CNS did not produce complete outputs: "
                + ", ".join(str(path) for path in missing)
            )
        self._normalize_paths(partial_outputs.values())
        for public, partial in partial_outputs.items():
            os.replace(partial, public)

    def _discard_partial_outputs(
        self, partial_outputs: Optional[dict[Path, Path]] = None
    ) -> None:
        """Remove unpublished CNS artifacts after a failed or retried job."""
        pending = partial_outputs or self._pending_partial_outputs
        for partial in pending.values():
            partial.unlink(missing_ok=True)
        self._pending_partial_outputs = {}

    def _output_path(self, output_file: Path) -> Path:
        """Resolve job outputs relative to the directory that created the job."""
        if output_file.is_absolute():
            return output_file
        return self.work_dir / output_file

    def has_cached_output_file(self) -> bool:
        """Return whether a CACHE record declares this job's PDB path.

        This is intentionally a path-only scheduler hint.  In particular, it
        must not build a canonical mapping or checksum on the scheduler's
        main thread.
        """
        if self.cache_context is None or not self.cache_context.source_indexes:
            return False
        try:
            self._validate_cache_outputs()
            output = self._absolute_output(self.output_pdb_files[0])
            relative = output.relative_to(self.cache_context.current_run.resolve())
            return any(
                source_cache_has_pdb_path(source_index, relative)
                for source_index in self.cache_context.source_indexes
            )
        except (OSError, ValueError):
            return False

    def _cache_mapping_and_checksum(self):
        mapping = getattr(self, "_cached_mapping", None)
        checksum = getattr(self, "_cached_job_checksum", None)
        if mapping is None or checksum is None:
            mapping = canonical_mapping_for_job(self)
            checksum = job_checksum(mapping)
            self._cached_mapping = mapping
            self._cached_job_checksum = checksum
        return mapping, checksum

    def _absolute_output(self, output: Path) -> Path:
        """Resolve a declared output independently of a worker's CWD."""
        return self._output_path(output).resolve()

    def _psf_output(self) -> Path | None:
        psf_outputs = [output for output in self.output_files if output.suffix == ".psf"]
        return self._absolute_output(psf_outputs[0]) if psf_outputs else None

    def _validate_cache_outputs(self) -> None:
        if len(self.output_pdb_files) != 1:
            raise CNSRunningError("CNS cache requires exactly one declared PDB output.")
        if len([path for path in self.output_files if path.suffix == ".psf"]) > 1:
            raise CNSRunningError("CNS cache requires at most one declared PSF output.")

    def _missing_outputs(self) -> list[Path]:
        """Return declared outputs that are not yet visible in the job directory."""
        return [
            path
            for path in self.output_files
            if not self._absolute_output(path).is_file()
        ]

    def cache_outputs_present(self) -> bool:
        """Return whether every predeclared cache artifact is visible."""
        self._validate_cache_outputs()
        return not self._missing_outputs()

    def prepare_cache_success_record(self, checksum: str | None = None, mapping=None):
        """Normalize visible outputs and calculate their cache record payload."""
        self._validate_cache_outputs()
        if checksum is None:
            mapping = canonical_mapping_for_job(self)
            checksum = job_checksum(mapping)
        self.normalize_outputs()
        psf_output = self._psf_output()
        result = result_checksum_for_paths(
            self._cache_output_names(),
            (
                self._absolute_output(self.output_pdb_files[0]),
                *((psf_output,) if psf_output is not None else ()),
            ),
        )
        return mapping, checksum, result, psf_output

    def append_cache_success_record(self, payload) -> None:
        """Append a precomputed cache payload from the writer thread."""
        mapping, checksum, result, psf_output = payload
        append_cache_record(
            self.cache_context,
            checksum,
            result,
            self._absolute_output(self.output_pdb_files[0]),
            psf_output,
        )
        if self.cache_debug:
            mapping = mapping or canonical_mapping_for_job(self)
            self._append_cache_debug_command(mapping, checksum, result)

    def write_cache_success_record(self) -> None:
        """Synchronously prepare and append a completed job's cache record."""
        self.append_cache_success_record(self.prepare_cache_success_record())

    def write_cache_failure_record(self) -> None:
        """Record a job whose complete declared artifact set never appeared."""
        self._validate_cache_outputs()
        mapping = canonical_mapping_for_job(self)
        checksum = job_checksum(mapping)
        append_cache_record(
            self.cache_context,
            checksum,
            "FAILED",
            self._absolute_output(self.output_pdb_files[0]),
            self._psf_output(),
        )
        if self.cache_debug:
            self._append_cache_debug_command(mapping, checksum, "FAILED")

    def cache_writer_completion(self):
        """Return the worker's already-required job identity to the writer.

        The writer still normalizes and checksums visible output bytes itself;
        only the identity used for the cache lookup is reused here.
        """
        return (
            self.cache_writer_id,
            getattr(self, "_cached_job_checksum", None),
            getattr(self, "_cached_mapping", None) if self.cache_debug else None,
        )

    def _cache_output_names(self) -> tuple[str, ...]:
        """Return fixed canonical names for the declared cache artifact shape."""
        return ("canonical-output.pdb",) + (
            ("canonical-output.psf",) if self._psf_output() is not None else ()
        )

    def _append_cache_debug_command(self, mapping, checksum: str, result: str) -> None:
        if self.cache_debug:
            command = stage_debug_synthesis(mapping, checksum).command
            append_debug_command(self.cache_context, checksum, result, command)

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
