import re
import shutil
import subprocess
import tempfile
import uuid
import zipfile
from enum import Enum
from pathlib import Path
import gzip

from haddock.libs.libsubprocess import CNSJob

PROXYINFO_CMD = "/opt/diracos/bin/dirac-proxy-info"
SUBMIT_CMD = "/opt/diracos/bin/dirac-wms-job-submit"
STATUS_CMD = "/opt/diracos/bin/dirac-wms-job-status"
GETOUTPUT_CMD = "/opt/diracos/bin/dirac-wms-job-get-output"


def ping_dirac() -> bool:
    """Ping the Dirac caserver to check if it's reachable."""
    result = subprocess.run([PROXYINFO_CMD], shell=True, capture_output=True)
    return result.returncode == 0


class JobStatus(Enum):
    WAITING = "Waiting"
    RUNNING = "Running"
    UNKNOWN = "Unknown"
    DONE = "Done"
    MATCHED = "Matched"
    COMPLETING = "Completing"
    FAILED = "Failed"

    @classmethod
    def from_string(cls, value):
        """Convert string to JobStatus enum."""
        value = value.strip().lower()
        for status in cls:
            if status.value.lower() == value:
                return status
        return cls.UNKNOWN


class GridJob:
    def __init__(self, cnsjob: CNSJob) -> None:
        self.name = str(uuid.uuid4())
        self.loc = Path(tempfile.mkdtemp(prefix="haddock_grid_"))
        self.wd = Path.cwd()
        self.input_f = self.loc / "input.inp"
        self.cns_out_f = cnsjob.output_file
        self.output_f = None
        self.id = None
        self.site = None
        self.stdout_f = None
        self.stderr_f = None
        self.output_f = None
        self.job_script = self.loc / "job.sh"
        self.jdl = self.loc / "job.jdl"
        self.status = JobStatus.UNKNOWN
        self.payload_fnames = []

        # NOTE: `cnsjob.input_file` can be a string or a path, handle the polymorfism here
        if isinstance(cnsjob.input_file, Path):
            self.input_f = cnsjob.input_file

        elif isinstance(cnsjob.input_file, str):
            with open(self.input_f, "w") as f:
                f.write(cnsjob.input_file)

        # CREATE JOB FILE
        self.create_job_script()

        # CREATE PAYLOAD
        self.prepare_payload(
            cns_exec_path=Path(cnsjob.cns_exec),
            cns_script_path=Path(cnsjob.envvars["MODULE"]),
            toppar_path=Path(cnsjob.envvars["TOPPAR"]),
        )

        # CREATE THE JDL
        self.create_jdl()

    def prepare_payload(
        self, cns_exec_path: Path, cns_script_path: Path, toppar_path: Path
    ):
        # CRAWL THE INPUT FILE AND REPLICATE THE FILE STRUCTURE
        self.process_input_f()
        #  CNS SCRIPTS
        for f in cns_script_path.glob("*"):
            self.payload_fnames.append(Path(f))
        # TOPPAR (including subdirectories)
        for f in toppar_path.rglob("*"):
            if f.is_file():  # Only add files, not directories
                self.payload_fnames.append(f)
        # CREATE PAYLOAD
        with zipfile.ZipFile(f"{self.loc}/payload.zip", "w") as z:
            z.write(self.input_f, arcname="input.inp")
            z.write(cns_exec_path, arcname="cns")
            for f in set(self.payload_fnames):
                # Preserve the relative path structure from toppar_path
                if f.is_relative_to(toppar_path):
                    relative_path = f.relative_to(toppar_path)
                    z.write(f, arcname=str(relative_path))
                else:
                    z.write(f, arcname=Path(f).name)

    def create_job_script(self):
        instructions = """#!/bin/bash
export MODULE=./
export TOPPAR=./
unzip payload.zip
./cns < input.inp > cns.log
"""
        # CREATE JOB SCRIPT
        with open(self.job_script, "w") as f:
            f.write(instructions)

    def process_input_f(self):
        """Read the `.inp` file, edit the paths and copy the files to the `loc`"""
        original_paths = []
        output_filename = None

        # FIXME: This works but it's very esoteric, refactor for clarity
        patterns = [
            (
                r"@@([^@\n]*\.(?:pdb|psf)[^@\n]*)",
                lambda m: (
                    original_paths.append(Path(m.group(1))),
                    f"@@{Path(m.group(1)).name}",
                )[1],
            ),
            (
                r'(\$input_p(?:db|sf)_filename[^=]*=\s*")([^"]*\.(?:pdb|psf)[^"]*)"',
                lambda m: (
                    original_paths.append(Path(m.group(2))),
                    f'{m.group(1)}{Path(m.group(2)).name}"',
                )[1],
            ),
        ]

        # Pattern to capture output filename
        output_pattern = r'\$output_pdb_filename\s*=\s*"?([^"\s\n]+\.pdb)"?'

        # Read the file first
        with open(self.input_f, "r") as f:
            lines = f.readlines()

        # Write the modified lines back
        with open(self.input_f, "w") as f:
            for line in lines:
                # Check for output filename before modifying the line
                output_match = re.search(output_pattern, line)
                if output_match:
                    output_filename = output_match.group(1)

                # original_line = line
                for pattern, replacement in patterns:
                    line = re.sub(pattern, replacement, line)
                # if line != original_line:
                #     print(f"old line: {original_line.strip()}")
                #     print(f"new line: {line.strip()}")
                #     print("-----------------------")
                f.write(line)

        for src_path in original_paths:
            dst_path = self.loc / src_path.name
            shutil.copy(src_path, dst_path)
            self.payload_fnames.append(dst_path)

        # Store the output filename for later use
        if output_filename:
            self.output_f = Path(output_filename).name
            print(f"Output file will be: {self.output_f}")

        return original_paths, output_filename

    def submit(self) -> None:
        result = subprocess.run(
            [SUBMIT_CMD, f"{self.loc}/job.jdl"],
            shell=False,
            capture_output=True,
            text=True,
            cwd=self.loc,
        )

        self.id = int(result.stdout.split()[-1])
        self.update_status()

    def retrieve_output(self):
        # TODO: Add error handling
        subprocess.run(
            [GETOUTPUT_CMD, str(self.id)],
            shell=False,
            capture_output=True,
            text=True,
            cwd=self.loc,
        )

        # NOTE: This is the output of the `job.sh`
        self.stdout_f = Path(f"{self.loc}/{self.id}/job.out")
        self.stderr_f = Path(f"{self.loc}/{self.id}/job.err")

        # Copy the output to the expected location
        src = Path(f"{self.loc}/{self.id}/{self.output_f}")
        dst = Path(self.wd / f"{self.output_f}")
        shutil.copy(src, dst)

        src = Path(f"{self.loc}/{self.id}/cns.log")
        with open(src, "rb") as f_in, gzip.open(f"{self.cns_out_f}.gz", "wb") as f_out:
            f_out.writelines(f_in)

        # ALL DONE, cleanup
        self.clean()

    def update_status(self):
        result = subprocess.run(
            [STATUS_CMD, str(self.id)], shell=False, capture_output=True, text=True
        )

        output_dict = self.parse_output(result.stdout)

        self.id = output_dict["JobID"]
        self.status = JobStatus.from_string(output_dict["Status"])
        self.site = output_dict.get("Site", "Unknown")

    def create_jdl(self):
        jdl_lines = [
            f'JobName = "{self.name}";',
            'Executable = "job.sh";',
            'Arguments = "";',
            'StdOutput = "job.out";',
            'StdError = "job.err";',
            'JobType = "WeNMR";',
            'Site = "EGI.SARA.nl";',
            'InputSandbox = {"job.sh", "payload.zip"};',
            'OutputSandbox = {"job.out", "job.err", "cns.log",'
            + f' "{self.output_f}"'
            + "};",
        ]

        jdl_string = "[\n    " + "\n    ".join(jdl_lines) + "\n]\n"

        with open(self.jdl, "w") as f:
            f.write(jdl_string)

    def clean(self):
        shutil.rmtree(self.loc)

    @staticmethod
    def parse_output(output_str: str) -> dict:
        items = output_str.replace(";", "")
        status_dict = {}
        for item in items.split(" "):
            if "=" in item:
                key, value = item.split("=", 1)
                status_dict[key.strip()] = value.strip()
        return status_dict

    def __repr__(self) -> str:
        return f"ID: {self.id} Name: {self.name} Status: {self.status.value} Site: {self.site}"


class GRIDScheduler:
    def __init__(self, tasks: list[CNSJob], params: dict):
        self.tasks = tasks
        self.params = params

    def run(self) -> None:
        """Execute the tasks."""
        queue = {}

        # Convert CNSJobs to GridJobs
        queue = {GridJob(t): False for t in self.tasks}

        # Submit jobs
        for job in queue.keys():
            job.submit()

        # Wait for jobs to finish
        complete = False
        while not complete:
            for job, done in queue.items():
                if not done:
                    job.update_status()
                    if job.status == JobStatus.DONE:
                        queue[job] = True
            complete = all(queue.values())

        # Retrieve outputs
        for job in queue.keys():
            job.retrieve_output()
