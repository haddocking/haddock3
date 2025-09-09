import os
import subprocess
import tempfile
import uuid
import zipfile
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

SUBMIT_CMD = "/opt/diracos/bin/dirac-wms-job-submit"
STATUS_CMD = "/opt/diracos/bin/dirac-wms-job-status"
GETOUTPUT_CMD = "/opt/diracos/bin/dirac-wms-job-get-output"


def ping_dirac() -> bool:
    """Ping the Dirac caserver to check if it's reachable."""
    cmd = "/opt/diracos/bin/dirac-proxy-info"
    result = subprocess.run([cmd], shell=True, capture_output=True)
    print(result.stdout.decode())
    return result.returncode == 0


class JobStatus(Enum):
    WAITING = "Waiting"
    RUNNING = "Running"
    UNKNOWN = "Unknown"
    DONE = "Done"
    MATCHED = "Matched"
    COMPLETING = "Completing"

    @classmethod
    def from_string(cls, value):
        """Convert string to JobStatus enum."""
        value = value.strip().lower()
        for status in cls:
            if status.value.lower() == value:
                return status
        return cls.UNKNOWN


class GridJob:

    def __init__(self, instructions: str) -> None:
        self.loc = tempfile.mkdtemp()
        self.id = None
        self.site = None
        self.stdout_f = None
        self.stderr_f = None
        self.output_f = None
        self.job = self.loc + "/job.sh"
        self.jdl = self.loc + "/job.jdl"
        self.status = JobStatus.UNKNOWN
        self.name = str(uuid.uuid4())

        # DEBUG:
        with open(self.job, "w") as f:
            f.write(instructions)

        with open(self.jdl, "w") as f:
            f.write(self.make_jdl())

    def run(self) -> None:
        result = subprocess.run(
            [SUBMIT_CMD, f"{self.loc}/job.jdl"],
            shell=False,
            capture_output=True,
            text=True,
            cwd=self.loc,
        )

        self.id = int(result.stdout.split()[-1])
        self.update_status()

    def retrieve_output(self) -> list[Path]:
        # TODO: Add error handling
        subprocess.run(
            [GETOUTPUT_CMD, str(self.id)],
            shell=False,
            capture_output=True,
            text=True,
            cwd=self.loc,
        )

        self.stdout_f = Path(f"{self.loc}/{self.id}/job.out")
        self.stderr_f = Path(f"{self.loc}/{self.id}/job.err")
        self.output_f = Path(f"{self.loc}/{self.id}/output.zip")

        output_files = []
        with zipfile.ZipFile(self.output_f, "r") as z:
            z.extractall(self.loc)
            for f in z.namelist():
                output_files.append(Path(f"{self.loc}/{self.id}/{f}"))

        print(f"Output files: {output_files}")
        return output_files

    def update_status(self):
        result = subprocess.run(
            [STATUS_CMD, str(self.id)], shell=False, capture_output=True, text=True
        )

        output_dict = self.parse_output(result.stdout)

        self.id = output_dict["JobID"]
        self.status = JobStatus.from_string(output_dict["Status"])
        self.site = output_dict.get("Site", "Unknown")

    def make_jdl(self) -> str:
        jdl_lines = [
            f'JobName = "{self.name}";',
            'Executable = "job.sh";',
            'Arguments = "";',
            'StdOutput = "job.out";',
            'StdError = "job.err";',
            'JobType = "WeNMR";',
            'InputSandbox = {"job.sh"};',
            'OutputSandbox = {"job.out", "job.err", "output.zip"};',
        ]
        return "[\n    " + "\n    ".join(jdl_lines) + "\n]"

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
