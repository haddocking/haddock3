import subprocess
import tempfile
import uuid
import shutil
import zipfile
from enum import Enum
from pathlib import Path
import glob

from haddock.libs.libsubprocess import CNSJob

SUBMIT_CMD = "/opt/diracos/bin/dirac-wms-job-submit"
STATUS_CMD = "/opt/diracos/bin/dirac-wms-job-status"
GETOUTPUT_CMD = "/opt/diracos/bin/dirac-wms-job-get-output"


def ping_dirac() -> bool:
    """Ping the Dirac caserver to check if it's reachable."""
    cmd = "/opt/diracos/bin/dirac-proxy-info"
    result = subprocess.run([cmd], shell=True, capture_output=True)
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
    def __init__(self, cnsjob: CNSJob, instructions: str) -> None:
        self.name = str(uuid.uuid4())
        self.loc = Path(tempfile.mkdtemp(prefix="haddock_grid_"))
        self.input_f = self.loc / "input.inp"
        self.expected_output_f = cnsjob.output_file
        self.id = None
        self.site = None
        self.stdout_f = None
        self.stderr_f = None
        self.output_f = None
        self.job_script = self.loc / "job.sh"
        self.jdl = self.loc / "job.jdl"
        self.status = JobStatus.UNKNOWN

        # TODO: Figure out all the files needed to be copied
        # Write the `.inp` file
        with open(self.input_f, "w") as f:
            f.write(cnsjob.input_file)

        # CREATE PAYLOAD
        with zipfile.ZipFile(f"{self.loc}/payload.zip", "w") as z:
            z.write(self.input_f, arcname="input.inp")
            z.write(cnsjob.cns_exec, arcname="cns")

        # CREATE JOB SCRIPT
        with open(self.job_script, "w") as f:
            f.write(instructions)

        # CREATE THE JDL
        with open(self.jdl, "w") as f:
            f.write(self.make_jdl())

    def submit(self) -> None:
        result = subprocess.run(
            [SUBMIT_CMD, f"{self.loc}/job.jdl"],
            shell=False,
            capture_output=True,
            text=True,
            cwd=self.loc,
        )

        # print(self.input_file)
        print(result.stdout)

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

        output_files = glob.glob(f"{self.loc}/{self.id}/*")

        print(f"Output files: {output_files}")
        return output_files

    def update_status(self):
        result = subprocess.run(
            [STATUS_CMD, str(self.id)], shell=False, capture_output=True, text=True
        )

        print(result.stdout)

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
            'InputSandbox = {"job.sh", "payload.zip"};',
            'OutputSandbox = {"job.out", "job.err",'
            + f' "{self.expected_output_f}"'
            + "};",
        ]
        return "[\n    " + "\n    ".join(jdl_lines) + "\n]\n"

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
