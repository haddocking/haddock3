"""Running CNS scripts"""
import os
import shlex
import subprocess

from haddock.core.defaults import CNS_EXE, NUM_CORES
from haddock.core.exceptions import CNSRunningError, JobRunningError
from haddock.libs.libparallel import Scheduler


class Job:
    """A job to be executed by the engine"""
    def __init__(self, input, output, executable, *args):
        self.input = input
        self.output = output
        self.executable = executable
        self.args = args

    def run(self):
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
    """A CNS job script"""
    def __init__(self, input_file, output_file, cns_folder='.',
                 cns_exec=CNS_EXE):
        """
        :param input_file: input CNS script
        :param output_file: CNS output
        :cns_folder: absolute execution path
        :cns_exec: CNS binary including absolute path
        """
        self.input_file = input_file
        self.output_file = output_file
        self.cns_folder = cns_folder
        self.cns_exec = cns_exec

    def run(self):
        """Run this CNS job script"""
        with open(self.input_file) as inp:
            with open(self.output_file, 'w+') as outf:
                env = {'RUN': self.cns_folder}
                p = subprocess.Popen(self.cns_exec,
                                     stdin=inp,
                                     stdout=outf,
                                     close_fds=True,
                                     env=env)
                out, error = p.communicate()
                p.kill()
        if error:
            raise CNSRunningError(error)
        return out
