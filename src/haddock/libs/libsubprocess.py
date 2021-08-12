"""Running CNS scripts"""
import subprocess
from haddock.error import CNSRunningError
from haddock.libs.libparallel import Scheduler
from haddock.defaults import CNS_EXE, NUM_CORES


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
