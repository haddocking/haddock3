import os
import subprocess
import shutil
from haddock.modules.error import CNSError
from haddock.config import CNS_container, host_data_folder, container_data_folder


class CNS:
    """Wraps interacting with the CNS core facilities"""

    def __init__(self):
        self.stdout = None
        self.cns_exec = CNS_container
        if not os.path.exists(host_data_folder):
            os.mkdir(os.path.abspath(host_data_folder))

    def commit(self, input_file):
        """Pass the input defined in input_file to CNS"""
        shutil.copyfile(input_file, os.path.join(host_data_folder, os.path.basename(input_file)))
        p = subprocess.Popen([self.cns_exec.replace('input_file', os.path.join(container_data_folder, input_file))],
                             #stdin=os.path.join(container_data_folder, input_file), 
                             stdout=self.stdout, 
                             close_fds=True)
        out, error = p.communicate()

        p.kill()

        if error:
            raise CNSError(error)

        return out
