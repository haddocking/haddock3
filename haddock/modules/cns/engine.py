import os
import subprocess
from haddock.modules.error import CNSError
from haddock.config import CNS_bin

class CNS:
    """Wraps interacting with the CNS core facilities"""

    def __init__(self):
        self.stdout = None
        self.cns_exec = CNS_bin

    def commit(self, input_file):
        """Pass the input defined in input_file to CNS"""
        print(self.cns_exec)
        print(input_file)
        p = subprocess.Popen([self.cns_exec],
                             stdin=input_file, stdout=self.stdout, close_fds=True)
        out, error = p.communicate()

        p.kill()

        if error:
            raise CNSError(error)

        return out
