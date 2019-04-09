import os
import subprocess


class CNS:

    """CNS communicator class"""

    def __init__(self):
        self.stdout = None
        try:
            self.cns_exec = os.environ['CNS_EXE']
        except KeyError:
            print('System variable CNS_EXE not defined')
            exit()

    def commit(self, input_file):

        """Pass the input to CNS"""

        p = subprocess.Popen([self.cns_exec],
                             stdin=input_file, stdout=self.stdout, close_fds=True)

        out, error = p.communicate()

        p.kill()

        if error:
            print(f'Error: {error}')
