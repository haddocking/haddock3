import os
from haddock.modules.cns.engine import CNS


class InputComposer(object):

    def __init__(self):
        """Initialize the composer

        define legacy_protocls and toppar loc"""

        self.cmd = ''
        self.legacy_protocol_path = '/'.join(
            os.path.dirname(os.path.realpath(__file__)).split('/')[:-2]) + '/legacy_protocols'
        self.toppar_path = '/'.join(
            os.path.dirname(os.path.realpath(__file__)).split('/')[:-2]) + '/toppar'

    def load_topology(self, topology_file):

        """Load topology"""

        self.cmd += '\n'
        self.cmd += 'topology\n'
        self.cmd += f'  @@{self.toppar_path}/{topology_file}\n'
        self.cmd += 'end\n'

    def load_parameter(self, parameter_file):

        """Load parameters"""

        self.cmd += '\n'
        self.cmd += 'parameter\n'
        self.cmd += f'  @@{self.toppar_path}/{parameter_file}\n'
        self.cmd += 'end\n'

    def load_coordinates(self, coordinate_file):

        """Define a segment and load coordinates"""

        self.cmd += '\n'
        self.cmd += 'segment\n'
        self.cmd += '  chain\n'
        self.cmd += f'    coordinates @@{coordinate_file}\n'
        self.cmd += '  end\n'
        self.cmd += 'end\n'

    def save_psf(self, output_file):

        """Output the resulting coordinates as psf"""

        self.cmd += '\n'
        self.cmd += f'write structure output={output_file} end\n'

    def patch_types_cg(self):

        """Change chemical types based on sidechain"""

        self.cmd += self.load_cns_protocol(f'{self.legacy_protocol_path}/patch-types-cg.cns')

    def patch_bb_cg(self):

        """Apply patch for CG beads"""

        self.cmd += self.load_cns_protocol(f'{self.legacy_protocol_path}/patch-bb-cg.cns')

    def build_missing_atoms(self):

        """Build missing atoms

        WARNING: The full protocol is not yet supported
        """

        self.cmd += self.load_cns_protocol(f'{self.legacy_protocol_path}/build-missing.cns')

    def apply_protein_break(self):

        """Detect protein chain breaks"""

        self.cmd += self.load_cns_protocol(f'{self.legacy_protocol_path}/prot_break.cns')

    def minimize(self, steps):

        """Run minimization"""

        self.cmd += f'minimize powell nstep={steps} drop=10.0 nprint=10 end\n'

    @staticmethod
    def load_cns_protocol(protocol_file):

        """Load legacy cns protocols"""

        cmd = ''
        with open(protocol_file) as f:
            for line in f.read():
                cmd += line
        return cmd


class InputGenerator(InputComposer):

    """Generate input files"""

    def __init__(self):
        super().__init__()
        self.cns_input_file = None

    def build_missing(self, coordinate_file, parameter_file, topology_file, output_file):
        self.cns_input_file = output_file

        super().load_topology(topology_file)
        super().load_parameter(parameter_file)
        super().load_coordinates(coordinate_file)
        super().apply_protein_break()
        super().build_missing_atoms()
        super().minimize(steps=100)
        super().save_psf(coordinate_file.replace('.pdb', '.psf'))

        self.save_input()

    def generate_cg_psf(self, coordinate_file, parameter_file, topology_file, output_file):
        self.cns_input_file = output_file

        super().load_topology(topology_file)
        super().load_parameter(parameter_file)
        super().load_coordinates(coordinate_file)
        super().patch_types_cg()
        super().patch_bb_cg()
        super().save_psf(coordinate_file.replace('.pdb', '.psf'))

        self.save_input()

    def generate_psf(self, coordinate_file, parameter_file, topology_file, output_file):
        self.cns_input_file = output_file

        super().load_topology(topology_file)
        super().load_parameter(parameter_file)
        super().load_coordinates(coordinate_file)
        super().save_psf(coordinate_file.replace('.pdb', '.psf'))

        self.save_input()

    def save_input(self):

        """Save the input file"""

        with open(self.cns_input_file, 'w') as f:
            f.write(self.cmd)
            f.write('\nstop')
        f.close()

    def execute(self):

        """Execute this input in CNS"""

        self.sanity_check()

        cns = CNS()
        inp = open(self.cns_input_file)
        cns.commit(inp)

    def sanity_check(self):

        """Check for multi-layer CNS scripting"""

        with open(self.cns_input_file) as f:
            data = f.read()
            for line in data.split('\n'):
                if '.cns' in line:
                    if '!' not in line:
                        print(f'ERROR: Multi-layer CNS scripting detected in {self.cns_input_file}, not supported.')
                        exit()
                    else:
                        pass
                else:
                    pass
