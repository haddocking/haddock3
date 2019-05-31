import sys
from haddock.modules.cns.input import InputGenerator


def topology():

    if len(sys.argv) != 2:
        print("Usage: python generate-top.py protein.pdb")
        sys.exit(1)

    prot = sys.argv[1]

    top = InputGenerator()
    top.build_missing(coordinate_file=f'{prot}',
                      topology_file='protein-allhdg5-4.top',
                      parameter_file='protein-allhdg5-4.param',
                      output_file='gen_top.inp')
    top.execute()


if __name__ == '__main__':
    topology()
