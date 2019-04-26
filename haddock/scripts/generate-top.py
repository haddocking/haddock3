import sys

try:
    from haddock.modules.cns.input import InputGenerator
except ModuleNotFoundError:
    print('ERROR: Modules could not be found')
    print('Make sure haddock3 is in your $PYTHONPATH')
    sys.exit(1)


def main():

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
    main()
