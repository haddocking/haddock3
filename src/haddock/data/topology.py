"""HADDOCK topologies"""
from haddock.defaults import Default


class Topology:

    @staticmethod
    def get_supported():
        """Read the topology file and identify which data is supported"""
        supported = []
        with open(Default.TOPOLOGY_FILE) as input_handler:
            for line in input_handler:
                if "resi" in line[:4].casefold():
                    res = line.split()[1]
                    supported.append(res)
        return supported
