"""Class that represents a Haddock model"""


class HaddockModel:

    def __init__(self, pdb_f):
        self.energies = self._load_energies(pdb_f)

    @staticmethod
    def _load_energies(pdb_f):
        energy_dic = {}
        with open(pdb_f) as fh:
            for line in fh.readlines():
                if line.startswith('REMARK'):
                    # TODO: use regex to do this
                    if 'energies' in line:
                        energy_values = list(map(float, line.rstrip().split(':')[-1].split(',')))
                        total,bonds,angles,improper,dihe,vdw,elec,air,cdih,coup,rdcs,vean,dani,xpcs,rg = energy_values
                        energy_dic['total'] = total
                        energy_dic['bonds'] = bonds
                        energy_dic['angles'] = angles
                        energy_dic['improper'] = improper
                        energy_dic['dihe'] = dihe
                        energy_dic['vdw'] = vdw
                        energy_dic['elec'] = elec
                        energy_dic['air'] = air
                        energy_dic['cdih'] = cdih
                        energy_dic['coup'] = coup
                        energy_dic['rdcs'] = rdcs
                        energy_dic['vean'] = vean
                        energy_dic['dani'] = dani
                        energy_dic['xpcs'] = xpcs
                        energy_dic['rg'] = rg
                    if 'buried surface area' in line:
                        bsa = float(line.rstrip().split(':')[-1])
                        energy_dic['bsa'] = bsa
                    if 'Desolvation energy' in line:
                        desolv = float(line.rstrip().split(':')[-1])
                        energy_dic['desolv'] = desolv

        return energy_dic

    def calc_haddock_score(self, **weights):
        """Calculate the haddock score based on the weights and energies."""
        weighted_terms = []
        for key in weights:
            component_id = key.split('_')[1]
            weight = weights[key]
            value = self.energies[component_id]
            weighted_terms.append(value * weight)

        # the haddock score is simply the sum of the weighted terms
        haddock_score = sum(weighted_terms)
        return haddock_score
