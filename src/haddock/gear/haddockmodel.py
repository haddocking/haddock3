"""Represent an Haddock model."""

from haddock.core.typing import FilePath


class HaddockModel:
    """Represent HADDOCK model."""

    def __init__(self, pdb_f: FilePath) -> None:
        remarks = self._load_remarks(pdb_f)
        self.energies = self._load_energies(remarks)
        self.interface_energies = self._load_per_interface_energies(remarks)
    
    @staticmethod
    def _load_remarks(pdb_f: FilePath) -> list[str]:
        """Load remark lines from PDB file

        Parameters
        ----------
        pdb_f : FilePath
            Path to a PDB file

        Returns
        -------
        remarks : list[str]
            List of PDB lines containing remarks.
        """
        remarks: list[str] = []
        with open(pdb_f) as fh:
            for line in fh.readlines():
                if line.startswith('REMARK'):
                    remarks.append(line)
        return remarks

    @staticmethod
    def _load_energies(remarks: list[str]) -> dict[str, float]:
        """Load HADDOCK energy terms from PDB remarks.

        Parameters
        ----------
        remarks : list[str]
            List of PDB lines containing remarks.

        Returns
        -------
        energy_dic : dict[str, float]
            Dictionary of energy terms with their values
        """
        energy_dic: dict[str, float] = {}
        for line in remarks:
            # TODO: use regex to do this
            if 'energies' in line:
                energy_values = map(
                    float,
                    line.rstrip().split(':')[-1].split(',')
                    )
                total, bonds, angles, improper, dihe, vdw, elec, air, cdih, coup, rdcs, vean, dani, xpcs, rg = energy_values  # noqa: E501
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
            if 'Symmetry energy' in line:
                sym = float(line.rstrip().split(':')[-1])
                energy_dic['sym'] = sym

        return energy_dic

    @staticmethod
    def _load_per_interface_energies(
            remarks: list[str],
            ) -> dict[str, dict[str, float]]:
        """Read a pdb file and parse per interface scores.

        Parameters
        ----------
        remarks : list[str]
            List of PDB lines containing remarks.

        Returns
        -------
        interfaces_scores : dict[str, dict[str, float]]
            Dictionary holding per interfaces scores.
        """
        header = None
        interfaces_scores: dict[str, dict[str, float]] = {}
        for _ in remarks:
            if _.startswith("REMARK Interface"):
                s_ = _.strip().split()[2:]
                # Extract header
                if not header:
                    header = s_
                # Extract data
                else:
                    chain1 = s_[header.index("Chain1")]
                    chain2 = s_[header.index("Chain2")]
                    haddockscore = float(s_[header.index("HADDOCKscore")])
                    evdw = float(s_[header.index("Evdw")])
                    eelec = float(s_[header.index("Eelec")])
                    edesol = float(s_[header.index("Edesol")])
                    bsa = float(s_[header.index("BSA")])
                    # Combine chains together
                    chains_key = f"{chain1}_{chain2}"
                    # Hold data
                    interfaces_scores[chains_key] = {
                        "HADDOCKscore": haddockscore,
                        "vdw": evdw,
                        "elec": eelec,
                        "desolv": edesol,
                        "bsa": bsa,
                        }
        return interfaces_scores

    def calc_haddock_score(self, **weights: float) -> float:
        """Calculate the haddock score based on the weights and energies."""
        return self.calc_score(self.energies, **weights)

    @staticmethod
    def calc_score(energies: dict[str, float], **weights: float) -> float:
        """Compute sum of weighted energy terms.

        Parameters
        ----------
        energies : dict[str, float]
            Dict of energy values for each energy term

        Returns
        -------
        weighted_score : float
            Sum of weighted energy terms
        """
        weighted_score: float = 0.0
        for weight_name, weight_value in weights.items():
            component_id = weight_name.split("_")[1]
            try:
                weighted_score += energies[component_id] * weight_value
            except KeyError:
                continue
        return weighted_score
