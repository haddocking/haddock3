from pathlib import Path

ATOMS_TO_BE_MUTATED = ['C', 'N', 'CA', 'O', 'CB']
RES_CODES = dict([
            ("CYS", "C"),
            ("ASP", "D"),
            ("SER", "S"),
            ("GLN", "Q"),
            ("LYS", "K"),
            ("ILE", "I"),
            ("PRO", "P"),
            ("THR", "T"),
            ("PHE", "F"),
            ("ASN", "N"),
            ("GLY", "G"),
            ("HIS", "H"),
            ("LEU", "L"),
            ("ARG", "R"),
            ("TRP", "W"),
            ("ALA", "A"),
            ("VAL", "V"),
            ("GLU", "E"),
            ("TYR", "Y"),
            ("MET", "M"),
        ])


def mutate(pdb_f, target_chain, target_resnum, mut_resname):
    """Mutate a resnum to a resname."""
    mut_pdb_l = []
    resname = ''
    with open(pdb_f, 'r') as fh:
        for line in fh.readlines():
            if line.startswith('ATOM'):
                chain = line[21]
                resnum = int(line[22:26])
                atom_name = line[12:16].strip()
                if target_chain == chain and target_resnum == resnum:
                    if not resname:
                        resname = line[17:20].strip()
                    if atom_name in ATOMS_TO_BE_MUTATED:
                        # mutate
                        line = line[:17] + mut_resname + line[20:]
                        mut_pdb_l.append(line)
                else:
                    mut_pdb_l.append(line)
    mut_id = f'{RES_CODES[resname]}{target_resnum}{RES_CODES[mut_resname]}'
    mut_pdb_fname = Path(
        pdb_f.name.replace('.pdb', f'-{target_chain}_{mut_id}.pdb'))
    with open(mut_pdb_fname, 'w') as fh:
        fh.write(''.join(mut_pdb_l))
    return mut_pdb_fname


def add_delta_to_bfactor(pdb_f, bfactor_dic):
    """Add delta scores as b-factors."""
    output_pdb_f = pdb_f.replace('.pdb', '_bfactor.pdb')
    out_pdb_l = []
    with open(pdb_f, 'r') as fh:
        for line in fh.readlines():
            if line.startswith('ATOM'):
                chain = line[21]
                resnum = int(line[22:26])
                # bfactor = float(line[60:66])
                delta = 0.0
                ident = f'{chain}.{resnum}'
                if ident in bfactor_dic:
                    delta = bfactor_dic[ident]

                delta_str = f"{delta:.2f}".rjust(6, " ")
                line = line[:60] + delta_str + line[66:]
            out_pdb_l.append(line)
    with open(output_pdb_f, 'w') as out_fh:
        out_fh.write(''.join(out_pdb_l))

    return output_pdb_f
