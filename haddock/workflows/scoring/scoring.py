# convert the analysis scripts to python
import os
import subprocess
import itertools
import glob
import calc_fcc_matrix
import cluster_fcc


def load_structure(pdb_f):
    prot_dic = {}
    for line in open(pdb_f):
        if line.startswith(('ATOM', 'HETATM', 'ANISOU')):
            if 'HN' in line:
                new_line = line[:12] + ' H  ' + line[16:]
                line = new_line
            segid = line[72:76].strip()
            if segid in prot_dic:
                prot_dic[segid].append(line)
            else:
                prot_dic[segid] = [line]
    return prot_dic


def run_fastcontact(pdb_f):
    fastcontact_exe = '/Users/rodrigo/software/fastcontact/fastcontact'
    os.environ["FASTCONTACTDIR"] = "/Users/rodrigo/software/fastcontact/"
    felecv = .0
    fdesolv = .0
    freee = .0

    prot_dic = load_structure(pdb_f)

    for segid_x, segid_y in itertools.combinations(prot_dic, 2):

        prot_x = '{}.pdb'.format(segid_x)
        prot_y = '{}.pdb'.format(segid_y)
        open(prot_x,'w').write(''.join(prot_dic[segid_x]))
        open(prot_y, 'w').write(''.join(prot_dic[segid_y]))

        p = subprocess.Popen([fastcontact_exe, prot_x, prot_y], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        out = p.communicate()
        i, j, k = map(float, str(out[0]).split('\\n')[1].split())
        felecv += i
        fdesolv += j
        freee += k

        os.remove(prot_x)
        os.remove(prot_y)

    return felecv, fdesolv


def identify_chains(pdb):
    chain_l = []
    for l in open(pdb):
        if 'ATOM' in l[:4]:
            chain = l[21]
            if chain not in chain_l:
                chain_l.append(chain)
    return chain_l


# def run_dfire(pdb_f):
#     # TODO: properly compile and execute
#     dcomplex_exe = '/Users/rodrigo/software/dcomplex_single_file/dcomplex'
#     chain_l = identify_chains(pdb_f)
#
#     for chain_x, chain_y in itertools.combinations(chain_l, 2):
#         p = subprocess.Popen([dcomplex_exe, pdb_f, chain_x, chain_y], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
#         out = p.communicate()


# def run_molprob(pdb_f):
#     # TODO: install software and test
#     pass


def extract_energies(pdb_f):
    # extract information present on pdb file
    pdb_l = open(pdb_f).readlines()
    for i, line in enumerate(pdb_l):
        if 'total,bonds,angles,improper,dihe,vdw,elec,noe,cdih,coup,sani,vean,dani' in line:
            energies_l = map(float, pdb_l[i+1].split(':')[-1].split(','))
            total,bonds,angles,improper,dihe,vdw,elec,noe,cdih,coup,sani,vean,dani = energies_l
        if 'buried surface area' in line:
            bsa = float(line.split(':')[-1])
        if 'Desolvation energy' in line:
            desolv = float(line.split(':')[-1])
        if 'Binding energy' in line:
            binding = float(line.split(':')[-1])
        if 'violations.' in line:
            violation_l = map(float, line.split(':')[-1].split(','))
            noe_v, cdih_v, coup_v, sani_v, vean_v, dani_v = violation_l
    return (total,bonds,angles,improper,dihe,vdw,elec,noe,cdih,coup,sani,vean,dani,bsa,desolv,binding,noe_v,cdih_v,coup_v,sani_v,vean_v,dani_v)


def extract_haddock_score(pdb_f):
    hs = None
    if not os.path.isfile('file.list'):
        print('file.list not found')
        hs = float('nan')
    else:
        for line in open('file.list'):
            if pdb_f in line:
                hs = float(line.split('{')[-1][:-2])

    return hs


def calculate_rmsd(struc_dic):
    # TODO: implement multichain support
    profit_exe = '/programs/i386-mac/profit/3.1/profit'
    rmsd_d = {}

    # sort by haddock score
    hs_list = [(k, struc_dic[k]['haddock-score']) for k in struc_dic]
    sorted_hs_list = sorted(hs_list, key=lambda x: (-x[1], x[0]))
    sorted_hs_list.reverse()
    sorted_pdb_list = [e[0] for e in sorted_hs_list]

    refe = sorted_pdb_list[0]
    for mobi in sorted_pdb_list:
        cmd = 'refe {}\nmobi {}\nignore missing\natoms CA\nzone A*\nfit\nrzone B*\nquit'.format(refe, mobi)
        output = os.popen('echo "{}" | {}'.format(cmd, profit_exe))
        result = [l for l in output if 'RMS:' in l][-1]
        rmsd = float(result.split()[-1])
        rmsd_d[mobi] = rmsd

    return rmsd_d


def analyze_structures(pdb_list):

    structure_dic = {}
    for pdb in pdb_list:
        structure_dic[pdb] = {}

        fast_elec, fast_desol = run_fastcontact(pdb)
        # haddock_score = extract_haddock_score(pdb)
        total, bonds, angles, improper, dihe, vdw,  elec, noe, cdih, coup, sani, vean, dani, bsa, desolv, binding, noe_v, \
            cdih_v, coup_v, sani_v, vean_v, dani_v = extract_energies(pdb)
        # run_dfire(pdb)
        # run_molprob(pdb)
        # probe_score = float('nan')
        # probe_score2 = float('nan')

        structure_dic[pdb]['haddock-score'] = (1.0 * vdw) + (0.2 * elec) + (1.0 * desolv)
        structure_dic[pdb]['total'] = total
        structure_dic[pdb]['bonds'] = bonds
        structure_dic[pdb]['angles'] = angles
        structure_dic[pdb]['improper'] = improper
        structure_dic[pdb]['dihe'] = dihe
        structure_dic[pdb]['vdw'] = vdw
        structure_dic[pdb]['elec'] = elec
        # structure_dic[pdb]['air'] = air
        structure_dic[pdb]['cdih'] = cdih
        structure_dic[pdb]['coup'] = coup
        # structure_dic[pdb]['rdcs'] = rdcs
        structure_dic[pdb]['vean'] = vean
        structure_dic[pdb]['dani'] = dani
        # structure_dic[pdb]['xpcs'] = xpcs
        # structure_dic[pdb]['rg'] = rg
        structure_dic[pdb]['bsa'] = bsa
        structure_dic[pdb]['desolv'] = desolv
        structure_dic[pdb]['binding'] = binding
        structure_dic[pdb]['noe_v'] = noe_v
        structure_dic[pdb]['cdih_v'] = cdih_v
        structure_dic[pdb]['coup_v'] = coup_v
        structure_dic[pdb]['sani_v'] = sani_v
        structure_dic[pdb]['vean_v'] = vean_v
        structure_dic[pdb]['dani_v'] = dani_v
        structure_dic[pdb]['fastelec'] = fast_elec
        structure_dic[pdb]['fastdesol'] = fast_desol

    rmsd_dic = calculate_rmsd(structure_dic)
    for k in rmsd_dic:
        structure_dic[k]['rmsd'] = rmsd_dic[k]

    hs_list = [(k, structure_dic[k]['haddock-score']) for k in structure_dic]
    sorted_hs_list = sorted(hs_list, key=lambda x: (-x[1], x[0]))
    sorted_hs_list.reverse()
    sorted_pdb_list = [e[0] for e in sorted_hs_list]

    header = '#struc\trmsd_all\thaddock-score\tbsa\tEdesolv\tFastcontact-elec\tFastcontact-desolv\t' \
             'Dfire-Ebinding[kcal/mol]\tprobe-score\tprobe-score/A**2\tcombined-score'

    tbw = ''
    for pdb in sorted_pdb_list:
        tbw += f"{pdb}\t{structure_dic[pdb]['rmsd']:9.3f}\t{structure_dic[pdb]['haddock-score']:9.3f}\t" \
            f"{structure_dic[pdb]['bsa']:9.3f}\t{structure_dic[pdb]['fastelec']:9.3f}\t" \
            f"{structure_dic[pdb]['fastdesol']:9.3f}\t{float('nan')}\t{float('nan')}\t{float('nan')}\t{float('nan')}\n"

    out = open('structures.ranking','w')
    out.write(header)
    out.write(tbw)
    out.close()

    return structure_dic


def run_contacts(pdb_l, d_cutoff=5.0):
    contacts_exe = '/Users/rodrigo/software/contact_fcc'
    contact_l = []
    for pdb in pdb_l:

        outfile = os.path.join(os.path.dirname(pdb), "{}.contacts".format(pdb[:-4]))

        p = subprocess.Popen([contacts_exe, pdb, str(d_cutoff)], stdout=subprocess.PIPE)

        p_output = p.communicate()[0]
        contacts = set(map(int, sorted(list(set([l for l in p_output.decode('utf-8').split('\n')][:-1])))))
        contact_l.append(contacts)
        # set([int(l) for l in contacts])
        # f_contact = open(outfile, 'w')
        # f_contact.write('\n'.join(contacts))
        # f_contact.close()
    return contact_l


def calculate_matrix(pdb_l):

    fcc_matrix_name = 'target.fcc_matrix'
    c = run_contacts(pdb_l)
    m = calc_fcc_matrix.calculate_pairwise_matrix(c, True)
    f = open(fcc_matrix_name, 'w')
    calc_fcc_matrix._output_fcc(f.write, m, 50000) # is this use correct?
    f.close()

    return fcc_matrix_name


def cluster(matrix_f):

    cluster_size_l = [1,2,4]
    cutoff_l = [0.75, 0.65]
    for size in cluster_size_l:
        for cutoff in cutoff_l:
            pool = cluster_fcc.read_matrix(matrix_f, cutoff, 0.75)
            element_pool, clusters = cluster_fcc.cluster_elements(pool, size)

            out = open(f'cluster{cutoff}-{size}.out','w')
            cluster_fcc.output_clusters(out, clusters)
            out.close()


def main():

    pdb_l = glob.glob('*conv.pdb')
    structure_dic = analyze_structures(pdb_l)

    # calculate fcc matrix
    matrix_f = calculate_matrix(pdb_l)
    cluster(matrix_f)



if __name__ == '__main__':
    main()
