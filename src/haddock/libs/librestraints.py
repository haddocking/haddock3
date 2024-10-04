from pathlib import Path
import logging
import itertools
import re
import random
import sys

from Bio.PDB import PDBParser, NeighborSearch
from freesasa import Classifier, structureFromBioPDB, calc

# Scaling factors for relative ASA
# Calculated using extended ALA-X-ALA peptides
# Taken from NACCESS
# Taken from NACCESS
REL_ASA = {
    'total':
        {
            'ALA': 107.95,
            'CYS': 134.28,
            'ASP': 140.39,
            'GLU': 172.25,
            'PHE': 199.48,
            'GLY': 80.10,
            'HIS': 182.88,
            'ILE': 175.12,
            'LYS': 200.81,
            'LEU': 178.63,
            'MET': 194.15,
            'ASN': 143.94,
            'PRO': 136.13,
            'GLN': 178.50,
            'ARG': 238.76,
            'SER': 116.50,
            'THR': 139.27,
            'VAL': 151.44,
            'TRP': 249.36,
            'TYR': 212.76,
            'ASH': 140.39,
            'DDZ': 107.95,
            'GLH': 172.25,
            'CYM': 134.28,
            'CSP': 134.28,
            'CYF': 134.28,
            'CYC': 134.28,
            'CFE': 134.28,
            'NEP': 182.88,
            'ALY': 200.81,
            'MLZ': 200.81,
            'MLY': 200.81,
            'M3L': 200.81,
            'HYP': 136.13,
            'SEP': 116.50,
            'TOP': 139.27,
            'TYP': 212.76,
            'PTR': 212.76,
            'TYS': 212.76,
            'PNS': 116.50,
            },
    'bb':
        {
            'ALA': 38.54,
            'CYS': 37.53,
            'ASP': 37.70,
            'GLU': 37.51,
            'PHE': 35.37,
            'GLY': 47.77,
            'HIS': 35.80,
            'ILE': 37.16,
            'LYS': 37.51,
            'LEU': 37.51,
            'MET': 37.51,
            'ASN': 37.70,
            'PRO': 16.23,
            'GLN': 37.51,
            'ARG': 37.51,
            'SER': 38.40,
            'THR': 37.57,
            'VAL': 37.16,
            'TRP': 38.10,
            'TYR': 35.38,
            'ASH': 37.70,
            'DDZ': 38.54,
            'GLH': 37.51,
            'CYM': 37.53,
            'CYC': 37.53,
            'CSP': 37.53,
            'CYF': 37.53,
            'CFE': 37.53,
            'NEP': 35.80,
            'ALY': 37.51,
            'MLZ': 37.51,
            'MLY': 37.51,
            'M3L': 37.51,
            'HYP': 16.23,
            'SEP': 38.40,
            'TOP': 37.57,
            'TYP': 35.38,
            'PTR': 35.38,
            'TYS': 35.38,
            'PNS': 38.40,
            },
    'sc':
        {
            'ALA': 69.41,
            'CYS': 96.75,
            'ASP': 102.69,
            'GLU': 134.74,
            'PHE': 164.11,
            'GLY': 32.33,
            'HIS': 147.08,
            'ILE': 137.96,
            'LYS': 163.30,
            'LEU': 141.12,
            'MET': 156.64,
            'ASN': 106.24,
            'PRO': 119.90,
            'GLN': 140.99,
            'ARG': 201.25,
            'SER': 78.11,
            'THR': 101.70,
            'VAL': 114.28,
            'TRP': 211.26,
            'TYR': 177.38,
            'ASH': 102.69,
            'DDZ': 69.41,
            'GLH': 134.74,
            'CYM': 96.75,
            'CYC': 96.75,
            'CSP': 96.75,
            'CYF': 96.75,
            'CFE': 96.75,
            'NEP': 147.08,
            'ALY': 163.30,
            'MLZ': 163.30,
            'MLY': 163.30,
            'M3L': 163.30,
            'HYP': 119.90,
            'SEP': 78.11,
            'TOP': 101.70,
            'TYP': 177.38,
            'PTR': 177.38,
            'TYS': 177.38,
            'PNS': 78.11,
            }
    }
DEFAULT_PROBE_RADIUS = 1.4

def get_surface_resids(structure, cutoff=15):
    """
    Calls freesasa using its Python API and returns
    per-residue accessibilities.
    """
    asa_data, rsa_data, rel_main_chain, rel_side_chain = {}, {}, {}, {}
    _rsa = REL_ASA['total']
    _rsa_bb = REL_ASA['bb']
    _rsa_sc = REL_ASA['sc']

    #classifier = Classifier(config_path)
    classifier = Classifier()

    struct = structureFromBioPDB(structure, classifier, )
    result = calc(struct)

    # iterate over all atoms to get SASA and residue name
    for idx in range(struct.nAtoms()):
        atname = struct.atomName(idx).strip()
        resname = struct.residueName(idx)
        resid = int(struct.residueNumber(idx))
        chain = struct.chainLabel(idx)
        at_uid = (chain, resname, resid, atname)
        res_uid = (chain, resname, resid)

        asa = result.atomArea(idx)
        asa_data[at_uid] = asa
        # add asa to residue
        rsa_data[res_uid] = rsa_data.get(res_uid, 0) + asa

        if atname in ('C', 'N', 'O'):
            rel_main_chain[res_uid] = rel_main_chain.get(res_uid, 0) + asa
        else:
            rel_side_chain[res_uid] = rel_side_chain.get(res_uid, 0) + asa

    # convert total asa ro relative asa
    rsa_data.update((res_uid, asa / _rsa[res_uid[1]]) for res_uid, asa in rsa_data.items())
    rel_main_chain.update((res_uid, asa / _rsa_bb[res_uid[1]] * 100) for res_uid, asa in rel_main_chain.items())
    rel_side_chain.update((res_uid, asa / _rsa_sc[res_uid[1]] * 100) for res_uid, asa in rel_side_chain.items())

    # We format to fit the pipeline
    resid_access = {}
    for res_uid, access in rel_main_chain.items():
        resid_access[res_uid[2]] = {'side_chain_rel': rel_side_chain.get(res_uid), 'main_chain_rel': access}
    surface_resids = [r for r, v in resid_access.items() if v['side_chain_rel'] >= cutoff or
                      v['main_chain_rel'] >= cutoff]
    return surface_resids


def parse_actpass_file(actpass_file):
    """Parse actpass file
    
    Parameters
    ----------
    actpass_file: str or Path
        path to actpass_file

    Returns
    -------
    active: list
        list of active residues
    passive: list
        list of passive residues
    """
    
    if Path(actpass_file).exists() is False:
        raise Exception(f"actpass file {actpass_file} does not exist.")

    lines = open(actpass_file, "r").readlines()
    nlines = len(lines)
    if nlines != 2:
        raise Exception(f"actpass file {actpass_file} does not have two lines (counted {nlines})")
    active, passive = [[int(x) for x in line.split()] for line in lines]
    return active, passive


def active_passive_to_ambig(active1, passive1, active2, passive2, segid1='A', segid2='B'):
    """Convert active and passive residues to Ambiguous Interaction Restraints

    Parameters
    ----------
    active1 : list
        List of active residue numbers of the first segid

    passive1 : list
        List of passive residue numbers of the first segid

    passive2 : list
        List of passive residue numbers of the second segid

    active2 : list
        List of active residue numbers of the second segid
    
    active2 : list
        List of passive residue numbers of the second segid

    segid1 : string
        Segid to use for the first model

    segid2 : string
        Segid to use for the second model

    """

    all1 = active1 + passive1
    all2 = active2 + passive2

    for resi1 in active1:
        print('assign (resi {:d} and segid {:s})'.format(resi1, segid1))
        print('(')
        c = 0
        for resi2 in all2:
            print('       (resi {:d} and segid {:s})'.format(resi2, segid2))
            c += 1
            if c != len(all2):
                print('        or')

        print(') 2.0 2.0 0.0\n')
            
    for resi2 in active2:
        print('assign (resi {:d} and segid {:s})'.format(resi2, segid2))
        print('(')
        c = 0
        for resi1 in all1:
            print('       (resi {:d} and segid {:s})'.format(resi1, segid1))
            c += 1
            if c != len(all1):
                print('        or')

        print(') 2.0 2.0 0.0\n')


# Functions/Methods

def calc_euclidean(i, j):
	return ((j[0]-i[0])**2 + (j[1]-i[1])**2 + (j[2]-i[2])**2)**0.5


def read_structure(pdbf, exclude=None):
	"""
	Reads a PDB file and returns a list of parsed atoms
	"""
	_atoms = {'CA', 'P'}  # Alpha-Carbon (Prot), Backbone Phosphorous (DNA)
	_altloc = {' ', 'A'}

	if not exclude:
		exclude = set()
	else:
		exclude = set(exclude)

	res_list = []
	with open(pdbf, 'r') as pdb_handle:
		for line in pdb_handle:
			field = line[0:4]
			if field != 'ATOM':
				continue

			aname = line[12:16].strip()
			chain = line[21] if line[21].strip() else line[72:76].strip()  # chain ID or segID
			if chain not in exclude and aname in _atoms and line[16] in _altloc:
				resi = int(line[22:26])
				coords = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
				res_list.append((chain, resi, aname, coords))

	if not res_list:
		logging.critical('[!] PDB File seems empty or no CA/P atoms found: {0}'.format(pdbf))
		sys.exit(1)

	return res_list


def get_bodies(atom_lst, prot_threshold=4.0, dna_threshold=7.5):
	"""
	Determines gaps in an atom list following simple distance based criteria.
	Returns continuous fragments.
	"""

	bodies = []

	threshold = {'CA': prot_threshold, 'P': dna_threshold}

	body_start = 0
	i = None
	for i, atom in enumerate(atom_lst[1:], start=1):
		p_atom = atom_lst[i-1]

		chain, resi, aname, xyz = atom
		p_chain, p_resi, p_aname, p_xyz = p_atom

		if (chain == p_chain) and (aname == p_aname):  # Internal Gap
			d_xyz = calc_euclidean(xyz, p_xyz)
			if d_xyz >= threshold[aname]:
				logging.debug('[+++] (Internal) Body: {0}:{1}'.format(body_start, i-1))
				bodies.append((body_start, i-1))
				body_start = i  # Set new beginning

		elif (chain != p_chain) or (aname != p_aname):  # Different molecules/types
			logging.debug('[+++] Body: {0}:{1}'.format(body_start, i-1))
			bodies.append((body_start, i-1))
			body_start = i  # Set new beginning

	if not bodies:  # Single continuous molecule
		bodies.append((0, len(atom_lst)))
	else:
		logging.debug('[+++] Body: {0}:{1}'.format(body_start, i))
		bodies.append((body_start, i))  # Last body

	logging.info('[++] Found {0} bodies'.format(len(bodies)))

	return bodies


def build_restraints(bodies):
	"""
	Generates distance restraints to maintain the relative
	orientation of the different bodies during the simulations.

	Generates two unique restraints per pair of bodies.

	Each restraint is created using two random atoms on each body
	and using their exact euclidean distance as target distance.
	"""

	def pick_residues(body, max_trials=10):
		# Pick two random residues in each body
		# Make sure they are far apart from each other
		n_trials = 0
		while 1:
			try:
				res_i, res_ii = random.sample(body, 2)
			except ValueError:
				# Likely, sample size is 1
				logging.warning('[!] One-sized body found. This may lead to problems..')
				return body[0], body[0]

			logging.debug('[+++] Trial {0}: {1} & {2}'.format(n_trials, res_i, res_ii))
			if abs(res_i - res_ii) > 3:
				logging.info('[++] Picked residues {0} & {1}'.format(res_i, res_ii))
				return res_i, res_ii
			n_trials += 1
			if n_trials == max_trials:
				msg = '[!] Could not pick two unique distant residues in body after {0} tries'
				logging.info(msg.format(max_trials))
				return res_i, res_ii

	restraints = []

	n_bodies = range(len(bodies))
	combinations = itertools.combinations(n_bodies, 2)

	for pair_bodies in combinations:
		body_i, body_j = pair_bodies
		logging.debug('[+++] Restraining body {0} to body {1}'.format(body_i, body_j))

		st_body_i, en_body_i = bodies[body_i]
		st_body_j, en_body_j = bodies[body_j]
		res_i, res_ii = pick_residues(range(st_body_i, en_body_i+1))
		res_j, res_jj = pick_residues(range(st_body_j, en_body_j+1))

		logging.info('[++] Created restraint: {0}:{1} <--> {2}:{3}'.format(body_i, res_i, body_j, res_j))
		restraints.append((res_i, res_j))
		logging.info('[++] Created restraint: {0}:{1} <--> {2}:{3}'.format(body_i, res_ii, body_j, res_jj))
		restraints.append((res_ii, res_jj))

	return restraints


def generate_tbl(atom_lst, restraints):
	"""
	Makes a list of TBL-formatted restraints.

    Parameters
    ----------
    atom_lst : list
        List of atoms in the form (chain, resi, aname, coords)
    
    restraints : list
        List of restraints in the form (res_i, res_j)
	"""

	for r in restraints:
		i, j = r
		atom_i, atom_j = atom_lst[i], atom_lst[j]
		dist_ij = calc_euclidean(atom_i[3], atom_j[3])

		tbl = "assign (segid {0[0]} and resi {0[1]} and name {0[2]}) ".format(atom_i)
		tbl += "(segid {0[0]} and resi {0[1]} and name {0[2]}) ".format(atom_j)
		tbl += "{0:3.3f} 0.0 0.0".format(dist_ij)
		print(tbl)


def check_parenthesis(file):
    open_parenthesis = re.compile('[(]')
    close_parenthesis = re.compile('[)]')
    quotation_marks = re.compile('[\"]')
    opened = 0
    closed = 0
    quote = 0
    for _match in open_parenthesis.finditer(file):
        opened += 1
    for _match in close_parenthesis.finditer(file):
        closed += 1
    for _match in quotation_marks.finditer(file):
        quote += 1
    if opened != closed:
        raise Exception("Problem with TBL file parentheses ({:d} opening for {:d} "
                        "closing parentheses)".format(opened, closed))
    if quote % 2 != 0:
        raise Exception("Problem with TBL file, odd number of quotation marks "
                        "({:d} quotation marks)".format(quote))


def validate_tbldata(restraints, pcs=False):
    # List of selection keywords
    selectors = ('name', 'resn', 'atom', 'resi', 'attr', 'segi', 'chem', 'id',
                 'byres', 'not')
    connectors = ('and', 'or', 'not', 'byres')
    output = ""
    parentmatch = re.compile('[()]')
    # Global mode is activated outside any assign statement
    mode = "global"
    # Remove any carriage return/new line caracters
    lines = restraints.replace('\r', '').split("\n")
    # Line number
    lnr = 0
    # Temporary line storage for future output
    tmp_output = None
    tmp = None
    level = None
    s = None
    selections = None
    postselections = None
    numbers = None
    types = None
    lastassign = None
    for line in lines:
        lnr += 1
        # Take everything which is before putative "!" (comment) caracter
        if line.find("!") > -1:
            line = line[:line.find("!")]
        # Remove whitespaces and merge with previous line if in OR statement
        # for AIR restraints
        if mode != "format":
            line = line.strip()
        else:
            line = tmp + line.strip()
            mode = "postassign"
        # End of line
        if not len(line):
            continue
        # Check if all "" are closed
        if line.count('"') % 2 > 0:
            raise Exception('Unclosed " at line {:d} for line: {}'.format(lnr, line))
        if mode in ("global", "postglobal"):
            # Start of an assign statement
            if line.lower().startswith("assi"):
                if mode == "postglobal":
                    output += "\n!\n"
                    output += tmp_output
                mode = "assign"
                selections = []
                # assign is the only character of the line,check following lines
                if line.find("(") == -1:
                    continue
                # Reset temporary buffer
                tmp_output = None
            # Check for "OR" restraint
            elif mode == "postglobal" and line.lower().startswith("or"):
                mode = "postassign"
                selections = []
                line = line[len("or"):]
                # Case where "OR" statement is the only one on the line
                # (rest of selection at the next line)
                if line == "":
                    continue
            # We are not treating an assign params (postglobal) or at the
            # end (global) and no "OR" restraint is present -> ERROR
            else:
                raise Exception("Invalid TBL file: Unknown statement (line {:d}): {}".format(lnr, line))
        # We check if the selection is made over two lines
        # (thanks to "segid" keyword)
        if mode == "postassign":
            if line.count("segid") == 1:
                mode = "format"
                tmp = line
                continue

        matched = True
        while matched:
            matched = False
            # We are looking for parenthesis as start of the selections
            if mode in ("assign", "postassign"):
                # Ambiguous restraint for two different pairs of atoms
                # (ex: THR 1 B <-> ALA 2 A OR GLY 2 B <-> ASP 10 A )
                pos = line.find("(")
                if pos != -1:
                    matched = True
                    line = line[pos + 1:]
                    # We look for an "OR" selection
                    if mode == "postassign":
                        mode = "postsel"
                    else:
                        mode = "sel"
                    lastassign = lnr
                    s = ""
                    level = 1
            # Get the structural selections
            if mode in ("sel", "postsel"):
                # Detect opening and closing parenthesis
                for match in parentmatch.finditer(line):
                    if match.group() == "(":
                        level += 1
                    else:
                        level -= 1
                    # End of the parenthesis content
                    if level == 0:
                        if mode == "postsel":
                            mode = "postassign"
                        else:
                            mode = "assign"
                        matched = True
                        # print l[:match.start()]
                        # Get back the parentheses content and check
                        # for selection keywords
                        idx_connectors = []
                        sel = line[:match.start()]
                        # We can directly check for the 1st keyword presence
                        if len(sel) > 0:
                            syntax_ok = False
                            for s in selectors:
                                if sel.strip().startswith(s):
                                    syntax_ok = True
                            if not syntax_ok and not \
                                    sel.strip().startswith('byres') and not \
                                    sel.strip().startswith('not'):
                                raise Exception("1) Missing or wrong keyword in-term {} (stopped at line {:d})".
                                                format(sel, lnr))
                        # Get indexes of connectors (AND, OR, etc.)
                        for c in connectors:
                            idx_connectors.append([m.start()+len(c) for m in
                                                   re.finditer(c, sel)])
                        # Flatten the list
                        tmp_list = [item for sublist in idx_connectors for
                                    item in sublist]
                        idx_connectors = tmp_list
                        # Check that each connector is followed by a keyword
                        for i in idx_connectors:
                            new_sel = sel[i:].strip().strip("(").strip()
                            if len(new_sel) > 0:
                                syntax_ok = False
                                for s in selectors:
                                    if new_sel.startswith(s):
                                        syntax_ok = True
                                        break
                                if not syntax_ok:
                                    raise Exception("Missing or wrong keyword in-term {} (stopped at line {:d})".
                                                    format(sel, lnr))
                        s += sel
                        selections.append(s)
                        s = None
                        # Get the rest of the line
                        line = line[match.end():]
                        # Go to the process of the selection
                        break
                # No parenthesis, we get the whole line
                if level > 0:
                    if line != "":
                        s += line + "\n\t"
                    else:
                        # To avoid multiple blank lines
                        if repr(s) == "''":
                            # print repr(l), repr(s), selections
                            s += "\n\t"
        # TBD
        if mode in ("sel", "postsel"):
            continue
        # Selection parsed for "OR" line
        if mode == "postassign":
            if len(selections) != postselections:
                raise Exception("Invalid TBL file: wrong number of selections: in-term {:d}, cross-term {:d} "
                                "(stopped at line {:d})".format(postselections, len(selections), lnr))
            tmp_output += " or"
            for s in selections:
                tmp_output += "\t({})\n".format(s)
            # We let the possibility for other "OR"
            mode = "postglobal"
        if len(line) == 0:
            continue
        # Define the distance restraints type to adapt the parsing
        if mode == "assign":
            mode = "numbers"
            if pcs:
                if len(selections) == 5:
                    types = (" {:.3f}", " {:.3f}")
                else:
                    raise Exception("Invalid TBL file: wrong number of selections (must be 5 in PCS mode)")
            else:
                if len(selections) == 2:
                    types = (" {:.3f}", " {:.3f}", " {:.3f}")
                elif len(selections) == 4:
                    types = (" {:.3f}", " {:.3f}", " {:.3f}", " {:d}")
                elif len(selections) == 5:
                    raise Exception("Invalid TBL file: wrong number of selections (can be 5 only in PCS mode)")
                elif len(selections) == 6:
                    types = (" {:.3f}", " {:.3f}")
                else:
                    check_parenthesis(restraints)
                    raise Exception("Invalid TBL file: wrong number of selections (must be 2, 4 or 6)")
            postselections = len(selections)
            numbers = []
        # Distance restraints parsing
        if mode == "numbers":
            ll = line.split()
            for num in ll:
                if len(numbers) == len(types):
                    break
                numbers.append(float(num))
            if len(numbers) == len(types):
                tmp_output = "assign "
                for s in selections:
                    tmp_output += "\t({})\n".format(s)
                tmp_output = tmp_output[:-len("\n")]
                for n, t in zip(numbers, types):
                    tmp_output += t.format(n)
                tmp_output += "\n"
                mode = "postglobal"
    # "OR" lines have been parsed, we store the selections
    if mode == "postglobal":
        output += "!\n"
        output += tmp_output
        mode = "global"
    # If mode is not back to global, something has not been processed properly
    if mode != "global":
        raise Exception("Invalid TBL file: Malformed ASSIGN statement (line {:d}), use --quick to check for "
                        "putative syntax issues".format(lastassign))
    if not len(output.strip()):
        raise Exception("Invalid or empty TBL file")

    # Remove extra lines before each assign statement
    output = output.replace("\n\n", "\n")
    # Remove extra line at the beginning of the file
    if output.startswith("\n"):
        output = output.replace("\n", "", 1)
    return output

def passive_from_active_raw(structure, active, chain_id=None, surface=None, radius=6.5):
    """Get the passive residues.
    
    Parameters
    ----------

    structure : str
        path to the PDB file

    active : list
        List of active residues

    chain_id : str
        Chain ID
    
    surface : list
       List of surface residues.

    radius : float
        Radius from active residues
    """

    # Parse the PDB file
    if Path(structure).exists():
        p = PDBParser(QUIET=True)
        s = p.get_structure('pdb', structure)
    else:
        raise FileNotFoundError('File not found: {0}'.format(structure))

    try:
        if chain_id:
            atom_list = [a for a in s[0][chain_id].get_atoms()]
        else:
            atom_list = [a for a in s[0].get_atoms()]
    except KeyError as e:
        raise KeyError('Chain {0} does not exist in the PDB file {1}, please enter a proper chain id'.
              format(chain_id, structure)) from e

    act_atoms = [a.get_coord() for a in atom_list if a.parent.id[1] in active]

    try:
        if not surface:
            surface = get_surface_resids(s)
    except Exception as e:
        raise Exception("There was an error while calculating surface residues: {}".format(e)) from e

    ns = NeighborSearch(atom_list)
    neighbors = []
    for a in act_atoms:
        neighbors.append(ns.search(a, radius, "R"))  # HADDOCK used 6.5A as default

    passive_list = set()
    for n in neighbors:
        for r in n:
            passive_list.add(r.id[1])
    tmp = passive_list & set(surface)
    passive_list = tmp - set(active)
    return sorted(passive_list)
