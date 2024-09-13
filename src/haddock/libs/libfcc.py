"""FCC related functions

NOTE: This functions were ported directly from `https://github.com/haddocking/fcc`!
"""


class Element:
    """Defines a 'clusterable' Element"""

    __slots__ = ["name", "cluster", "neighbors"]

    def __init__(self, name):
        self.name = name
        self.cluster = 0
        self.neighbors = set()

    def add_neighbor(self, neighbor):
        """Adds another element to the neighbor list"""
        self.neighbors.add(neighbor)

    def assign_cluster(self, clust_id):
        """Assigns the Element to Cluster. 0 if unclustered"""
        self.cluster = clust_id


class Cluster:
    """Defines a Cluster. A Cluster is created with a name and a center (Element class)"""

    __slots__ = ["name", "center", "members"]

    def __init__(self, name, center):
        self.name = name
        self.center = center

        self.members = []

        self.populate()

    def __len__(self):
        return len(self.members) + 1  # +1 Center

    def populate(self):
        """
        Populates the Cluster member list through the
        neighbor list of its center.
        """

        name = self.name
        # Assign center
        ctr = self.center
        ctr.assign_cluster(name)

        mlist = self.members
        # Assign members
        ctr_nlist = (n for n in ctr.neighbors if not n.cluster)
        for e in ctr_nlist:
            mlist.append(e)
            e.assign_cluster(name)

    def add_member(self, element):
        """
        Adds one single element to the cluster.
        """
        line = self.members
        line.append(element)
        element.assign_cluster(self.name)


def cluster_elements(e_pool, threshold):
    """
    Groups Elements within a given threshold
    together in the same cluster.
    """

    cluster_list = []
    threshold -= 1  # Account for center
    ep = e_pool
    cn = 1  # Cluster Number
    while 1:
        # Clusterable elements
        ce = [e for e in ep if not ep[e].cluster]
        if not ce:  # No more elements to cluster
            break

        # Select Cluster Center
        # Element with largest neighbor list
        ctr_nlist, ctr = sorted(
            [(len([se for se in ep[e].neighbors if not se.cluster]), e) for e in ce]
        )[-1]

        # Cluster until length of remaining elements lists are above threshold
        if ctr_nlist < threshold:
            break

        # Create Cluster
        c = Cluster(cn, ep[ctr])
        cn += 1
        cluster_list.append(c)

    return ep, cluster_list


def output_clusters(handle, cluster):
    """Outputs the cluster name, center, and members."""

    write = handle.write

    for c in cluster:
        write("Cluster %s -> %s " % (c.name, c.center.name))
        for m in sorted(c.members, key=lambda k: k.name):
            write("%s " % m.name)
        write("\n")


def read_matrix(path, cutoff_param, strictness):
    """
    Reads in a four column matrix (1 2 0.123 0.456\n)
    and creates an dictionary of Elements.

    The strictness factor is a <float> that multiplies by the cutoff
    to produce a new cutoff for the second half of the matrix. Used to
    allow some variability while keeping very small interfaces from clustering
    with anything remotely similar.
    """

    cutoff_param = float(cutoff_param)
    partner_cutoff = float(cutoff_param) * float(strictness)

    elements = {}

    f = open(path, "r")
    for line in f:
        ref, mobi, d_rm, d_mr = line.split()
        ref = int(ref)
        mobi = int(mobi)
        d_rm = float(d_rm)
        d_mr = float(d_mr)

        # Create or Retrieve Elements
        if ref not in elements:
            r = Element(ref)
            elements[ref] = r
        else:
            r = elements[ref]

        if mobi not in elements:
            m = Element(mobi)
            elements[mobi] = m
        else:
            m = elements[mobi]

            # Assign neighbors
        if d_rm >= cutoff_param and d_mr >= partner_cutoff:
            r.add_neighbor(m)
        if d_mr >= cutoff_param and d_rm >= partner_cutoff:
            m.add_neighbor(r)

    f.close()

    return elements


def parse_contact_file(f_list, ignore_chain):
    """Parses a list of contact files."""

    if ignore_chain:
        contacts = [
            [int(line[0:5] + line[6:-1]) for line in open(con_f)]
            for con_f in f_list
            if con_f.strip()
        ]
    else:
        contacts = [
            set([int(line) for line in open(con_f)])
            for con_f in f_list
            if con_f.strip()
        ]

    return contacts


def calculate_fcc(list_a, list_b):
    """
    Calculates the fraction of common elements between two lists
    taking into account chain IDs
    """

    cc = len(list_a.intersection(list_b))
    cc_v = len(list_b.intersection(list_a))

    return cc, cc_v


def calculate_fcc_nc(list_a, list_b):
    """
    Calculates the fraction of common elements between two lists
    not taking into account chain IDs. Much Slower.
    """

    largest, smallest = sorted([list_a, list_b], key=len)
    ncommon = len([ele for ele in largest if ele in smallest])
    return ncommon, ncommon


def calculate_pairwise_matrix(contacts, ignore_chain):
    """Calculates a matrix of pairwise fraction of common contacts (FCC).
    Outputs numeric indexes.

    contacts: list_of_unique_pairs_of_residues [set/list]

    Returns pairwise matrix as an iterator, each entry in the form:
    FCC(cplx_1/cplx_2) FCC(cplx_2/cplx_1)
    """

    contact_lengths = []
    for con in contacts:
        try:
            ic = 1.0 / len(con)
        except ZeroDivisionError:
            ic = 0
        contact_lengths.append(ic)

    if ignore_chain:
        calc_fcc = calculate_fcc_nc
    else:
        calc_fcc = calculate_fcc

    for i in range(len(contacts)):

        for k in range(i + 1, len(contacts)):
            cc, cc_v = calc_fcc(contacts[i], contacts[k])
            fcc, fcc_v = cc * contact_lengths[i], cc * contact_lengths[k]
            yield i + 1, k + 1, fcc, fcc_v
