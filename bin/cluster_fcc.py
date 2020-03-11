#!/usr/bin/env python
# -*- coding: UTF-8  -*-

"""
Asymmetric Taylor-Butina Disjoint Clustering Algorithm.

Authors:
           RODRIGUES Joao
           TRELLET Mikael
           MELQUIOND Adrien
"""

class Element(object):
    """Defines a 'clusterable' Element"""

    __slots__ = ['name', 'cluster', 'neighbors']

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

class Cluster(object):
    """Defines a Cluster. A Cluster is created with a name and a center (Element class)"""

    __slots__ = ['name', 'center', 'members']

    def __init__(self, name, center):
        
        self.name = name
        self.center = center

        self.members = []
        
        self.populate()
        
    def __len__(self):
        return len(self.members)+1 # +1 Center
    
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
        l = self.members
        l.append(element)
        element.assign_cluster(self.name)

def read_matrix(path, cutoff, strictness):
    """ 
    Reads in a four column matrix (1 2 0.123 0.456\n) 
    and creates an dictionary of Elements.
    
    The strictness factor is a <float> that multiplies by the cutoff 
    to produce a new cutoff for the second half of the matrix. Used to
    allow some variability while keeping very small interfaces from clustering
    with anything remotely similar.
    """

    cutoff = float(cutoff)
    partner_cutoff = float(cutoff) * float(strictness)
    
    elements = {}

    f = open(path, 'r')
    for line in f:
        ref, mobi, dRM, dMR = line.split()
        ref = int(ref)
        mobi = int(mobi)
        dRM = float(dRM)
        dMR = float(dMR)

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
        if dRM >= cutoff and dMR >= partner_cutoff:
            r.add_neighbor(m)
        if dMR >= cutoff and dRM >= partner_cutoff:
            m.add_neighbor(r)

    f.close()

    return elements

def remove_true_singletons(element_pool):
    """ Removes from the pool elements without any neighbor """
    
    ep = element_pool

    ts = set([e for e in ep if not ep[e].neighbors])

    # Remove ts from everybody's neighbor list
    ts_e = set(ep[e] for e in ts)
    for e in element_pool:
        ep[e].neighbors = ep[e].neighbors.difference(ts_e)

    # Remove ts from pool
    for e in ts:
        del ep[e]

    return (ts, ep)

def cluster_elements(element_pool, threshold):
    """ 
    Groups Elements within a given threshold 
    together in the same cluster.
    """
    
    clusters = []
    threshold -= 1 # Account for center
    # with open('cli.dat', 'w') as fh:

    # fh.close()
    ep = element_pool
    cn = 1 # Cluster Number
    while 1:
        # Clusterable elements
        ce = [e for e in ep if not ep[e].cluster]
        if not ce: # No more elements to cluster
            break
        
        # Select Cluster Center
        # Element with largest neighbor list
        ctr_nlist, ctr = sorted([(len([se for se in ep[e].neighbors if not se.cluster]), e) for e in ce])[-1]

        # Cluster until length of remaining elements lists are above threshold
        if ctr_nlist < threshold:
            break
        
        # Create Cluster
        c = Cluster(cn, ep[ctr])
        cn += 1
        clusters.append(c)

    return (ep, clusters)

def output_clusters(handle, clusters):
    """Outputs the cluster name, center, and members."""

    write = handle.write

    for c in clusters:
        write( "Cluster %s -> %s " %(c.name, c.center.name) )
        for m in sorted(c.members):
            write( "%s " %m.name )
        write("\n")

if __name__ == "__main__":

    import optparse
    import sys
    from time import time, ctime
    import os

    USAGE="%s <matrix file> <threshold [float]> [options]" %os.path.basename(sys.argv[0])

    parser = optparse.OptionParser(usage=USAGE)
    parser.add_option('-o', '--output', dest="output_handle", action='store', type='str',
                    default=sys.stdout,
                    help='Output File [STDOUT]')
    parser.add_option('-c', '--cluster-size', dest="clus_size", action="store", type="int",
                    default=4, 
                    help="Minimum number of elements in a cluster [4]")
    parser.add_option('-s', '--strictness', dest="strictness", action="store", type='float',
                    default=0.75,
                    help="Multiplier for cutoff for M->R inclusion threshold. [0.75 or effective cutoff of 0.5625]")


    (options, args) = parser.parse_args()

    if sys.version_info[0:2] < (2,6):
        cur_version = "%s.%s" %sys.version_info[0:2]
        sys.stderr.write("- Python version not supported (%s). Please use 2.5 or newer.\n" %cur_version )
        sys.exit(1)
    if len(args) != 2:
        sys.stderr.write("- Invalid number of arguments: %i\n" %len(args))
        sys.stderr.write("USAGE: %s\n" %USAGE)
        sys.exit(1)
    
    fmatrix, cutoff = args
    cutoff = float(cutoff)
    
    # Read Matrix
    sys.stderr.write("+ BEGIN: %s\n" %ctime())
    t_init = time()

    try:
        pool = read_matrix(fmatrix, cutoff, options.strictness)
    except IOError:
        sys.stderr.write("File not found: %s\n" %fmatrix)
        sys.exit(1)

    sys.stderr.write("+ Read %ix%i distance matrix in %i seconds\n" %(len(pool), len(pool), int(time()-t_init)))
    
    # ts, pool = remove_true_singletons(pool)
    # sys.stderr.write("+ Detected %i True Singletons\n" %len(ts))

    # Cluster
    element_pool, clusters = cluster_elements(pool, options.clus_size)

    # Output Clusters
    o = options.output_handle
    if isinstance(o, str):
        o_handle = open(o, 'w')
    else:
        o_handle = o

    sys.stderr.write("+ Writing %i Clusters\n" %len(clusters))
    output_clusters(o_handle, clusters)
    if isinstance(o, str):
        o_handle.close()

    total_elements = len(element_pool)
    clustered = sum([len(c) for c in clusters])
    # Calculate coverage
    clust_coverage = clustered*100/float(total_elements)
    sys.stderr.write("+ Coverage %3.2f%% (%i/%i)\n" %(clust_coverage, clustered, total_elements))
    t_elapsed = time()-t_init
    sys.stderr.write( "+ END: %s [%3.2f seconds]\n" %(ctime(), t_elapsed))
