#!/usr/bin/env python

"""
Outputs the names of the cluster members based on a clustering output file
and the original file listing the PDB files

Authors:
           RODRIGUES Joao
"""

import os
import sys

USAGE = "python %s <cluster_x.out> <file.nam>" %os.path.basename(sys.argv[0])

def read_clusters(path):
    """
    Reads clusters from a FCC output file.
    """

    clusters = []
    cl_file = open(path, 'r')
    for line in cl_file:
        # Cluster 8 -> 193 141 142 144 151 168 171 172 178
        models = map(int, line.split()[3:])
        clusters.append(models)

    return clusters

def read_list(path):
    """
    Reads a list containing one file per line.
    Returns an index of line number - line content
    """

    with open(path, 'r') as fhandle:
        fdata = {}
        for nline, line in enumerate(fhandle):
            if not line.strip():
                continue
            # Remove extension
            fdata[nline+1] = '.'.join(line.strip().split('.')[:-1])

    return fdata

def cross_data(clusters, flist):
    """
    Matches names in flist to the numbers in clusters.
    """

    named_clusters = []
    for cl in clusters:
        ncl = [flist[s] for s in cl]
        named_clusters.append(ncl)

    return named_clusters            

if __name__ == '__main__':

    if len(sys.argv[1:]) != 2:
        print USAGE
        sys.exit(1)


    cluster_file = os.path.abspath(sys.argv[1])
    pdblist_file = os.path.abspath(sys.argv[2])
    
    try:
        cl_list = read_clusters(cluster_file)
    except IOError:
        sys.stderr.write('Error: file not found (%s)\nAborting..\n' %cluster_file)
        sys.exit(1)

    try:
        pdb_list = read_list(pdblist_file)
    except IOError:
        sys.stderr.write('Error: file not found (%s)\nAborting..\n' %pdblist_file)
        sys.exit(1)

    named_clusters = cross_data(cl_list, pdb_list)
    
    # Output
    for i, nc in enumerate(named_clusters):
        print "Cluster %i -> %s" %(i+1, ' '.join(nc))
