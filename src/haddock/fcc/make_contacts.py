#!/usr/bin/env python
# -*- coding: UTF-8  -*-

"""
Script to calculate contact lists on PDB files.
Requires external executable to calculate the contacts!

Authors:
           RODRIGUES Joao
           TRELLET Mikael
           MELQUIOND Adrien
"""

from multiprocessing import Process
from subprocess import Popen, PIPE

def _calculate_contacts(executable, pdbfile, d_cutoff, filter_selection=None, extension='.contacts'):
    """
    Outputs a list of contacts based on vector analysis
    of the PDB file.
    
    Arguments:
    executable  - path to contact calculation program
    pdbfile     - path to PDB-formatted file (.pdb extension)
    d_cutoff    - minimal distance in A to consider a contact (float)
    filter_selection - list of identifiers to filter contacts (list of strings)
    """
    
    pdbname = os.path.basename(pdbfile)[:-4]

    p = Popen([executable, pdbfile, d_cutoff], stdout=PIPE)
    p_output = p.communicate()[0]
    contacts = sorted(list(set([l for l in p_output.split('\n')][:-1])))
    
    # Filter contacts
    if filter_selection:
        contacts = filter(lambda x: x[5] in filter_selection and x[-1] in filter_selection, contacts)
    # 

    outfile = os.path.join(os.path.dirname(pdbfile), "%s%s" %(pdbname, extension))
    with open(outfile, 'w') as o:
        o.write('\n'.join(contacts))

    return 0
    
if __name__ == '__main__':
    
    import optparse
    import os, sys

    USAGE = "%s [-f structures.txt] [-n 4] [-c 5.0] file1.pdb file2.pdb"

    parser = optparse.OptionParser(usage=USAGE)
    parser.add_option('-c', '--cutoff', dest="d_cutoff", action='store', type='string',
                        default="5.0", 
                        help='Distance cutoff to evaluate contacts. [default: 5.0A]')
    parser.add_option('-f', '--file', dest="input_file", action='store', type='string',
                        help='Input file (one file path per line)')
    parser.add_option('-n', '--nproc', dest="nproc", action='store', type='string',
                        default=1, 
                        help='Number of simultaneous processes to launch in each round. [default: 1]')
    parser.add_option('-e', '--exec', dest="executable", action='store', type='string',
                        default='%s/../src/contact_fcc' %os.path.dirname(sys.argv[0]), 
                        help='Path to the executable C++ program to calculate contacts [default: ../fcc/src/contact_fcc]')
    parser.add_option('-s', '--selection', dest="selection", action='store', type='string',
                        default=None, 
                        help='Filter contacts based on their segids. [Default: No filtering. All chains] [Example: A,C]')

    (options, args) = parser.parse_args()

    if options.input_file:
        args = [name.strip() for name in open(options.input_file) if name.strip()]
    
    if not args:
        print "No files provided. Exiting"
        print USAGE
        sys.exit(1)
    
    # Convert to full paths
    args = map(os.path.abspath, args)
    
    nproc = int(options.nproc)
    cutoff = options.d_cutoff

    executable = options.executable
    if not os.path.exists(executable):
        print "Path not found: %s" %os.path.abspath(executable)
        sys.exit(1)
    executable = os.path.abspath(executable)

    if options.selection:
        filter_selection = set(options.selection.split(','))
        representative = open(args[0])
        repr_chains = dict([(j,str(i)) for i,j in enumerate(sorted(set([l[72] for l in representative if l.startswith('ATOM')])), start=1)])
        filter_selection = map(repr_chains.get, filter_selection)
        representative.close()
        oextension = '.contacts-'+''.join(options.selection.split(','))
    else:
        filter_selection = None
        oextension = '.contacts'

    queue = []
    
    while 1:
        
        arg = args.pop()        
        # Create Process for arg
        p = Process(target=_calculate_contacts, args=(executable, arg, cutoff, filter_selection, oextension))
        queue.append(p)
        
        if (len(queue) == nproc) or (not args and len(queue)):

            for job in queue:
                job.start()
            for job in queue: # Waiting for job to finish
                job.join()
            queue = []
        
        if not args and not queue:
            break
    
    print "Finished"
