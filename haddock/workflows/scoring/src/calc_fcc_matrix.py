#!/usr/bin/env python


#   Authors:TRELLET Mikael                  
#           RODRIGUES Joao                  
#           MELQUIOND Adrien                
#                                           
# + Undocumented feature: empty contact lists are treated as 0 contacts. (230412 JR)

"""
Calculates a matrix of fraction of common contacts between two or more structures.

WARNING!!
Structures are named numerically in the matrix by the order in which they are read.
By default we sort them by the last number (see _fit*pdb files in it1/analysis/ for example)
but a specific order (say, file.nam) is also possible using the -f option and a custom file.
"""

# Contact Parsing routines
def parse_contact_file(f_list, ignore_chain):
    """Parses a list of contact files."""
    
    if ignore_chain:
        contacts = [ [ int(l[0:5]+l[6:-1]) for l in open(f)] for f in f_list if f.strip()]
    else:
        contacts = [ set([ int(l) for l in open(f)]) for f in f_list if f.strip()]

    return contacts

# FCC Calculation Routine
def calculate_fcc(listA, listB):
    """
    Calculates the fraction of common elements between two lists
    taking into account chain IDs
    """
    
    cc = len(listA.intersection(listB))
    cc_v = len(listB.intersection(listA))
    
    return (cc, cc_v)
    
def calculate_fcc_nc(listA, listB):
    """
    Calculates the fraction of common elements between two lists
    not taking into account chain IDs. Much Slower.
    """
    
    largest,smallest = sorted([listA, listB], key=len)
    ncommon = len([ele for ele in largest if ele in smallest])
    return (ncommon, ncommon)

# Matrix Calculation

def calculate_pairwise_matrix(contacts, ignore_chain):
    """ Calculates a matrix of pairwise fraction of common contacts (FCC).
        Outputs numeric indexes.

        contacts: list_of_unique_pairs_of_residues [set/list]
        
        Returns pairwise matrix as an iterator, each entry in the form:
        FCC(cplx_1/cplx_2) FCC(cplx_2/cplx_1)
    """

    #contact_lengths = [1/float(len(i)) for i in contacts]
    contact_lengths = []
    for c in contacts:
      try:
        ic = 1.0/len(c)
      except ZeroDivisionError:
        ic = 0
      contact_lengths.append(ic)
    
    if ignore_chain:
        calc_fcc = calculate_fcc_nc
    else:
        calc_fcc = calculate_fcc
    
    for i in range(len(contacts)):

        for k in range(i+1, len(contacts)):
            cc, cc_v = calc_fcc(contacts[i], contacts[k])
            fcc, fcc_v = cc*contact_lengths[i], cc*contact_lengths[k]
            yield (i+1, k+1, fcc, fcc_v)
            # yield (i+1, k+1, len(contacts[i].intersection(contacts[k]))*contact_lengths[i], len(contacts[k].intersection(contacts[i]))*contact_lengths[k])

def _output_fcc(output, values, f_buffer):

    c = 0
    buf = []
    for i in values:
        buf.append(i)
        if len(buf) == f_buffer:
            output( ''.join(["%s %s %1.3f %1.3f\n" %(i[0],i[1],i[2],i[3]) for i in buf]) )
            buf = []
    output( ''.join(["%s %s %1.3f %1.3f\n" %(i[0],i[1],i[2],i[3]) for i in buf]) )
    
if __name__ == '__main__':
    
    import optparse
    import sys
    from time import time, ctime
    import os
    
    USAGE = "%s <contacts file 1> <contacts file 2> ... [options]\n" %os.path.basename(sys.argv[0])
    
    parser = optparse.OptionParser(usage=USAGE)
    parser.add_option('-o', '--output', dest="output_file", action='store', type='string',
                        default=sys.stdout,
                        help='Output File [default: STDOUT]')
    parser.add_option('-f', '--file', dest="input_file", action='store', type='string',
                        help='Input file (one contact file name per line)')
    parser.add_option('-b', '--buffer_size', dest="buffer_size", action='store', type='string',
                        default=50000, 
                        help='Buffer size for writing output. Number of lines to cache before writing to file [default: 50000]')
    parser.add_option('-i', '--ignore_chain', dest="ignore_chain_char", action='store_true',
                        help='Ignore chain character in residue code. Use for homomeric complexes.')
    parser.add_option('-H', '--haddock-run', dest="HADDOCK", action='store_true',
                        help='Touches a final MTX_DONE file. Only relevant for internals of the HADDOCK software.')
 
    (options, args) = parser.parse_args()
    
    if options.input_file: # take order as is
        args = [name.strip() for name in open(options.input_file)]
    else: # Assume _fit files and sort numerically
        args = sorted(args, key=lambda x: int(x.split('_')[-1].split('.')[0]))

    if len(args) < 2:
        sys.stderr.write("- Provide (at least) two structures to calculate a matrix. You provided %s.\n" %len(args))
        sys.stderr.write(USAGE)
        sys.exit(1)

    sys.stderr.write("+ BEGIN: %s\n" %ctime())
    if options.ignore_chain_char:
        sys.stderr.write("+ Ignoring chains. Expect a considerable slowdown!!\n")
        exclude_chains = True
    else:
        exclude_chains = False
        
    t_init = time()
    sys.stderr.write("+ Parsing %i contact files\n" %len(args))

    c = parse_contact_file(args, exclude_chains)
    
    m = calculate_pairwise_matrix(c, exclude_chains)
    
    if isinstance(options.output_file, str):
        f = open(options.output_file, 'w')
    else:
        f = options.output_file

    sys.stderr.write("+ Calculating Matrix\n") # Matrix is calculated when writing. Generator property.
    sys.stderr.write("+ Writing matrix to %s\n" %f.name)
    _output_fcc(f.write, m, options.buffer_size)
        
    if isinstance(options.output_file, str):
        f.close()
    t_elapsed = time()-t_init
    sys.stderr.write("+ END: %s [%6.2f seconds elapsed]\n" %(ctime(), t_elapsed))
    if options.HADDOCK:
        rundir = os.path.dirname(options.output_file)
        touchFile = open(os.path.join(rundir, 'MTX_DONE'), 'w')
        touchFile.close()
