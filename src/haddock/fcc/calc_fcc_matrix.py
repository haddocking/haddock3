#!/usr/bin/env python
# -*- coding: UTF-8  -*-

"""
Calculates a matrix of fraction of common contacts between two or more structures.

Authors:
        RODRIGUES Joao
        TRELLET Mikael
        MELQUIOND Adrien
"""


# Contact Parsing routines
def parse_contact_file(f_list, ignore_chain):
    """Parses a list of contact files."""

    if ignore_chain:
        contacts = [
            [int(l[0:5] + l[6:-1]) for l in open(f)] for f in f_list if f.strip()
        ]
    else:
        contacts = [set([int(l) for l in open(f)]) for f in f_list if f.strip()]

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

    largest, smallest = sorted([listA, listB], key=len)
    ncommon = len([ele for ele in largest if ele in smallest])
    return (ncommon, ncommon)


# Matrix Calculation


def calculate_pairwise_matrix(contacts, ignore_chain):
    """Calculates a matrix of pairwise fraction of common contacts (FCC).
    Outputs numeric indexes.

    contacts: list_of_unique_pairs_of_residues [set/list]

    Returns pairwise matrix as an iterator, each entry in the form:
    FCC(cplx_1/cplx_2) FCC(cplx_2/cplx_1)
    """

    contact_lengths = []
    for c in contacts:
        try:
            ic = 1.0 / len(c)
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
            yield (i + 1, k + 1, fcc, fcc_v)


def _output_fcc(output, values, f_buffer):

    buf = []
    for i in values:
        buf.append(i)
        if len(buf) == f_buffer:
            output(
                "".join(["%s %s %1.3f %1.3f\n" % (i[0], i[1], i[2], i[3]) for i in buf])
            )
            buf = []
    output("".join(["%s %s %1.3f %1.3f\n" % (i[0], i[1], i[2], i[3]) for i in buf]))
