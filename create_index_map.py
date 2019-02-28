#!/usr/bin/env python

import MDAnalysis as mda
from MDAnalysis.analysis import align
from MDAnalysis.analysis.rms import rmsd
import numpy as np


def rotate(Umobile, Ureference, selection):
    """Rotate the a system onto a reference.

    Args:
        Umobile (Universe): Universe that should be rotated to fit a reference
        Ureference (Universe): Universe that is used as a reference for RMSD alignment
        selection (str): The selection that will be used for the alignment.

    Returns:
        float: RMSD of the aligned structures across the selection

    """
    mobile0 = Umobile.select_atoms(selection).positions - Umobile.atoms.center_of_mass()
    ref0 = Ureference.select_atoms(selection).positions - Ureference.atoms.center_of_mass()
    R, rmsd = align.rotation_matrix(mobile0, ref0)
    Umobile.atoms.translate(-Umobile.select_atoms(selection).center_of_mass())
    Umobile.atoms.rotate(R)
    Umobile.atoms.translate(Ureference.select_atoms(selection).center_of_mass())
    return rmsd

def create_index_map(U1, U2):
    """Output a dictionary of map between U1 and U2. The keys will be indices from U1
    and the values will the corresponding indices from U2. The two structures should be
    aligned and almost exactly the same.

    Args:
        U1 (Universe): First structure (should be the gro file)
        U2 (Universe): Second structure (should be the dms file)

    Returns:
        dict: a dictionary of map between U1 and U2. The keys will be indices from U1
              and the values will the corresponding indices from U2
    """

    def find_index(atom, candidates):
        """Given an atom in U1, find the index in U2
        
        Args:
            atom (Atom): Atom from U1
            candidates (list): List of atoms from U2

        Returns:
            int: Index from U2 that matches atom
        """
        N = 0
        for index, a in enumerate(candidates):
            sep = np.sqrt(sum([i**2 for i in (a.position - atom.position)]))
            if sep < 0.01:
                if N > 2:
                    first_half = list(candidates[0:index])
                    second_half = list(candidates[index:])
                    candidates[0:len(second_half)] = second_half[:]
                    candidates[len(second_half):] = first_half[:]
                    return candidates.pop(0).index
                return candidates.pop(index).index
            N += 1
        raise ValueError

    mapping = {}
    
    candidates = [a for a in U2.atoms]

    for a1 in U1.atoms:
        index1 = a1.index
        index2 = find_index(a1, candidates)
        mapping[index1] = index2

    return mapping

if __name__ == "__main__":
    import sys
    import json

    if "-h" in sys.argv:
        help_msg = """example usage: create_index_map.py reference_structure.gro reference_structure.dms "protein" """
        print(help_msg)
        exit(1)
    
    Umobile = mda.Universe(sys.argv[1])
    Ureference = mda.Universe(sys.argv[2])

    # it really doesn't matter which structure is rotated
    rmsd_alignment = rotate(Umobile, Ureference, sys.argv[3])
    if rmsd_alignment > 0.1:
        print(rmsd_alignment)
        raise ValueError

    mapping = create_index_map(Umobile, Ureference)

    if not len(set(mapping.keys())) == len(mapping.keys()):
        raise ValueError("Repeated indices in the keys")

    if not len(set([mapping[i] for i in mapping.keys()])) == len(mapping.keys()):
        raise ValueError("Repeated indices in the values")

    with open("map.json", 'w') as f:
        f.write(json.dumps(mapping))

    print("Keys are from {0} and values are from {1}".format(sys.argv[1], sys.argv[2]))
