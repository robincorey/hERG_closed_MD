#!/usr/bin/env python3

import sys
import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array
import numpy as np

# Residues of interest
RESIDUES = [624, 652, 660]

# Chain connectivity defining the square
CHAIN_PAIRS = [
    ("A", "C"),
    ("B", "D"),
]


def min_residue_distance(universe, chain1, chain2, resid):
    """Return minimum atom-atom distance between the same residue
    in two chains."""
    
    sel1 = universe.select_atoms(f"segid {chain1} and resid {resid}")
    sel2 = universe.select_atoms(f"segid {chain2} and resid {resid}")

    if len(sel1) == 0:
        raise ValueError(f"Residue {resid} not found in chain {chain1}")

    if len(sel2) == 0:
        raise ValueError(f"Residue {resid} not found in chain {chain2}")

    dists = distance_array(sel1.positions, sel2.positions)

    return np.min(dists)


def main():

    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} structure.pdb")
        sys.exit(1)

    pdb_file = sys.argv[1]

    u = mda.Universe(pdb_file)

    for resid in RESIDUES:

        print(f"\nResidue {resid}")
#        print("-" * 40)

        for chain1, chain2 in CHAIN_PAIRS:

            dmin = min_residue_distance(
                u,
                chain1,
                chain2,
                resid,
            )

            print(
                f"{chain1}-{chain2}: {dmin:8.3f}"
            )


if __name__ == "__main__":
    main()
