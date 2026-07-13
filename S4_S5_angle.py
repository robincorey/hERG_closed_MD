#!/usr/bin/env python3

import sys
import numpy as np
import MDAnalysis as mda

CHAINS = ["A", "B", "C", "D"]

# Vector definitions
V1_START = 118
V1_END   = 140

V2_START = 175
V2_END   = 150


def get_ca_position(u, chain, resid):

    sel = u.select_atoms(
        f"chainID {chain} and resid {resid} and name CA"
    )

    if len(sel) != 1:
        raise ValueError(
            f"Could not find unique CA atom for residue "
            f"{resid} in chain {chain}"
        )

    return sel.positions[0]


def angle_between(v1, v2):

    v1 = v1 / np.linalg.norm(v1)
    v2 = v2 / np.linalg.norm(v2)

    cos_theta = np.dot(v1, v2)

    # numerical safety
    cos_theta = np.clip(cos_theta, -1.0, 1.0)

    theta = np.degrees(np.arccos(cos_theta))

    return theta


def chain_angle(u, chain):

    p518 = get_ca_position(u, chain, V1_START)
    p540 = get_ca_position(u, chain, V1_END)

    p550 = get_ca_position(u, chain, V2_START)
    p575 = get_ca_position(u, chain, V2_END)

    v1 = p540 - p518
    v2 = p575 - p550

    return angle_between(v1, v2)


def main():

    if len(sys.argv) != 2:
        print(f"Usage: {sys.argv[0]} structure.pdb")
        sys.exit(1)

    u = mda.Universe(sys.argv[1])

    for chain in CHAINS:

        angle = chain_angle(u, chain)

        print(
            f"{angle:7.2f}"
        )


if __name__ == "__main__":
    main()
