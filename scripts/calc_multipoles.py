import argparse

import numpy as np
from openmm.app import ForceField

from OGEN.files.reader import read_molecule
from OGEN.openmm_api import get_charges_and_positions, mol_to_openmm
from OGEN.utils import get_coordinates, get_dipole, get_quadrupole, get_octopole

parser = argparse.ArgumentParser(
    prog='multipoles',
    description='Calculate dipole moment from structure and Openmm forcefield')
parser.add_argument('coords', help='file storing molecule coordinates')
parser.add_argument('ff', help='openmm-compatible forcefield')
args = parser.parse_args()

structure_path: str = args.coords
ff_path: str = args.ff

structure = read_molecule(structure_path)
ff = ForceField(ff_path)

topology = mol_to_openmm(structure)

atomic_coordinates = get_coordinates(structure)

charges, positions = get_charges_and_positions(topology, atomic_coordinates,
                                               ff)
dipole = get_dipole(positions, charges)
quadrupole = get_quadrupole(positions, charges)
# octapole = get_octopole(positions, charges)

print(np.linalg.norm(dipole))
print(dipole)

print(np.linalg.norm(quadrupole))
print(quadrupole)

# print(np.linalg.norm(octapole))
# print(octapole)
