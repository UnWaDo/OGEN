import sys, os
import argparse
import numpy as np
import itertools
import pandas as pd
import matplotlib.pyplot as plt
from openmm import LangevinMiddleIntegrator
from openmm.app import PDBFile, ForceField, Modeller, Simulation
from openmm.unit import kilocalorie_per_mole, picosecond, picoseconds, kelvin, angstrom, elementary_charge


D_TO_CM = 1e-21 / 299792458
AU_ANSGSTROM_TO_CM = 1.6e-19*1e-10
BOHR_TO_ANGSTROM = 0.529177249


def get_charges_and_positions(structure_path, forcefield):
    mol = PDBFile(structure_path)
    model = Modeller(mol.topology, mol.getPositions())
    model.addExtraParticles(forcefield)
    system = forcefield.createSystem(model.topology)
    integrator = LangevinMiddleIntegrator(
        300 * kelvin,
        1 / picosecond,
        0.004 * picoseconds
    )
    sim = Simulation(model.topology, system, integrator)
    sim.context.setPositions(model.getPositions())
    forces = {system.getForce(index).__class__.__name__: system.getForce(
        index) for index in range(system.getNumForces())}
    nonbonded_force = forces['NonbondedForce']
    charges = np.array([nonbonded_force.getParticleParameters(i)[0] / elementary_charge for i in range(nonbonded_force.getNumParticles())])
    positions = np.array(sim.context.getState(getPositions=True).getPositions() / angstrom)
    return charges, positions


parser = argparse.ArgumentParser(
    prog = 'gen_cube',
    description = 'Generate Gaussian cube file based on existing one recalculating ESP from MM'
)
parser.add_argument('existing_cube', help='cube file with points')
parser.add_argument('pdb', help='pdb file storing molecule topology')
parser.add_argument('ff', help='openmm-compatible forcefield')
parser.add_argument('--out', help='Name of out cube file', default='mm.cube')
args = parser.parse_args()

cube_path = str(args.existing_cube)
cube_out_path = str(args.out)
structure_path = str(args.pdb)
ff_path = str(args.ff)

ff = ForceField(ff_path)
charges, positions = get_charges_and_positions(structure_path, ff)

with open(cube_path, 'r') as cube_in:
    cube_in.readline()
    cube_in.readline()
    data = cube_in.readline().split()
    n_atoms = int(data[0])
    origin = np.array([float(d) for d in data[1:]]) * BOHR_TO_ANGSTROM
    amounts = []
    axis = []
    for i in range(3):
        data = cube_in.readline().split()
        amounts.append(int(data[0]))
        axis.append([float(d) for d in data[1:]])
    axis = np.array(axis) * BOHR_TO_ANGSTROM
    atoms = []
    for i in range(n_atoms):
        data = cube_in.readline().split()
        atoms.append(np.array([int(data[0])] + [float(d) * BOHR_TO_ANGSTROM for d in data[2:]]))
atoms = np.array(atoms)

with open(cube_out_path, 'w+') as cube_out:
    cube_out.write('First comment line\n Electrostatic potential from Total SCF Density\n')
    cube_out.write('{:6d} {:11.6f} {:11.6f} {:11.6f}\n'.format(n_atoms, *(origin / BOHR_TO_ANGSTROM)))
    for i in range(3):
        cube_out.write('{:6d} {:11.6f} {:11.6f} {:11.6f}\n'.format(amounts[i], *(axis[i] / BOHR_TO_ANGSTROM)))
    for i in range(n_atoms):
        cube_out.write('{:6d} {:11.6f} {:11.6f} {:11.6f} {:11.6f}\n'.format(
            int(atoms[i][0]),
            atoms[i][0],
            *(atoms[i][1:] / BOHR_TO_ANGSTROM))
        )
    nums_in_line = []
    for idx in itertools.product(*[range(a) for a in amounts]):
        coords = origin + axis @ np.array(idx)
        distances = np.linalg.norm(positions - coords, axis=1) / BOHR_TO_ANGSTROM
        nums_in_line.append((1 / distances) @ charges)
        if len(nums_in_line) == 6 or idx[2] == amounts[2] - 1:
            out = ''.join(['%13.5E' % n for n in nums_in_line])
            cube_out.write(out + '\n')
            nums_in_line = []
