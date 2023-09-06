import sys, os
import argparse
import numpy as np
import itertools
import pandas as pd
import matplotlib.pyplot as plt
from openmm import LangevinMiddleIntegrator
from openmm.app import PDBFile, ForceField, Modeller, Simulation
from openmm.unit import kilocalorie_per_mole, picosecond, picoseconds, kelvin, angstrom, elementary_charge


from OGEN.utils import get_dipole


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
    description = 'Calculate dipole moment from structure and Openmm forcefield'
)
parser.add_argument('pdb', help='pdb file storing molecule coordinates')
parser.add_argument('ff', help='openmm-compatible forcefield')
args = parser.parse_args()

structure_path = str(args.pdb)
ff_path = str(args.ff)

ff = ForceField(ff_path)
charges, positions = get_charges_and_positions(structure_path, ff)
dipole = get_dipole(positions, charges)

print(np.linalg.norm(dipole))
print(dipole)
