import argparse
import numpy as np

from openmm import LangevinMiddleIntegrator
from openmm.app import PDBFile, ForceField, Modeller, Simulation
from openmm.unit import kilocalorie_per_mole, picosecond, picoseconds, kelvin, angstrom, elementary_charge

from OGEN.utils.Constants import BOHR_TO_ANGSTROM


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


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog = 'calc_quadrupole',
        description = 'Calculate quadrupole moment from structure and Openmm forcefield'
    )
    parser.add_argument('pdb', help='pdb file storing molecule coordinates')
    parser.add_argument('ff', help='openmm-compatible forcefield')
    parser.add_argument('esps', help='file with ESP points in numpy format')
    args = parser.parse_args()

    sample_points = np.loadtxt(args.esps)
    esps = sample_points[:,3]
    points = sample_points[:,:-1] * BOHR_TO_ANGSTROM

    ff = ForceField(args.ff)
    charges, positions = get_charges_and_positions(args.pdb, ff)

    invr = np.zeros((len(points), len(positions)))
    for i in range(invr.shape[0]):
        for j in range(invr.shape[1]):
            invr[i, j] = 1/np.linalg.norm(points[i]-positions[j])
    invr *= BOHR_TO_ANGSTROM # convert to atomic units

    x_esp = np.einsum('ij, j -> i', invr, charges)
    print(np.square(x_esp - esps).sum())
