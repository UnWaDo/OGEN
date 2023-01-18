import sys
from openmm.app import ForceField, PDBReporter
from openmm import LangevinMiddleIntegrator
from openmm.app import PDBFile, ForceField, Modeller, Simulation
from typing import Dict, List, Tuple
from openmm.unit import kilocalorie_per_mole, picosecond, picoseconds, kelvin
import os


def calc_energy(structure_path: str, ffs: List[str]) -> float:
    ff = ForceField(*ffs)

    mol = PDBFile(structure_path)
    model = Modeller(mol.topology, mol.getPositions())
    model.addExtraParticles(ff)
    sys = ff.createSystem(model.topology)
    integrator = LangevinMiddleIntegrator(
        300 * kelvin,
        1 / picosecond,
        0.004 * picoseconds
    )
    sim = Simulation(model.topology, sys, integrator)
    sim.context.setPositions(model.getPositions())
    name, _ = os.path.splitext(os.path.basename(structure_path))
    if not os.path.exists('reports'):
        os.makedirs('reports')
    reporter = PDBReporter('reports/%s.pdb' % name, 1)
    reporter.report(sim, sim.context.getState(getPositions=True))
    energies = sim.context.getState(getEnergy=True)
    return (energies.getKineticEnergy() + energies.getPotentialEnergy()) / kilocalorie_per_mole

if __name__ == '__main__':
    if len(sys.argv) < 3:
        print('Not enough arguments, pass .pdb file and .xml files (e.g. file.pdb xmls/*.xml)', file=sys.stderr)
        exit(1)
    pdb = sys.argv[1]
    xmls = sys.argv[2:]
    print(calc_energy(pdb, xmls))
