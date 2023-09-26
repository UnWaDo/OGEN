import sys
from openmm.app import ForceField, PDBReporter, Element
from openmm import LangevinMiddleIntegrator, CustomNonbondedForce
from openmm.app import PDBFile, ForceField, Modeller, Simulation, Topology
from typing import Dict, List, Tuple
from openmm.unit import kilocalorie_per_mole, picosecond, picoseconds, kelvin
import os
import argparse
from FilesProcessing import parse_pdb, Molecule


def to_openmm(mol: Molecule):
    topology = Topology()
    chain = topology.addChain()
    residue = topology.addResidue('UNK', chain)

    atoms = []
    for atom in mol.atoms:
        atoms.append(topology.addAtom(
            atom.type, Element.getBySymbol(atom.element),
            residue, atom.num))

    for bond in mol.get_bonds():
        topology.addBond(
            atoms[mol.atoms.index(bond[0])],
            atoms[mol.atoms.index(bond[1])]
        )
    positions = [a.coords / 10 for a in mol.atoms]

    return topology, positions


def OPLS_LJ(system):
    forces = {system.getForce(index).__class__.__name__: system.getForce(
        index) for index in range(system.getNumForces())}
    nonbonded_force = forces['NonbondedForce']
    lorentz = CustomNonbondedForce(
        '4*epsilon*((sigma/r)^12-(sigma/r)^6); sigma=sqrt(sigma1*sigma2); epsilon=sqrt(epsilon1*epsilon2)')
    lorentz.setNonbondedMethod(nonbonded_force.getNonbondedMethod())
    lorentz.addPerParticleParameter('sigma')
    lorentz.addPerParticleParameter('epsilon')
    lorentz.setCutoffDistance(nonbonded_force.getCutoffDistance())
    system.addForce(lorentz)
    LJset = {}
    for index in range(nonbonded_force.getNumParticles()):
        charge, sigma, epsilon = nonbonded_force.getParticleParameters(index)
        LJset[index] = (sigma, epsilon)
        lorentz.addParticle([sigma, epsilon])
        nonbonded_force.setParticleParameters(
            index, charge, sigma, epsilon * 0)
    for i in range(nonbonded_force.getNumExceptions()):
        (p1, p2, q, sig, eps) = nonbonded_force.getExceptionParameters(i)
        # ALL THE 12,13 and 14 interactions are EXCLUDED FROM CUSTOM NONBONDED
        # FORCE
        lorentz.addExclusion(p1, p2)
        if eps._value != 0.0:
            sig14 = (LJset[p1][0] * LJset[p2][0]) ** 0.5
            eps14 = (LJset[p1][1] * LJset[p2][1]) ** 0.5
            nonbonded_force.setExceptionParameters(i, p1, p2, q, sig14, eps)
    return system


def calc_energy(topology: Topology, positions: List[Tuple[float, float, float]],
                ffs: List[str], name: str, opls_lj = False) -> float:
    ff = ForceField(*ffs)

    model = Modeller(topology, positions)
    model.addExtraParticles(ff)
    sys = ff.createSystem(model.topology)
    if opls_lj:
        sys = OPLS_LJ(sys)
    integrator = LangevinMiddleIntegrator(
        300 * kelvin,
        1 / picosecond,
        0.004 * picoseconds
    )
    sim = Simulation(model.topology, sys, integrator)
    sim.context.setPositions(model.getPositions())
    if not os.path.exists('reports'):
        os.makedirs('reports')
    reporter = PDBReporter('reports/%s.pdb' % name, 1)
    reporter.report(sim, sim.context.getState(getPositions=True))
    energies = sim.context.getState(getEnergy=True)
    return (energies.getKineticEnergy() + energies.getPotentialEnergy()) / kilocalorie_per_mole

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('xml_dir', help='directory with xml files')
    parser.add_argument('out', help='name of the output file. .csv will be added by default')
    parser.add_argument('pdbs', nargs='+', help='space-separated list of .pdb files')
    args = parser.parse_args()

    pdbs = args.pdbs
    xmls = [os.path.join(args.xml_dir, f) for f in os.listdir(args.xml_dir)]
    opls_lj = not ('qube' in args.xml_dir)
    energies = {}
    for pdb in pdbs:
        try:
            name, _ = os.path.splitext(os.path.basename(pdb))

            mol = PDBFile(pdb)
            _, components, _ = parse_pdb(pdb)
            topologies = []
            positions = []
            for component in components:
                t, p = to_openmm(component)
                topologies.append(t)
                positions.append(p)
            e = calc_energy(mol.topology, mol.getPositions(),
                                         xmls, name, opls_lj)
            energies[name] = e
            diss_energy = e

            for i, (t, p) in enumerate(zip(topologies, positions)):
                component_name = '%s_mol%02d' % (name, i + 1)
                e = calc_energy(t, p, xmls,
                                component_name, opls_lj)
                energies[component_name] = e
                diss_energy -= e
            energies['%s_diss' % name] = diss_energy
        except Exception as e:
            print('Error during calculation of %s. Error is %s' % (pdb, repr(e)), file=sys.stderr)
            exit(1)
    
    out = args.out
    if not out.endswith('.csv'):
        out += '.csv'
    output = []
    for name, energy in energies.items():
        output.append('%s,%f\n' % (name, energy))
    with open(out, 'w+') as o:
        o.writelines(output)
