import io
import os
import re
from datetime import datetime
from string import Template
from typing import Dict, List, Tuple

import numpy as np
import numpy.typing as npt
from openbabel import openbabel as ob
from openbabel import pybel
from openmm import *
from openmm.app import ForceField, Modeller, PDBFile, Simulation, Topology, Element
from openmm.unit import kelvin, kilocalorie_per_mole, picosecond, picoseconds, elementary_charge, angstrom

conf_dt = np.dtype([('A', np.int16), ('B', np.int16), ('E', np.float64)])
ERR_DIFF_SIZE = 'Energies and references has different sizes: %d and %d'
ERR_UNMATCHING = 'Energies and reference has unmatching points'


def create_ff(
    template: Template,
    coords: List[float],
    charges: List[float]
) -> ForceField:
    '''
        coords: list of coords used in template [must be in template units]
        charges: list of charges
    '''
    data = {'charge_%d' % i: c for i, c in enumerate(charges)}
    data['date'] = datetime.today().strftime('%Y_%m_%d')
    data.update({'coord_%d' % i: c for i, c in enumerate(coords)})
    ff = template.substitute(data)
    return ForceField(io.StringIO(ff))


def write_pdb(ff: ForceField, structure: str, dest: str) -> None:
    mol = PDBFile(structure)
    model = Modeller(mol.topology, mol.getPositions())
    model.addExtraParticles(ff)
    with open(dest, 'w') as d:
        PDBFile.writeFile(model.getTopology(), model.getPositions(), file=d)


def write_pdbs(
    ff: ForceField,
    structures: Dict[Tuple[int, int], str]
) -> None:
    if not os.path.isdir('opt_pdb'):
        os.mkdir('opt_pdb')
    for name in structures:
        write_pdb(
            ff,
            structures[name],
            'opt_pdb/%d_%d.pdb' % (name[0], name[1])
        )


def calc_energy(ff: ForceField, structure: str) -> float:
    mol = PDBFile(structure)
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
    energies = sim.context.getState(getEnergy=True)
    return (
        energies.getKineticEnergy() + energies.getPotentialEnergy()
    ) / kilocalorie_per_mole


def calc_energies(
    ff: ForceField,
    structures: Dict[Tuple[int, int], str]
) -> npt.ArrayLike:
    energies = [(
        name[0],
        name[1],
        calc_energy(ff, structures[name])
    ) for name in structures]
    energies = np.array(energies, dtype=conf_dt)
    for i, e in enumerate(energies):
        if np.isnan(e['E']):
            energies[i]['E'] = np.nanmax(energies['E']) * 1000
    return energies


def calc_error(
    energies: npt.ArrayLike,
    reference: npt.ArrayLike
) -> Tuple[float, float]:
    if len(energies) != len(reference):
        raise Exception(ERR_DIFF_SIZE % (len(energies), len(reference)))
    reference.sort(order=['A', 'B'])
    energies.sort(order=['A', 'B'])
    if (not np.array_equal(reference['A'], energies['A']) or
            not np.array_equal(reference['B'], energies['B'])):
        raise Exception(ERR_UNMATCHING)
    min_idx = np.argmin(reference['E'])
    calc = energies['E'] - energies['E'][min_idx]
    ref = reference['E'] - reference['E'][min_idx]
    errors = np.abs(calc - ref)
    return errors.mean(), errors.max()


def load_structures(folder: str) -> Dict[Tuple[int, int], str]:
    angles = re.compile(r'(\d+)_(\d+)')
    files = os.listdir(folder)
    structures = {}
    for f in files:
        m = angles.search(f)
        if m is None:
            raise Exception('Invalid file %s' % f)
        structures[(
            int(m.group(1)),
            int(m.group(2))
        )] = os.path.join(folder, f)
    return structures


def load_reference(file: str) -> npt.ArrayLike:
    return np.loadtxt(file, dtype=conf_dt, delimiter=',')


def prepare_parameters(
    template_path: str,
    structures_path: str,
    reference_path: str
) -> Tuple[
        Template,
        Dict[Tuple[int, int], str],
        npt.ArrayLike
]:
    with open(template_path, 'r') as t:
        template = ''.join(t.readlines())
    template = Template(template)
    structures = load_structures(structures_path)
    reference = load_reference(reference_path)
    return template, structures, reference


def parse_structure(path: str) -> Tuple[List[str], npt.ArrayLike]:
    symbols = []
    coords = []
    with open(path, 'r') as sf:
        for line in sf:
            p = line.split()
            symbols.append(p[0])
            coords.append([float(a) for a in p if not a.isalpha()])
    return symbols, np.array(coords)

def get_charges_and_positions(topology: Topology, positions: np.ndarray,
                              forcefield: ForceField):
    """
        positions: np.ndarray in angstroms
    """
    model = Modeller(topology, positions / 10)
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


def mol_to_openmm(mol: pybel.Molecule) -> Topology:
    topology = Topology()
    chain = topology.addChain()
    residue = topology.addResidue('UNK', chain)

    atoms = []
    for atom in mol:
        atoms.append(
            topology.addAtom(atom.type,
                Element.getByAtomicNumber(atom.atomicnum), residue, atom.idx))

    for bond in ob.OBMolBondIter(mol.OBMol):
        begin = bond.GetBeginAtomIdx()
        end = bond.GetEndAtomIdx()
        topology.addBond(atoms[begin - 1], atoms[end - 1])

    return topology
