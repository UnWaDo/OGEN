import argparse
import os

import numpy as np
from OGEN.ForceField import (AtomType, FFAtom, FFBond, ForceField, Torsion,
                             PeriodicTorsion, HarmonicAngle, HarmonicBond, NonbondedAtom,
                             Residue)
from OGEN.ForceField.VSites import AverageTwo, LocalCoords
from OGEN.MultiWFN import SpaceFunctions, calculate_rsf_at_points, get_multipoles
from OGEN.OPLS import calc_opls_parameters
from OGEN.Points.la_selectors import are_colinear
from OGEN.Points.selectors import gen_oscs
from OGEN.RESP import fit_charges, gen_points
from OGEN.files.reader import read_molecule
from OGEN.files.writer import to_xyz
from openbabel import pybel
from rdkit import Chem


HELP_FCHK = 'path to .fchk (or other wavefunction) file'
HELP_RSF = 'real-space function for off-site charges generation. Default is ELF'
HELP_OSCS_COUNT = 'amount of OSCs to put on oxygen and sulfur atoms. Default is 3, which puts 2 OSCs on carbonyl oxygen and 1 OSCs elsewhere'
HELP_FLUORINE = 'whether to generate OSCs along the C−F bond. Default is False'
HELP_NO_OSCS = 'toggle this parameter to disable OSCs generation'
HELP_REUSE = 'whether to store files, generated during the procedure. Can speed up recalculations but makes folder a bit messy'
HELP_NO_RESP = 'disables OSCs calculations and charges evaluation using RESP'
HELP_USE_BFGS = 'whether to use BFGS to solve for charges instead of standard RESP'
HELP_NO_RING = 'whether to disable OSC in the middle of the aromatic rings'
HELP_USE_MULTIPOLES = 'whether to use multipoles to fit charges'
HELP_WEIGHT_DIPOLE = 'which weight to use for dipole'
HELP_WEIGHT_QUADRUPOLE = 'which weight to use for quadrupole'
HELP_ONLY_CHARGES = 'disables calculation of parameters other than charges'

parser = argparse.ArgumentParser()
parser.add_argument('fchk', help=HELP_FCHK)
parser.add_argument('--rsf', default='ELF', help=HELP_RSF, choices=['ELF', 'Lap', 'LOL'])
parser.add_argument('--oscs-count', type=int, default=3, help=HELP_OSCS_COUNT, choices=[1, 2, 3])
parser.add_argument('--fluorine', action='store_true', help=HELP_FLUORINE)
parser.add_argument('--no-oscs', action='store_true', help=HELP_NO_OSCS)
parser.add_argument('--no-reuse', action='store_true', help=HELP_REUSE)
parser.add_argument('--no-resp', action='store_true', help=HELP_NO_RESP)
parser.add_argument('--use-bfgs', action='store_true', help=HELP_USE_BFGS)
parser.add_argument('--no-ring', action='store_true', help=HELP_NO_RING)
parser.add_argument('--use-multipoles', action='store_true', help=HELP_USE_MULTIPOLES)
parser.add_argument('--dipole-weight', type=float, default=1, help=HELP_WEIGHT_DIPOLE)
parser.add_argument('--quadrupole-weight', type=float, default=1e-2, help=HELP_WEIGHT_QUADRUPOLE)
parser.add_argument('--only-charges', action='store_true', help=HELP_ONLY_CHARGES)
args = parser.parse_args()

wavefunction_file = os.path.abspath(args.fchk)
name, _ = os.path.splitext(os.path.basename(wavefunction_file))

pybel_molecule = read_molecule(wavefunction_file)

working_dir = f'OGEN_{name}'
if not os.path.exists(working_dir):
    os.mkdir(working_dir)
os.chdir(working_dir)

mol_file = f'{name}.mol'
pybel_molecule.write('mol', mol_file, overwrite=True)
mol = Chem.MolFromMolFile(mol_file, removeHs=False, sanitize=False)
conformer = mol.GetConformer()

if args.no_oscs or args.no_resp:
    atom_oscs = []
    ring_oscs = []
else:
    atom_oscs, ring_oscs = gen_oscs(
        fchk_file = wavefunction_file,
        mol = mol,
        mode_rsf = args.rsf,
        mode_number = args.oscs_count,
        mode_fluorine = args.fluorine,
        reuse=True
    )
    if args.no_ring:
        ring_oscs = []

if args.no_oscs or args.no_resp:
    xyz_name = f'{name}_no.xyz'
else:
    xyz_name = (f'{name}_{args.rsf}{args.oscs_count}'
        f'{"F" if args.fluorine else ""}'
        f'{"NR" if args.no_ring else ""}.xyz')

to_xyz(xyz_name, mol, atom_oscs, ring_oscs, name)

points_filename = f'{name}.esp'

if os.path.exists(points_filename):
    points = np.loadtxt(points_filename)

elif not args.no_resp:
    points = gen_points(
        [conformer.GetAtomPosition(i) for i in range(len(mol.GetAtoms()))],
        [a.GetSymbol().upper() for a in mol.GetAtoms()]
    )
    points = calculate_rsf_at_points(
        wavefunction_file,
        points,
        SpaceFunctions.ESP,
    )
    np.savetxt(points_filename, points)

if args.only_charges:
    disp = [
        {'sigma': 1, 'epsilon': 0, 'charge': 0}
        for _ in range(
            len(mol.GetAtoms()) + len(atom_oscs) + len(ring_oscs)
        )
    ]
    bonds_energy = []
    angles_energy = []
    tors_energy = []
else:
    disp, bonds_energy, angles_energy, tors_energy = calc_opls_parameters(mol)


if not args.no_resp:

    multipoles = {}
    if args.use_multipoles:
        multipoles = get_multipoles(wavefunction_file, not args.no_reuse)

    qf, errors = fit_charges(
        symbols = [a.GetSymbol().upper() for a in mol.GetAtoms()],
        coords = np.array([conformer.GetAtomPosition(i) for i in range(len(mol.GetAtoms()))]),
        sample_points = points,
        extra = [o for _, o in atom_oscs + ring_oscs],
        use_bfgs=args.use_bfgs,
        multipoles=multipoles,
        multipoles_weights={
            'dipoles': args.dipole_weight,
            'quadrupoles': args.quadrupole_weight,
        },
    )
    atom_charges = qf[1][:len(mol.GetAtoms())]
    oscs_charges = qf[1][len(mol.GetAtoms()):]
else:
    atom_charges = [d['charge'] for d in disp]
    oscs_charges = []

ff = ForceField()
xml_atoms = []
for a in mol.GetAtoms():
    i = a.GetIdx()
    c = conformer.GetAtomPosition(i)
    ff.atom_types.append(AtomType(
        name = 'at%s_%02d' % (a.GetSymbol(), i),
        class_name = 'cl%s_%02d' % (a.GetSymbol(), i),
        element = a.GetSymbol(),
        mass = a.GetMass()
    ))
    ff.nonbonded_forces.append(NonbondedAtom(
        charge = atom_charges[i],
        sigma = disp[i]['sigma'],
        epsilon = disp[i]['epsilon'],
        atom_type = ff.atom_types[-1]
    ))
    xml_atoms.append(FFAtom(
        name = '%s%d' % (a.GetSymbol(), i),
        atom_type = ff.atom_types[-1]
    ))
xml_bonds = []
for b in mol.GetBonds():
    a1, a2 = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
    xml_bonds.append(FFBond(indices=(a1, a2)))
for be in bonds_energy:
    ff.bond_forces.append(HarmonicBond(
        length = be['r_eq'],
        k = be['k'],
        classes = tuple(ff.atom_types[i].class_name for i in be['atoms'])
    ))
for ae in angles_energy:
    ff.angle_forces.append(HarmonicAngle(
        angle = ae['a_eq'],
        k = ae['k'],
        classes = tuple(ff.atom_types[i].class_name for i in ae['atoms'])
    ))
for te in tors_energy:
    if not sum(te['vs']):
        continue
    ff.torsion_forces.append(PeriodicTorsion(
        params = [Torsion(
            phase = np.pi if i % 2 else 0,
            k = v,
            periodicity = i + 1
        ) for i, v in enumerate(te['vs'])],
        classes = tuple(ff.atom_types[i].class_name for i in te['atoms']),
        is_improper = te['improper']
    ))

ff.residues.append(Residue('UNL', xml_atoms, xml_bonds))

for i, (ai, osc) in enumerate(atom_oscs):
    neighb = [a.GetIdx() for a in mol.GetAtomWithIdx(ai).GetNeighbors()]
    neighb = [ai] + neighb
    atom_coords = [conformer.GetAtomPosition(n) for n in neighb]
    if len(neighb) == 2 and not are_colinear([osc] + atom_coords):
        neighb.extend([
            a.GetIdx()
                for a in mol.GetAtomWithIdx(neighb[-1]).GetNeighbors()
                    if a.GetIdx() not in neighb
        ])
        atom_coords = [conformer.GetAtomPosition(n) for n in neighb]
    if len(neighb) == 2:
        vs_coords = AverageTwo.from_coordinates(
            osc,
            atoms = neighb,
            atom_coords = atom_coords
        )
    else:
        vs_coords = LocalCoords.from_coordinates(
            vs = osc,
            atoms = neighb,
            atom_coords = atom_coords
        )
    ff.add_virtual_site(
        charge = oscs_charges[i],
        coords = vs_coords,
        res = ff.residues[0]
    )
for i, (cycle, osc) in enumerate(ring_oscs):
    atom_coords = [conformer.GetAtomPosition(c) for c in cycle]
    vs_coords = LocalCoords.from_coordinates(
        vs = osc,
        atoms = cycle,
        atom_coords = atom_coords
    )
    ff.add_virtual_site(
        charge = oscs_charges[len(atom_oscs) + i],
        coords = vs_coords,
        res = ff.residues[0]
    )

suffixes = []
if args.no_oscs:
    suffixes.append('no')
elif args.no_resp:
    suffixes.append('opls')
else:
    suffixes.append('%s%d%s%s' % (
        args.rsf, args.oscs_count,
        'F' if args.fluorine else '',
        'NR' if args.no_ring else ''
    ))

if args.use_multipoles:
    def make_suffix(weight: float) -> str:
        if weight == 0:
            return '0'
        power = int(np.log10(weight))
        if power >= 0:
            return f'p{power}'
        return f'n{-power}'

    suffixes.append('mult-%s-%s' % (
        make_suffix(args.dipole_weight),
        make_suffix(args.quadrupole_weight)
    ))

elif args.use_bfgs:
    suffixes.append('bfgs')


xml_name = '_'.join([name] + suffixes) + '.xml'
ff.to_xml(xml_name)
