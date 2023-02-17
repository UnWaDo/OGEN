import argparse
import os

import numpy as np
from OGEN.ForceField import (AtomType, FFAtom, FFBond, ForceField, Torsion,
                             PeriodicTorsion, HarmonicAngle, HarmonicBond, NonbondedAtom,
                             Residue)
from OGEN.ForceField.VSites import AverageTwo, LocalCoords
from OGEN.MultiWFN import SpaceFunctions, calculate_rsf_at_points
from OGEN.OPLS import calc_opls_parameters
from OGEN.Points.la_selectors import are_colinear
from OGEN.Points.selectors import gen_oscs
from OGEN.RESP import fit_charges, gen_points
from openbabel import pybel
from rdkit import Chem


HELP_FCHK = 'path to .fchk file'
HELP_RSF = 'real-space function for off-site charges generation. Default is ELF'
HELP_OSCS_COUNT = 'amount of OSCs to put on oxygen and sulfur atoms. Default is 3, which puts 2 OSCs on carbonyl oxygen and 1 OSCs elsewhere'
HELP_FLUORINE = 'whether to generate OSCs along the C−F bond. Default is False'
HELP_NO_OSCS = 'toggle this parameter to disable OSCs generation'
HELP_REUSE = 'whether to store files, generated during the procedure. Can speed up recalculations but makes folder a bit messy'

parser = argparse.ArgumentParser()
parser.add_argument('fchk', help=HELP_FCHK)
parser.add_argument('--rsf', default='ELF', help=HELP_RSF, choices=['ELF', 'Lap', 'LOL'])
parser.add_argument('--oscs-count', type=int, default=3, help=HELP_OSCS_COUNT, choices=[1, 2, 3])
parser.add_argument('--fluorine', action='store_true', help=HELP_FLUORINE)
parser.add_argument('--no-oscs', action='store_true', help=HELP_NO_OSCS)
parser.add_argument('--no-reuse', action='store_true', help=HELP_REUSE)
args = parser.parse_args()

fchk_file = os.path.abspath(args.fchk)
pybel_fchk = next(pybel.readfile('fchk', fchk_file))
name, _ = os.path.splitext(os.path.basename(fchk_file))

working_dir = 'OGEN_%s' % name
if not os.path.exists(working_dir):
    os.mkdir(working_dir)
os.chdir(working_dir)

mol_file = '%s.mol' % name
pybel_fchk.write('mol', mol_file, overwrite=True)
mol = Chem.MolFromMolFile(mol_file, removeHs=False, sanitize=False)
conformer = mol.GetConformer()

if args.no_oscs:
    atom_oscs = []
    ring_oscs = []
else:
    atom_oscs, ring_oscs = gen_oscs(
        fchk_file = fchk_file,
        mol = mol,
        mode_rsf = args.rsf,
        mode_number = args.oscs_count,
        mode_fluorine = args.fluorine,
        reuse=True
    )

if args.no_oscs:
    xyz_name = '%s_no.xyz' % name
else:
    xyz_name = '%s_%s%d%s.xyz' % (name, args.rsf, args.oscs_count, 'F' if args.fluorine else '')
with open(xyz_name, 'w+') as xyz:
    xyz.write('%d\n%s\n' % (len(mol.GetAtoms()) + len(atom_oscs) + len(ring_oscs), name))
    confomer = mol.GetConformers()[0]
    for a in mol.GetAtoms():
        coords = confomer.GetAtomPosition(a.GetIdx())
        xyz.write('%2s %20.8e %20.8e %20.8e\n' % (
            a.GetSymbol(),
            coords[0],
            coords[1],
            coords[2]
        ))
    for _, coords in atom_oscs + ring_oscs:
        xyz.write('%2s %20.8e %20.8e %20.8e\n' % (
            'X',
            coords[0],
            coords[1],
            coords[2]
        ))

points_filename = '%s.esp' % name
if os.path.exists(points_filename):
    points = np.loadtxt(points_filename)
else:
    points = gen_points(
        [conformer.GetAtomPosition(i) for i in range(len(mol.GetAtoms()))],
        [a.GetSymbol().upper() for a in mol.GetAtoms()]
    )
    points = calculate_rsf_at_points(
        fchk_file,
        points,
        SpaceFunctions.ESP,
    )
    np.savetxt(points_filename, points)

qf, errors = fit_charges(
    symbols = [a.GetSymbol().upper() for a in mol.GetAtoms()],
    coords = [conformer.GetAtomPosition(i) for i in range(len(mol.GetAtoms()))],
    sample_points = points,
    extra = [o for _, o in atom_oscs + ring_oscs]
)
atom_charges = qf[1][:len(mol.GetAtoms())]
oscs_charges = qf[1][len(mol.GetAtoms()):]

disp, bonds_energy, angles_energy, tors_energy = calc_opls_parameters(mol)

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
    atom_coords = [confomer.GetAtomPosition(n) for n in neighb]
    if len(neighb) == 2 and not are_colinear([osc] + atom_coords):
        neighb.extend([
            a.GetIdx()
                for a in mol.GetAtomWithIdx(neighb[-1]).GetNeighbors()
                    if a.GetIdx() not in neighb
        ])
        atom_coords = [confomer.GetAtomPosition(n) for n in neighb]
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
    atom_coords = [confomer.GetAtomPosition(c) for c in cycle]
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

if args.no_oscs:
    xml_name = '%s_no.xml' % name
else:
    xml_name = '%s_%s%d%s.xml' % (name, args.rsf, args.oscs_count, 'F' if args.fluorine else '')
ff.to_xml(xml_name)
