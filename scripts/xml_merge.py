import argparse
from copy import copy
import os
import sys

from OGEN.ForceField import FFAtom, AtomType, ForceField, NonbondedAtom, VirtualSite, Residue


def get_types(res: Residue):
    return [atom.type for atom in res.atoms]


def get_types_no_oscs(res: Residue):
    vsites_idx = [v.index for v in res.vsites]
    return [atom.type for i, atom in enumerate(res.atoms) if i not in vsites_idx]


def copy_disp(ff: ForceField, ff_other: ForceField):
    mol = ff.residues[0].to_mol()
    mol_other = ff_other.residues[0].to_mol()
    try:
        atoms_map, = mol.GetSubstructMatches(mol_other)
    except ValueError:
        raise Exception('Dispersion file contain different molecule')
    if len(atoms_map) != mol.GetNumAtoms() or len(atoms_map) != mol_other.GetNumAtoms():
        raise Exception('Dispersion file contain different molecule')
    atom_types1 = get_types_no_oscs(ff.residues[0])
    atom_types2 = get_types_no_oscs(ff_other.residues[0])
    atom_types1 = [atom_types1[i] for i in atoms_map]

    for nb_other in ff_other.nonbonded_forces:
        try:
            at = atom_types1[atom_types2.index(nb_other.type)]
        except ValueError:
            continue
        for nb in ff.nonbonded_forces:
            if nb.type == at:
                nb.epsilon = nb_other.epsilon
                nb.sigma = nb_other.sigma
                break
    ff.nonbonded_combination = ff_other.nonbonded_combination
    ff.lj_scale = ff_other.lj_scale

def copy_charges(ff: ForceField, ff_other: ForceField):
    mol = ff.residues[0].to_mol()
    mol_other = ff_other.residues[0].to_mol()
    try:
        atoms_map, = mol.GetSubstructMatches(mol_other)
    except ValueError:
        raise Exception('Dispersion file contain different molecule')
    if len(atoms_map) != mol.GetNumAtoms() or len(atoms_map) != mol_other.GetNumAtoms():
        raise Exception('Dispersion file contain different molecule')
    atom_types1 = get_types_no_oscs(ff.residues[0])
    atom_types2 = get_types_no_oscs(ff_other.residues[0])
    atom_types1 = [atom_types1[i] for i in atoms_map]

    for i, vs in enumerate(ff_other.residues[0].vsites):
        new_type = AtomType(name='osc%02d' % (i + 1), class_name='osc%02d_cl' % (i + 1))
        vs_atom = ff_other.residues[0].atoms[vs.index]
        ff.atom_types.append(new_type)
        atom_types1.append(new_type)
        atom_types2.append(vs_atom.type)
        ff.residues[0].atoms.append(FFAtom(name='X%02d' % (i+1), atom_type=new_type))
        coordinates = copy(vs.coordinates)
        coordinates.atoms = []
        for a in vs.coordinates.atoms:
            other_a = ff_other.residues[0].atoms[a]
            idx = atom_types2.index(other_a.type)
            for j, ff_a in enumerate(ff.residues[0].atoms):
                if ff_a.type == atom_types1[idx]:
                    break
            else:
                Exception()
            coordinates.atoms.append(j)
        ff.residues[0].vsites.append(VirtualSite(
            index = len(ff.residues[0].atoms) - 1,
            coordinates = coordinates
        ))
    for nb_other in ff_other.nonbonded_forces:
        at = atom_types1[atom_types2.index(nb_other.type)]
        for nb in ff.nonbonded_forces:
            if nb.type == at:
                nb.charge = nb_other.charge
                break
        else:
            ff.nonbonded_forces.append(NonbondedAtom(
                charge = nb_other.charge,
                epsilon = nb_other.epsilon,
                sigma = nb_other.sigma,
                atom_type = at
            ))
    ff.coulomb_scale = ff_other.coulomb_scale


OUTPUT_DEFAULT = 'output'
parser = argparse.ArgumentParser()
parser.add_argument('bonded', help='path to .xml forcefield file with bonded parameters')
parser.add_argument('disp', help='path to .xml forcefield file with dispersion parameters')
parser.add_argument('charges', help='path to .xml forcefield file with charges')
parser.add_argument('--output', '-o', required=False, help='path to resulting file. default is output/*xml*')

args = parser.parse_args()
if not os.path.exists(args.bonded):
    print('Input file %s do not exist' % args.bonded, file=sys.stderr)
    exit(1)
if not os.path.exists(args.disp):
    print('Input file %s do not exist' % args.disp, file=sys.stderr)
    exit(1)
if not os.path.exists(args.charges):
    print('Input file %s do not exist' % args.charges, file=sys.stderr)
    exit(1)
if args.output is None:
    args.output = os.path.join(OUTPUT_DEFAULT, os.path.basename(args.bonded))
else:
    path, ext = os.path.splitext(args.output)
    if ext != '.xml':
        args.output = os.path.join(args.output, os.path.basename(args.bonded))
output_dir = os.path.dirname(args.output)
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

name, _ = os.path.splitext(os.path.basename(args.bonded))
ff = ForceField.from_file(args.bonded)
ff.remove_virtual_sites()

copy_disp(ff, ForceField.from_file(args.disp))
copy_charges(ff, ForceField.from_file(args.charges))
ff.to_xml(args.output)
