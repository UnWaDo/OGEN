import argparse
import os
from typing import List, Tuple
from openbabel import pybel
from rdkit import Chem
from rdkit.Chem.rdchem import Mol as RDMol
import numpy as np

from OGEN.Points.selectors import select_atomic_points, select_aromatic_points


def gen_oscs(
    fchk_file: str,
    mol: RDMol,
    mode: str,
    reuse = False
) -> List[Tuple[int, np.ndarray]]:
    atom_oscs = []
    for a in mol.GetAtoms():
        if a.GetAtomicNum() <= 6:
            continue
        a_points = select_atomic_points(
            fchk_path = fchk_file,
            mol = mol,
            atom = a,
            mode = mode,
            reuse = reuse
        )
        atom_oscs.extend([(a.GetIdx(), ap) for ap in a_points])
    ring_oscs = []
    for cycle in Chem.GetSSSR(mol):
        crit_point = select_aromatic_points(fchk_file, mol, cycle, mode, reuse)
        ring_oscs.append(([c for c in cycle], crit_point[0]))
    return atom_oscs, ring_oscs


parser = argparse.ArgumentParser()
parser.add_argument('fchk', help='path to .fchk file')
parser.add_argument('--no-oscs', action='store_false', help='whether to generate OSCs')
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

atom_oscs, ring_oscs = gen_oscs(fchk_file, mol, 'ELF_2_X', reuse=True)
with open('%s.xyz' % name, 'w+') as xyz:
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


# for b in mol.GetBonds():
#     print(b.GetBeginAtomIdx(), b.GetEndAtomIdx())
# if args.mode != 'no_oscs':
#     atom_oscs, ring_oscs = gen_oscs(args.fchk, atoms, args.mode)
# else:
#     atom_oscs = []
#     ring_oscs = []
# mol2 = Chem.MolFromSmiles('OC(=O)C')

# print([i for i in range(mol.GetNumAtoms())])
# print(mol.GetSubstructMatches(mol2))

# print(Chem.MolToMolBlock(mol))
# for f in add_files:
#     os.remove(f)