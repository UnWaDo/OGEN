from typing import List
from rdkit.Chem.rdchem import Mol as RDMol

from ..Points import AtomOSC, RingOSC


def to_xyz(path: str,
           mol: RDMol,
           atom_oscs: List[AtomOSC] = [],
           ring_oscs: List[RingOSC] = [],
           name: str = ''):

    with open(path, 'w+') as xyz:
        atoms_number = len(mol.GetAtoms()) + len(atom_oscs) + len(ring_oscs)

        xyz.write(f'{atoms_number}\n{name}\n')

        confomer = mol.GetConformers()[0]
        for a in mol.GetAtoms():
            coords = confomer.GetAtomPosition(a.GetIdx())
            xyz.write(f'{a.GetSymbol():2s} '
                      f'{coords[0]:20.8e} '
                      f'{coords[1]:20.8e} '
                      f'{coords[2]:20.8e}\n')

        for _, coords in atom_oscs + ring_oscs:
            xyz.write('X  '
                      f'{coords[0]:20.8e} '
                      f'{coords[1]:20.8e} '
                      f'{coords[2]:20.8e}\n')
