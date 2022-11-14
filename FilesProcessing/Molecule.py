from itertools import combinations, product
from queue import Queue
from typing import List, Tuple, Union
import numpy as np
from . import Atom


class Molecule:
    def __init__(self, atoms: List[Atom], bonds: List[List[int]]):
        '''
            bonds: List[List[int]]
                Represents molecule's bonds. The first atom given is connected to
                    other atoms. Identification is performed based on atom.num
        '''
        self.atoms = sorted(atoms, key=lambda a: a.element)
        self.brutto = {
            e: len([a for a in self.atoms if a.element == e]) 
                for e in set([a.element for a in self.atoms])
        }
        mapping = [a.num for a in self.atoms]
        self.bonds = np.zeros((len(self.atoms),) * 2, dtype=int)
        for bonded in bonds:
            first = mapping.index(bonded[0])
            others = [mapping.index(b) for b in bonded[1:]]
            for o in set(others):
                multiplicity = others.count(o)
                self.bonds[first][o] = multiplicity
                self.bonds[o][first] = multiplicity
        self.path_matrix = None

    @staticmethod
    def create_from_elements(elements: List[str], bonds: List[List[int]]):
        atoms = [Atom(
            element = e,
            type = '%s%d' % (e, i + 1),
            num = i + 1
        ) for i, e in enumerate(elements)]
        bonds = [[i + 1] + [j + 1 for j in bonds[i]] for i in range(len(elements))]
        return Molecule(atoms, bonds)

    def get_neighbors(self, atom: Atom, sort_by_num=False):
        if atom not in self.atoms:
            return None
        neighbours = []
        for i, v in enumerate(self.bonds[self.atoms.index(atom)]):
            if v != 0:
                neighbours.append(self.atoms[i])
        if sort_by_num:
            return sorted(neighbours, key=lambda x: x.num)
        return neighbours

    def create_hetatm_section(self):
        return ''.join(['%s\n' % str(a) for a in sorted(self.atoms, key=lambda x: x.num)])

    def create_conect_section(self):
        conect = []
        for i, a in enumerate(self.atoms):
            total_bonds = sum(self.bonds[i])
            if total_bonds > 4:
                # TODO: make two lines with <= 4 bonds in each
                raise Exception("Unimplemented generationf of CONECT section for n > 4")
            c = 'CONECT%5d' % a.num
            for j, n in enumerate(self.bonds[i]):
                if n == 0:
                    continue
                c += ('%5d' % self.atoms[j].num) * n
            conect.append((a.num, c))
        return ''.join(['%s\n' % c[1] for c in sorted(conect, key=lambda x: x[0])])

    def calc_min_path(self, atom_i: int) -> np.ndarray:
        used = np.zeros(self.bonds.shape[0], dtype=int)
        dist = np.zeros(self.bonds.shape[0], dtype=int)
        used[atom_i] = 1
        q = Queue()
        q.put(atom_i)
        while not q.empty():
            e = q.get()
            for i, x in enumerate(self.bonds[e]):
                if x == 0 or used[i] == 1:
                    continue
                q.put(i)
                used[i] = 1
                dist[i] = dist[e] + 1
        return dist

    def calc_path_matrix(self) -> np.ndarray:
        path_matrix = np.zeros(self.bonds.shape, dtype=int)
        for i in range(path_matrix.shape[0]):
            d = self.calc_min_path(i)
            for j in range(path_matrix.shape[0]):
                path_matrix[i][j] = d[j]
                path_matrix[j][i] = d[j]
        self.path_matrix = path_matrix
        return self.path_matrix

    def __eq__(self, other: 'Molecule'):
        if self.map_atoms(other) is None:
            return False
        return True

    def map_atoms(self, other: 'Molecule') -> Union[None, List[Tuple[Atom, Atom]]]:
        if len(self.brutto) != len(other.brutto):
            return None
        for el in self.brutto:
            if other.brutto.get(el) is None:
                return None
            if self.brutto[el] != other.brutto[el]:
                return None
        if self.path_matrix is None:
            self.calc_path_matrix()
        if other.path_matrix is None:
            other.calc_path_matrix()
        self_conns = self.path_matrix.copy()
        other_conns = other.path_matrix.copy()
        c = 0
        for j in sorted(self.brutto.keys()):
            self_conns[:, c : c + self.brutto[j]].sort()
            other_conns[:, c : c + other.brutto[j]].sort()
            c += self.brutto[j]
        molecules_map = []
        other_ids = [i for i in range(len(other.atoms))]
        for i in range(len(self_conns)):
            for k, j in enumerate(other_ids):
                if np.array_equal(self_conns[i], other_conns[j]):
                    molecules_map.append((self.atoms[i], other.atoms[j]))
                    other_ids.pop(k)
                    break
            else:
                return None
        return molecules_map

    def split_residues(self) -> List['Molecule']:
        if self.path_matrix is None:
            self.calc_path_matrix()
        submolecules = [(
            [self.atoms[0]], [self.bonds[0]]
        )]
        for i in range(1, len(self.atoms)):
            if self.path_matrix[0][i] != 0:
                submolecules[0][0].append(self.atoms[i])
                submolecules[0][1].extend(self.bonds[i])
            else:
                if len(submolecules) == 1:
                    submolecules.append(([], []))
                submolecules[-1][0].append(self.atoms[i])
                submolecules[-1][1].extend(self.bonds[i])
        if len(submolecules) == 1:
            return [self]
        molecules = []
        for submolecule in submolecules:
            molecules.append(Molecule(submolecule[0], submolecule[1]).split_residues())
        molecules = sum(molecules, [])
        for i, molecule in enumerate(molecules):
            for atom in molecule.atoms:
                atom.res = i + 1
        return molecules

    def get_bonds(self) -> List[Tuple[Atom, Atom]]:
        bonds = []
        for a in sorted(self.atoms, key=lambda x: x.num):
            for n in self.get_neighbors(a, True):
                if n.num < a.num:
                    continue
                bonds.append((a, n))
        return bonds

    def get_angles(self) -> List[Tuple[Atom, Atom, Atom]]:
        angles = []
        for atom in self.atoms:
            neighbors = self.get_neighbors(atom, True)
            if len(neighbors) < 2:
                continue
            for a1, a2 in combinations(neighbors, 2):
                angles.append((a1, atom, a2))
        return angles

    def get_torsions(self) -> List[Tuple[Atom, Atom, Atom, Atom, bool]]:
        propers = []
        impropers = []
        for bond in self.get_bonds():
            a2, a3 = bond
            starts = [a for a in self.get_neighbors(a2, True) if a != a3]
            if len(starts) < 1:
                continue
            ends = [a for a in self.get_neighbors(a3, True)
                if a not in starts and a != a2
            ]
            if len(ends) < 1:
                continue
            for a1, a4 in product(starts, ends):
                propers.append((a1, a2, a3, a4, False))
        for a in sorted(self.atoms, key = lambda x: x.num):
            neighbors = self.get_neighbors(a, True)
            if len(neighbors) != 3:
                continue
            a2, a3, a4 = neighbors
            impropers.append((a, a2, a3, a4, True))
            impropers.append((a, a3, a4, a2, True))
            impropers.append((a, a4, a2, a3, True))
        return propers + impropers
