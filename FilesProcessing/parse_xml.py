import math
from typing import Dict, List, Tuple, Union

import numpy as np
from . import Molecule
from . import Atom as MAtom
import xml.etree.ElementTree as ET


class AtomType:
    def __init__(self,
        name: str,
        class_name: str,
        element: Union[str, None] = None,
        mass: float = 0
    ):
        self.name = name
        self.class_name = class_name
        self.element = element
        self.mass = mass

    def to_xml(self, parent: ET.Element = None) -> ET.Element:
        attrib = {
            'name': self.name,
            'class': self.class_name,
            'mass': str(self.mass)
        }
        if self.element is not None:
            attrib['element'] = self.element
        if parent is None:
            return ET.Element('Type', attrib=attrib)
        return ET.SubElement(parent, 'Type', attrib=attrib)

    @staticmethod
    def from_xml(element: ET.Element):
        return AtomType(
            name = element.attrib['name'],
            class_name = element.attrib['class'],
            element = element.attrib.get('element'),
            mass = element.attrib['mass'],
        )

    @staticmethod
    def create_atom_types(
        atom_types: List['AtomType'],
        parent: ET.Element = None
    ) -> ET.Element:
        if parent is not None:
            atom_types_el = ET.SubElement(parent, 'AtomTypes')
        else:
            atom_types_el = ET.Element('AtomTypes')
        for atom_type in atom_types:
            atom_type.to_xml(atom_types_el)
        return atom_types_el

    @staticmethod
    def parse_atom_types(atom_types: ET.Element) -> List['AtomType']:
        return [
            AtomType.from_xml(atom_type)
                for atom_type in atom_types.iterfind('Type')
        ]

class Atom:
    def __init__(
        self,
        name: str,
        atom_type: AtomType = None,
        type_name: str = '',
    ):
        self.name = name
        if atom_type is None:
            self.type = None
            self._type_name = type_name
        else:
            self.type = atom_type

    def to_xml(self, parent: ET.Element = None) -> ET.Element:
        attrib = {
            'name': self.name,
            'type': self.type_name,
        }
        if parent is None:
            return ET.Element('Atom', attrib=attrib)
        return ET.SubElement(parent, 'Atom', attrib=attrib)

    @property
    def type_name(self) -> str:
        if self.type is None:
            return self._type_name
        return self.type.name

    @staticmethod
    def from_xml(element: ET.Element, types: List[AtomType] = None):
        if types is None:
            return Atom(
                name = element.attrib['name'],
                type_name = element.attrib['type'],
            )
        for t in types:
            if t.name == element.attrib['type']:
                return Atom(
                    name = element.attrib['name'],
                    atom_type = t
                )
        return Atom(
            name = element.attrib['name'],
            type_name = element.attrib['type'],
        )

class Bond:
    def __init__(
        self,
        indices: Tuple[int, int] = None,
        atoms: Tuple[Union[Atom, str], Union[Atom, str]] = None
    ):
        self.atom1, self.atom2 = None, None
        self.from_id, self.to_id = None, None
        if atoms is not None:
            if type(atoms[0]) == str:
                self.atom1, self.atom2 = (Atom(a) for a in atoms)
            else:
                self.atom1, self.atom2 = atoms
        elif indices is not None:
            self.from_id, self.to_id = indices

    def to_xml(self, parent: ET.Element = None) -> ET.Element:
        attrib = {}
        if self.atom1 is not None:
            attrib['atomName1'] = self.atom1.name
            attrib['atomName2'] = self.atom2.name
        elif self.from_id is not None:
            attrib['from'] = str(self.from_id)
            attrib['to'] = str(self.to_id)
        if parent is None:
            return ET.Element('Bond', attrib=attrib)
        return ET.SubElement(parent, 'Bond', attrib=attrib)

    @staticmethod
    def from_xml(element: ET.Element, atoms: List[Atom] = None):
        if element.attrib.get('atomName1') is not None:
            id = (element.attrib['atomName1'], element.attrib['atomName2'])
        else:
            id = (int(element.attrib['from']), int(element.attrib['to']))
        if atoms is None:
            if element.attrib.get('atomName1') is not None:
                return Bond(atoms=id)
            else:
                return Bond(indices=id)
        if element.attrib.get('atomName1') is not None:
            return Bond(
                atoms = tuple(
                    a for a in atoms
                        if a.name in [
                            element.attrib['atomName1'],
                            element.attrib['atomName2']
                        ]
                )
            )
        return Bond(
            atoms = (
                atoms[int(element.attrib['from'])],
                atoms[int(element.attrib['to'])]
            )
        )


class VSCoords:
    atoms: List[int] = []

    def to_dict(self) -> Dict[str, str]:
        pass

    @staticmethod
    def from_xml(element: ET.Element) -> 'VSCoords':
        pass

    @staticmethod
    def from_coordinates(
        vs: np.ndarray,
        atoms: List[int],
        atom_coords: List[np.ndarray] # in Angstrom
    ) -> 'VSCoords':
        pass

    def to_coordinate(self, atom_coords: List[np.ndarray]) -> np.ndarray:
        pass


class LocalCoords(VSCoords):
    def __init__(self,
        atoms: List[int],
        pos: Tuple[float, float, float],
        wo: List[float],
        wx: List[float],
        wy: List[float]
    ):
        self.atoms = atoms
        self.pos = np.array(pos, dtype=float)
        self.wo = np.array(wo, dtype=float)
        self.wx = np.array(wx, dtype=float)
        self.wy = np.array(wy, dtype=float)

    def to_dict(self) -> Dict[str, str]:
        attrib = {}
        attrib.update({'atom%d' % (i + 1): str(a) for i, a in enumerate(self.atoms)})
        attrib.update({'p%d' % (i + 1): str(p) for i, p in enumerate(self.pos)})
        attrib.update({'wo%d' % (i + 1): str(w) for i, w in enumerate(self.wo)})
        attrib.update({'wx%d' % (i + 1): str(w) for i, w in enumerate(self.wx)})
        attrib.update({'wy%d' % (i + 1): str(w) for i, w in enumerate(self.wy)})
        return attrib

    @staticmethod
    def from_xml(element: ET.Element) -> 'LocalCoords':
        atoms_n = len([k for k in element.attrib if k.startswith('atom')])
        atoms = [int(element.attrib['atom%d' % i]) for i in range(1, atoms_n + 1)]
        wo = [float(element.attrib['wo%d' % i]) for i in range(1, atoms_n + 1)]
        wx = [float(element.attrib['wx%d' % i]) for i in range(1, atoms_n + 1)]
        wy = [float(element.attrib['wy%d' % i]) for i in range(1, atoms_n + 1)]
        return LocalCoords(
            atoms = atoms,
            pos = tuple(float(element.attrib['p%d' % i]) for i in range(1, 4)),
            wo = wo,
            wx = wx,
            wy = wy
        )

    @staticmethod
    def from_coordinates(
        vs: np.ndarray,
        atoms: List[int],
        atom_coords: List[np.ndarray] # in Angstrom
    ) -> 'LocalCoords':
        ws = np.zeros((3, len(atoms)))
        ws[0][0] = 1
        ws[1][0], ws[1][1] = -1, 1
        ws[2][0], ws[2][2] = -1, 1
        atom_coords = np.array(atom_coords)

        origin = np.einsum('ij, i -> j', atom_coords, ws[0])
        dirs = np.zeros((3, 3))
        dirs[:, :2] = np.einsum('ij, ki -> jk', atom_coords, ws[1:])
        dirs[:, 2] = np.cross(dirs[:, 0], dirs[:, 1])
        dirs[:, 1] = np.cross(dirs[:, 2], dirs[:, 0])
        dirs /= np.linalg.norm(dirs, axis=0)
        local_coords = np.linalg.solve(dirs, vs - origin) / 10
        return LocalCoords(
            atoms = atoms,
            pos = local_coords,
            wo = ws[0],
            wx = ws[1],
            wy = ws[2]
        )

    def to_coordinate(self, atom_coords: List[np.ndarray]) -> np.ndarray:
        atom_coords = np.array(atom_coords)
        origin = np.einsum('ij, i -> j', atom_coords, self.wo)
        dirs = np.zeros((3, 3))
        dirs[:, :2] = np.einsum('ij, ki -> jk', atom_coords, [self.wx, self.wy])
        dirs[:, 2] = np.cross(dirs[:, 0], dirs[:, 1])
        dirs[:, 1] = np.cross(dirs[:, 2], dirs[:, 0])
        dirs /= np.linalg.norm(dirs, axis=0)
        return origin + np.einsum(
                'ij, j -> i',
                dirs,
                self.pos
        ) * 10

class AverageTwo(VSCoords):
    def __init__(self,
        atoms: Tuple[int, int],
        weights: Tuple[float, float]
    ):
        self.atoms = atoms
        self.weights = np.array(weights, dtype=float)

    def to_dict(self) -> Dict[str, str]:
        attrib = {}
        attrib.update({'atom%d' % (i + 1): str(a) for i, a in enumerate(self.atoms)})
        attrib.update({'weight%d' % (i + 1): str(p) for i, p in enumerate(self.weights)})
        return attrib

    @staticmethod
    def from_xml(element: ET.Element) -> 'AverageTwo':
        return AverageTwo(
            atoms = tuple(int(element.attrib['atom%d' % i]) for i in range(1, 3)),
            weights = tuple(float(element.attrib['weight%d' % i]) for i in range(1, 3))
        )

    @staticmethod
    def from_coordinates(
        vs: np.ndarray,
        atoms: Tuple[int, int],
        atom_coords: Tuple[np.ndarray, np.ndarray]
    ) -> 'AverageTwo':
        v0 = vs - atom_coords[0]
        v1 = vs - atom_coords[1]
        n0 = np.linalg.norm(v0)
        n1 = np.linalg.norm(v1)
        n01 = np.linalg.norm(atom_coords[1] - atom_coords[0])
        metric = np.dot(v0, v1) / n0 / n1
        if not math.isclose(abs(metric), 1, abs_tol=1e-2):
            raise Exception('Atoms and virtual site are not colinear')
        if math.isclose(n01 + n1, n0):
            weights = (-n1 / n01, 1 + n1 / n01)
        elif math.isclose(n01 + n0, n1):
            weights = (1 + n0 / n01, -n0 / n01)
        else:
            weights = (n1 / n01, 1 - n1 / n01)
        return AverageTwo(
            atoms = atoms,
            weights = weights
        )

    def to_coordinate(self, atom_coords: Tuple[np.ndarray, np.ndarray]) -> np.ndarray:
        return atom_coords[0] * self.weights[0] + atom_coords[1] * self.weights[1]


class VirtualSite:
    def __init__(
        self,
        index: int,
        coordinates: VSCoords
    ):
        self.index = int(index)
        self.coordinates = coordinates

    def to_xml(self, parent: ET.Element = None) -> ET.Element:
        attrib = self.coordinates.to_dict()
        attrib['index'] = str(self.index)
        if type(self.coordinates) == LocalCoords:
            attrib['type'] = 'localCoords'
        elif type(self.coordinates) == AverageTwo:
            attrib['type'] = 'average2'
        if parent is None:
            return ET.Element('VirtualSite', attrib=attrib)
        return ET.SubElement(parent, 'VirtualSite', attrib=attrib)

    @staticmethod
    def from_xml(element: ET.Element):
        if element.attrib['type'] == 'localCoords':
            vscoords = LocalCoords
        elif element.attrib['type'] == 'average2':
            vscoords = AverageTwo
        return VirtualSite(
            index = int(element.attrib['index']),
            coordinates = vscoords.from_xml(element),
        )

class Residue:
    def __init__(self,
        name: str,
        atoms: List[Atom],
        bonds: List[Bond],
        vsites: List[VirtualSite] = []
    ):
        self.name = name
        self.atoms = atoms
        self.bonds = []
        for b in bonds:
            if b.from_id is not None:
                b = Bond(atoms=(self.atoms[b.from_id], self.atoms[b.to_id]))
            self.bonds.append(b)
        self.bonds = sorted(
            self.bonds,
            key = lambda b: (
                self.atoms.index(b.atom1),
                self.atoms.index(b.atom2)
            )
        )
        self.vsites = vsites

    def to_xml(self, parent: ET.Element = None) -> ET.Element:
        attrib = {'name': self.name}
        if parent is None:
            element = ET.Element('Residue', attrib=attrib)
        else:
            element = ET.SubElement(parent, 'Residue', attrib=attrib)
        for atom in self.atoms:
            atom.to_xml(element)
        for bond in self.bonds:
            bond.to_xml(element)
        for vsite in self.vsites:
            vsite.to_xml(element)
        return element

    @staticmethod
    def from_xml(element: ET.Element, types: List[AtomType] = None):
        atoms = [Atom.from_xml(atom, types=types) for atom in element.iterfind('Atom')]
        bonds = [Bond.from_xml(bond, atoms=atoms) for bond in element.iterfind('Bond')]
        vsites = [VirtualSite.from_xml(vsite) for vsite in element.iterfind('VirtualSite')]
        return Residue(
            name = element.attrib['name'],
            atoms = atoms,
            bonds = bonds,
            vsites = vsites
        )

    @staticmethod
    def create_residues(
        residues: List['Residue'],
        parent: ET.Element = None
    ) -> ET.Element:
        if parent is not None:
            residues_el = ET.SubElement(parent, 'Residues')
        else:
            residues_el = ET.Element('Residues')
        for residue in residues:
            residue.to_xml(residues_el)
        return residues_el

    @staticmethod
    def parse_residues(residues: ET.Element, types: List[AtomType] = None) -> List['Residue']:
        return [
            Residue.from_xml(residue, types=types)
                for residue in residues.iterfind('Residue')
        ]

    def to_mol(self) -> Molecule:
        vsites_idx = [v.index for v in self.vsites]
        atoms = [a for i, a in enumerate(self.atoms) if i not in vsites_idx]
        bonds = [[i + 1] for i in range(len(atoms))]
        for b in self.bonds:
            if b.from_id is not None:
                a1 = atoms.index(self.atoms[b.from_id])
                a2 = atoms.index(self.atoms[b.to_id])
            else:
                a1, a2 = atoms.index(b.atom1), atoms.index(b.atom2)
            bonds[a1].append(a2 + 1)
            bonds[a2].append(a1 + 1)
        return Molecule(
            atoms = [MAtom(
                element = a.type.element if a.type is not None else 'X',
                type = a.name,
                num = i + 1
            ) for i, a in enumerate(atoms)],
            bonds = bonds
        )


class HarmonicBond:
    def __init__(self,
        length: float,
        k: float,
        classes: Tuple[str, str]
    ):
        self.length = float(length)
        self.k = float(k)
        self.classes = classes

    def to_xml(self, parent: ET.Element = None) -> ET.Element:
        attrib = {
            'length': str(self.length),
            'k': str(self.k),
            'class1': self.classes[0],
            'class2': self.classes[1]
        }
        if parent is None:
            return ET.Element('Bond', attrib=attrib)
        return ET.SubElement(parent, 'Bond', attrib=attrib)

    @staticmethod
    def from_xml(element: ET.Element):
        return HarmonicBond(
            length = element.attrib['length'],
            k = element.attrib['k'],
            classes = (element.attrib['class1'], element.attrib['class2'])
        )

    @staticmethod
    def create_harmonic_bond(
        harmonic_bonds: List['HarmonicBond'],
        parent: ET.Element = None
    ) -> ET.Element:
        if parent is not None:
            harmonic_bonds_el = ET.SubElement(parent, 'HarmonicBondForce')
        else:
            harmonic_bonds_el = ET.Element('HarmonicBondForce')
        for bond in harmonic_bonds:
            bond.to_xml(harmonic_bonds_el)
        return harmonic_bonds_el

    @staticmethod
    def parse_harmonic_bonds(bonds: ET.Element) -> List['HarmonicBond']:
        return [
            HarmonicBond.from_xml(bond)
                for bond in bonds.iterfind('Bond')
        ]

class HarmonicAngle:
    def __init__(self,
        angle: float,
        k: float,
        classes: Tuple[str, str, str]
    ):
        self.angle = float(angle)
        self.k = float(k)
        self.classes = classes

    def to_xml(self, parent: ET.Element = None) -> ET.Element:
        attrib = {
            'angle': str(self.angle),
            'k': str(self.k),
            'class1': self.classes[0],
            'class2': self.classes[1],
            'class3': self.classes[2]
        }
        if parent is None:
            return ET.Element('Angle', attrib=attrib)
        return ET.SubElement(parent, 'Angle', attrib=attrib)

    @staticmethod
    def from_xml(element: ET.Element):
        return HarmonicAngle(
            angle = element.attrib['angle'],
            k = element.attrib['k'],
            classes = tuple(element.attrib['class%d' % i] for i in range(1, 4))
        )

    @staticmethod
    def create_harmonic_angle(
        harmonic_angles: List['HarmonicAngle'],
        parent: ET.Element = None
    ) -> ET.Element:
        if parent is not None:
            harmonic_angles_el = ET.SubElement(parent, 'HarmonicAngleForce')
        else:
            harmonic_angles_el = ET.Element('HarmonicAngleForce')
        for angle in harmonic_angles:
            angle.to_xml(harmonic_angles_el)
        return harmonic_angles_el

    @staticmethod
    def parse_harmonic_angles(angles: ET.Element) -> List['HarmonicBond']:
        return [
            HarmonicAngle.from_xml(angle)
                for angle in angles.iterfind('Angle')
        ]

class Torsion:
    def __init__(self,
        phase: float,
        k: float,
        periodicity: int
    ):
        self.phase = float(phase)
        self.k = float(k)
        self.periodicity = int(periodicity)

    def to_dict(self, i: int = 1) -> Dict[str, str]:
        return {
            'phase%d' % i: str(self.phase),
            'k%d' % i: str(self.k),
            'periodicity%d' % i: str(self.periodicity),
        }

    @staticmethod
    def from_xml(element: ET.Element) -> List['Torsion']:
        torsions = []
        i = 1
        while element.attrib.get('phase%d' % i) is not None:
            torsions.append(Torsion(
                phase = element.attrib['phase%d' % i],
                k = element.attrib['k%d' % i],
                periodicity = element.attrib['periodicity%d' % i],
            ))
            i += 1
        return torsions

class PeriodicTorsion:
    def __init__(self,
        params: List[Torsion],
        classes: Tuple[str, str, str, str],
        is_improper: bool = False
    ):
        self.params = params
        self.classes = classes
        self.is_improper = is_improper

    def to_xml(self, parent: ET.Element = None) -> ET.Element:
        attrib = {
            'class1': self.classes[0],
            'class2': self.classes[1],
            'class3': self.classes[2],
            'class4': self.classes[3]
        }
        for i, torsion in enumerate(self.params):
            attrib.update(torsion.to_dict(i + 1))
        tag = 'Improper' if self.is_improper else 'Proper'
        if parent is None:
            return ET.Element(tag, attrib=attrib)
        return ET.SubElement(parent, tag, attrib=attrib)

    @staticmethod
    def from_xml(element: ET.Element):
        return PeriodicTorsion(
            params = Torsion.from_xml(element),
            classes = tuple(element.attrib['class%d' % i] for i in range(1, 5)),
            is_improper = element.tag == 'Improper'
        )

    @staticmethod
    def create_periodic_torsion(
        periodic_torsions: List['PeriodicTorsion'],
        parent: ET.Element = None
    ) -> ET.Element:
        if parent is not None:
            periodic_torsions_el = ET.SubElement(parent, 'PeriodicTorsionForce')
        else:
            periodic_torsions_el = ET.Element('PeriodicTorsionForce')
        for torsion in periodic_torsions:
            torsion.to_xml(periodic_torsions_el)
        return periodic_torsions_el

    @staticmethod
    def parse_periodic_torsions(torsions: ET.Element) -> List['PeriodicTorsion']:
        return [
            PeriodicTorsion.from_xml(torsion)
                for torsion in torsions.iterfind('Proper')
        ] + [
            PeriodicTorsion.from_xml(torsion)
                for torsion in torsions.iterfind('Improper')
        ]

class NonbondedAtom:
    def __init__(self,
        charge: float,
        sigma: float,
        epsilon: float,
        atom_type: AtomType = None,
        type_name: str = ''
    ):
        self.charge = float(charge)
        self.sigma = float(sigma)
        self.epsilon = float(epsilon)
        if atom_type is None:
            self.type = None
            self._type_name = type_name
        else:
            self.type = atom_type

    def to_xml(self, parent: ET.Element = None) -> ET.Element:
        attrib = {
            'type': self.type_name,
            'charge': str(self.charge),
            'sigma': str(self.sigma),
            'epsilon': str(self.epsilon),
        }
        if parent is None:
            return ET.Element('Atom', attrib=attrib)
        return ET.SubElement(parent, 'Atom', attrib=attrib)

    @property
    def type_name(self) -> str:
        if self.type is None:
            return self._type_name
        return self.type.name

    @staticmethod
    def from_xml(element: ET.Element, types: List[AtomType] = None):
        if types is None:
            return NonbondedAtom(
                type_name = element.attrib['type'],
                charge = element.attrib['charge'],
                sigma = element.attrib['sigma'],
                epsilon = element.attrib['epsilon'],
            )
        for t in types:
            if t.name == element.attrib['type']:
                return NonbondedAtom(
                    atom_type = t,
                    charge = element.attrib['charge'],
                    sigma = element.attrib['sigma'],
                    epsilon = element.attrib['epsilon'],
                )
        return NonbondedAtom(
            type_name = element.attrib['type'],
            charge = element.attrib['charge'],
            sigma = element.attrib['sigma'],
            epsilon = element.attrib['epsilon'],
        )

    @staticmethod
    def create_nonbonded(
        nonbonded_atoms: List['NonbondedAtom'],
        parent: ET.Element = None
    ) -> ET.Element:
        if parent is not None:
            nonbonded_atoms_el = ET.SubElement(parent, 'NonbondedForce')
        else:
            nonbonded_atoms_el = ET.Element('NonbondedForce')
        for atom in nonbonded_atoms:
            atom.to_xml(nonbonded_atoms_el)
        return nonbonded_atoms_el

    @staticmethod
    def parse_nonbonded(nonbonded: ET.Element, types: List[AtomType] = None) -> List['NonbondedAtom']:
        return [
            NonbondedAtom.from_xml(atom, types=types)
                for atom in nonbonded.iterfind('Atom')
        ]

class ForceField:
    atom_types: List[AtomType] = []
    residues: List[Residue] = []
    bond_forces: List[HarmonicBond] = []
    angle_forces: List[HarmonicAngle] = []
    torsion_forces: List[PeriodicTorsion] = []
    nonbonded_forces: List[NonbondedAtom] = []

    torsion_ordering: str = None
    nonbonded_combination: str = None
    coulomb_scale: float = 0.5
    lj_scale: float = 0.5

    @staticmethod
    def from_file(filename: str) -> 'ForceField':
        xml_ff = ET.parse(filename).getroot()
        ff = ForceField()
        ff.atom_types = AtomType.parse_atom_types(xml_ff.find('AtomTypes'))
        ff.residues = Residue.parse_residues(xml_ff.find('Residues'), types=ff.atom_types)
        ff.bond_forces = HarmonicBond.parse_harmonic_bonds(xml_ff.find('HarmonicBondForce'))
        ff.angle_forces = HarmonicAngle.parse_harmonic_angles(xml_ff.find('HarmonicAngleForce'))
        torsion_el = xml_ff.find('PeriodicTorsionForce')
        if torsion_el is None:
            ff.torsion_forces = []
        else:
            ff.torsion_ordering = torsion_el.attrib.get('ordering')
            ff.torsion_forces = PeriodicTorsion.parse_periodic_torsions(torsion_el)
        nonbonded_el = xml_ff.find('NonbondedForce')
        ff.coulomb_scale = float(nonbonded_el.attrib['coulomb14scale'])
        ff.lj_scale = float(nonbonded_el.attrib['lj14scale'])
        ff.nonbonded_combination = nonbonded_el.attrib.get('combination')
        ff.nonbonded_forces = NonbondedAtom.parse_nonbonded(nonbonded_el, types=ff.atom_types)
        return ff

    def to_xml(self, filename: str):
        ff = ET.Element('ForceField')
        AtomType.create_atom_types(self.atom_types, ff)
        Residue.create_residues(self.residues, ff)
        HarmonicBond.create_harmonic_bond(self.bond_forces, ff)
        HarmonicAngle.create_harmonic_angle(self.angle_forces, ff)
        tors = PeriodicTorsion.create_periodic_torsion(self.torsion_forces, ff)
        if self.torsion_ordering is not None:
            tors.attrib['ordering'] = self.torsion_ordering
        nonb = NonbondedAtom.create_nonbonded(self.nonbonded_forces, ff)
        nonb.attrib['coulomb14scale'] = str(self.coulomb_scale)
        nonb.attrib['lj14scale'] = str(self.lj_scale)
        if self.nonbonded_combination is not None:
            nonb.attrib['combination'] = self.nonbonded_combination
        tree = ET.ElementTree(ff)
        ET.indent(tree, '')
        with open(filename, 'wb+') as f:
            tree.write(f)

    def remove_virtual_sites(self):
        unused_types = [at for at in self.atom_types]
        for res in self.residues:
            vsites = [v.index for v in res.vsites]
            res.atoms = [a for i, a in enumerate(res.atoms) if i not in vsites]
            for a in res.atoms:
                try:
                    unused_types.pop(unused_types.index(a.type))
                except ValueError:
                    pass
            res.vsites = []
        self.nonbonded_forces = [nonb for nonb in self.nonbonded_forces
            if nonb.type not in unused_types]
        self.atom_types = [at for at in self.atom_types if at not in unused_types]
    
    def add_virtual_site(self, charge: float, coords: VSCoords, res: Residue):
        at_id = len(self.atom_types)
        at = AtomType(
            name = 'v-site%d' % (at_id),
            class_name = 'X%d' % (at_id),
            element = None,
            mass = 0
        )
        self.atom_types.append(at)
        res.atoms.append(Atom(
            name = 'X%02d' % at_id,
            atom_type = at
        ))
        res.vsites.append(VirtualSite(
            index = len(res.atoms) - 1,
            coordinates = coords
        ))
        self.nonbonded_forces.append(NonbondedAtom(
            charge = charge,
            sigma = 1,
            epsilon = 0,
            atom_type = at
        ))


def merge_xmls(files: List[str], output: Union[str, None] = None) -> None:
    xmls = [(file) for file in files]
