import math
import xml.etree.ElementTree as ET
from typing import Dict, List, Tuple

import numpy as np


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
