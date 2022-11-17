import xml.etree.ElementTree as ET
from typing import List

from .AtomType import AtomType


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