import xml.etree.ElementTree as ET
from typing import List, Tuple


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
