import xml.etree.ElementTree as ET
from typing import List, Union


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
