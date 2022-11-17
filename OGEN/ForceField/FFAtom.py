import xml.etree.ElementTree as ET
from typing import List

from .AtomType import AtomType


class FFAtom:
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
            return FFAtom(
                name = element.attrib['name'],
                type_name = element.attrib['type'],
            )
        for t in types:
            if t.name == element.attrib['type']:
                return FFAtom(
                    name = element.attrib['name'],
                    atom_type = t
                )
        return FFAtom(
            name = element.attrib['name'],
            type_name = element.attrib['type'],
        )
