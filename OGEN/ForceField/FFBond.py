import xml.etree.ElementTree as ET
from typing import List, Tuple, Union

from .FFAtom import FFAtom


class FFBond:
    def __init__(
        self,
        indices: Tuple[int, int] = None,
        atoms: Tuple[Union[FFAtom, str], Union[FFAtom, str]] = None
    ):
        self.atom1, self.atom2 = None, None
        self.from_id, self.to_id = None, None
        if atoms is not None:
            if type(atoms[0]) == str:
                self.atom1, self.atom2 = (FFAtom(a) for a in atoms)
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
    def from_xml(element: ET.Element, atoms: List[FFAtom] = None):
        if element.attrib.get('atomName1') is not None:
            id = (element.attrib['atomName1'], element.attrib['atomName2'])
        else:
            id = (int(element.attrib['from']), int(element.attrib['to']))
        if atoms is None:
            if element.attrib.get('atomName1') is not None:
                return FFBond(atoms=id)
            else:
                return FFBond(indices=id)
        if element.attrib.get('atomName1') is not None:
            return FFBond(
                atoms = tuple(
                    a for a in atoms
                        if a.name in [
                            element.attrib['atomName1'],
                            element.attrib['atomName2']
                        ]
                )
            )
        return FFBond(
            atoms = (
                atoms[int(element.attrib['from'])],
                atoms[int(element.attrib['to'])]
            )
        )
