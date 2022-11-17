import xml.etree.ElementTree as ET
from typing import Dict, List, Tuple


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
