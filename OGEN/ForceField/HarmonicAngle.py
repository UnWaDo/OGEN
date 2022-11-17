import xml.etree.ElementTree as ET
from typing import List, Tuple


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
    def parse_harmonic_angles(angles: ET.Element) -> List['HarmonicAngle']:
        return [
            HarmonicAngle.from_xml(angle)
                for angle in angles.iterfind('Angle')
        ]
