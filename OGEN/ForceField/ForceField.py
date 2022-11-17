import xml.etree.ElementTree as ET
from typing import List

from .AtomType import AtomType
from .FFAtom import FFAtom
from .HarmonicAngle import HarmonicAngle
from .HarmonicBond import HarmonicBond
from .NonbondedAtom import NonbondedAtom
from .PeriodicTorsion import PeriodicTorsion
from .Residue import Residue
from .VSites import VirtualSite, VSCoords


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
        res.atoms.append(FFAtom(
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
