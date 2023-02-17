import xml.etree.ElementTree as ET
from typing import List

from rdkit.Chem.rdchem import Atom
from rdkit.Chem.rdchem import Mol as RDMol
from rdkit.Chem.rdchem import RWMol as RDMolEdit

from .. import Elements
from .AtomType import AtomType
from .FFAtom import FFAtom
from .FFBond import FFBond
from .VSites import VirtualSite


class Residue:
    def __init__(self,
        name: str,
        atoms: List[FFAtom],
        bonds: List[FFBond],
        vsites: List[VirtualSite] = []
    ):
        self.name = name
        self.atoms = atoms
        self.bonds = []
        for b in bonds:
            if b.from_id is not None:
                b = FFBond(atoms=(self.atoms[b.from_id], self.atoms[b.to_id]))
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
        atoms = [FFAtom.from_xml(atom, types=types) for atom in element.iterfind('Atom')]
        bonds = [FFBond.from_xml(bond, atoms=atoms) for bond in element.iterfind('Bond')]
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

    def to_mol(self) -> RDMol:
        vsites_idx = [v.index for v in self.vsites]
        mol = RDMolEdit()
        atoms = [a for i, a in enumerate(self.atoms) if i not in vsites_idx]
        for a in atoms:
            mol.AddAtom(Atom(Elements[a.type.element].value))
        for b in self.bonds:
            if b.from_id is not None:
                a1 = atoms.index(self.atoms[b.from_id])
                a2 = atoms.index(self.atoms[b.to_id])
            else:
                a1, a2 = atoms.index(b.atom1), atoms.index(b.atom2)
            mol.AddBond(a1, a2)
        return mol.GetMol()

