from enum import Enum
from typing import List
from openbabel import pybel
import numpy as np


BOHR_TO_ANGSTROM = 0.529177


class Elements(Enum):
    H = 1
    C = 6
    N = 7
    O = 8
    F = 9
    S = 16
    Cl = 17
    Br = 35
    I = 53


ATOMIC_RADII = {Elements.H: 0.32, Elements.C: 0.75,
        Elements.N: 0.71, Elements.O: 0.63, Elements.F: 0.64,
        Elements.S: 1.03, Elements.Cl: 0.99, Elements.Br: 1.14, Elements.I: 1.33}


class ParseError(Exception):
    pass



def to_obmol(
        symbols: List[Elements],
        coords: List[np.ndarray],
        charge: int
    ) -> pybel.Molecule:
    mol = pybel.ob.OBMol()

    for i in range(len(symbols)):
        atom = pybel.ob.OBAtom()
        atom.SetAtomicNum(symbols[i].value)
        atom.SetVector(*coords[i])

        mol.AddAtom(atom)
    mol.SetTotalSpinMultiplicity(1)
    mol.SetTotalCharge(charge)
    mol.ConnectTheDots()
    mol.PerceiveBondOrders()

    return pybel.Molecule(mol)
