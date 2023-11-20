import os
from typing import List

import numpy as np
from openbabel import pybel

from ..utils import Elements
from .from_wfx import from_wfx


def read_molecule(path: str) -> pybel.Molecule:
    _, extension = os.path.splitext(path)

    if extension == '.wfx':
        return to_obmol(*from_wfx(path))

    return next(pybel.readfile(extension[1:], path))


def to_obmol(symbols: List[Elements], coords: List[np.ndarray],
             charge: int) -> pybel.Molecule:
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
