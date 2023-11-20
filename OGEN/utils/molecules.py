import numpy as np
from openbabel import pybel
from openbabel import openbabel as ob


def get_coordinates(molecule: pybel.Molecule):
    return np.array([a.coords for a in molecule])

def get_elements(molecule: pybel.Molecule):
    return list(ob.GetSymbol(a.atomicnum) for a in molecule)
