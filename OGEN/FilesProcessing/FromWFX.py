import numpy as np
from openbabel import pybel

from .Processer import ParseError, Elements, to_obmol, BOHR_TO_ANGSTROM


def from_wfx(path: str) -> pybel.Molecule:
    elements = []
    coords = []
    charge = 0
    with open(path) as file:
        line = file.readline()
        while line and line != '<Number of Nuclei>\n':
            line = file.readline()

        if not line:
            return None

        atoms_number = int(file.readline())

        line = file.readline()
        while line and line != '<Net Charge>\n':
            line = file.readline()

        if not line:
            return None

        charge = int(file.readline())

        line = file.readline()
        while line and line != '<Atomic Numbers>\n':
            line = file.readline()

        if not line:
            return None

        for i in range(atoms_number):
            elements.append(Elements(int(file.readline())))

        line = file.readline()
        while line and line != '<Nuclear Cartesian Coordinates>\n':
            line = file.readline()

        if not line:
            return None

        for i in range(atoms_number):
            coords.append(np.array([
                float(t) for t in file.readline().split()
            ]) * BOHR_TO_ANGSTROM)

    return to_obmol(elements, coords, charge)
