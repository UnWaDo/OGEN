from typing import List, Tuple

import numpy as np

from ..utils.constants import BOHR_TO_ANGSTROM, Elements


def from_wfx(path: str) -> Tuple[List[Elements], List[np.ndarray], int]:
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

    return elements, coords, charge
