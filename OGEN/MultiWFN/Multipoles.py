import os
import subprocess
from typing import Dict, List
import numpy as np
from .CriticalPoints import MWFN_EXECUTABLE
from ..utils import D_TO_CM, AU_ANSGSTROM_TO_CM, BOHR_TO_ANGSTROM, D_TO_AU_ANGSTROM

MWFN_OTHER_FUNCTION_3 = '300\n'
MWFN_MULTIPOLES = '5\n'
MWFN_RETURN = '0\n'
MWFN_EXIT = 'q\n'

MULTIPOLES_FILENAME = 'multipole.txt'

HEADER = ' Calculating electric dipole, quadruple, octopole and Hexadecapole moment integral matrix...'
DIPOLE_LINE = ' Dipole moment (a.u.):'
QUADRUPOLE_LINE = ' Quadrupole moments (Traceless Cartesian form):'
OCTOPOLE_LINE = ' Octopole moments (Cartesian form):'
HEXADECAPOLE_LINE = ' Hexadecapole moments:'


def generate_multipoles(fchk_path: str, multipole_file: str = MULTIPOLES_FILENAME):

    with open(multipole_file, 'w') as file:

        subprocess.run([MWFN_EXECUTABLE, fchk_path],
                       input=(MWFN_OTHER_FUNCTION_3 + MWFN_MULTIPOLES +
                              MWFN_RETURN + MWFN_EXIT),
                       text=True,
                       stdout=file)


def get_dipole(path: str) -> np.ndarray:

    with open(path, 'r') as file:

        for line in file:
            if line.startswith(HEADER):
                break

        for line in file:
            if line.startswith(DIPOLE_LINE):
                break
        dipole = np.array([float(c) for c in line.split()[-3:]])

    return dipole * BOHR_TO_ANGSTROM / D_TO_AU_ANGSTROM


def get_quadrupole(path: str) -> np.ndarray:

    with open(path, 'r') as file:

        for line in file:
            if line.startswith(HEADER):
                break

        for line in file:
            if line.startswith(QUADRUPOLE_LINE):
                break
        lines = [file.readline() for i in range(3)]
        elements = [[float(c) for c in line.split()[1::2]] for line in lines]

        quadrupole = np.array(elements)

    return quadrupole * BOHR_TO_ANGSTROM**2 / D_TO_AU_ANGSTROM


def get_multipoles(
        fchk_path: str,
        preserve_files=True,
        multipole_file: str = MULTIPOLES_FILENAME) -> Dict[str, np.ndarray]:

    if not preserve_files or not os.path.exists(multipole_file):
        generate_multipoles(fchk_path, multipole_file)

    multipoles = {
        'dipole': get_dipole(multipole_file),
        'quadrupole': get_quadrupole(multipole_file)
    }

    if not preserve_files:
        os.remove(multipole_file)

    return multipoles
