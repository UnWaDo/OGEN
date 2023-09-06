import os
import subprocess
from typing import Dict, List
import numpy as np
from .CriticalPoints import MWFN_EXECUTABLE
from ..utils import D_TO_CM, AU_ANSGSTROM_TO_CM, BOHR_TO_ANGSTROM, D_TO_AU_ANGSTROM


MWFN_ATOMIC_SPACE = '15\n'
MWFN_MULTIPOLES = '2\n'
MWFN_TO_FILE = '2\n'
MWFN_RETURN = '0\n'
MWFN_EXIT = 'q\n'

MULTIPOLES_FILENAME = 'multipole.txt'


HEADER = '              *****  Molecular dipole and multipole moments  *****'
DIPOLE_LINE = ' Molecular dipole moment (a.u.):'
QUADRUPOLE_LINE = ' Molecular quadrupole moments (Traceless Cartesian form):'
OCTOPOLE_LINE = ' Molecular octopole moments (Cartesian form):'


def generate_multipoles(fchk_path: str):
    subprocess.run([
        MWFN_EXECUTABLE,
        fchk_path
    ], input = (
        MWFN_ATOMIC_SPACE + MWFN_MULTIPOLES + MWFN_TO_FILE + MWFN_RETURN + MWFN_EXIT
    ), text=True, stdout=subprocess.DEVNULL)


def get_dipole(path: str) -> np.ndarray:
    with open(path, 'r') as file:
        line = file.readline()
        while line and not line.startswith(HEADER):
            line = file.readline()

        while line and not line.startswith(DIPOLE_LINE):
            line = file.readline()
        dipole = np.array([float(c) for c in line.split()[-3:]])

    return dipole * BOHR_TO_ANGSTROM / D_TO_AU_ANGSTROM


def get_quadrupole(path: str) -> np.ndarray:
    with open(path, 'r') as file:
        line = file.readline()
        while line and not line.startswith(HEADER):
            line = file.readline()

        while line and not line.startswith(QUADRUPOLE_LINE):
            line = file.readline()
        lines = [file.readline() for i in range(3)]
        elements = [[float(c) for c in line.split()[1::2]] for line in lines]

        quadrupole = np.array(elements)

    return quadrupole * BOHR_TO_ANGSTROM ** 2 / D_TO_AU_ANGSTROM


def get_multipoles(fchk_path: str, preserve_files = True) -> Dict[str, np.ndarray]:
    if not preserve_files or not os.path.exists(MULTIPOLES_FILENAME):
        generate_multipoles(fchk_path)

    multipoles = {
        'dipole': get_dipole(MULTIPOLES_FILENAME),
        'quadrupole': get_quadrupole(MULTIPOLES_FILENAME)
    }

    if not preserve_files:
        os.remove(MULTIPOLES_FILENAME)

    return multipoles
