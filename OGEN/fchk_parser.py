from enum import Enum
from typing import List, TextIO, Tuple

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


def get_fchk_values(file: TextIO, header: str, dtype: type = float) -> np.ndarray:
    for line in file:
        if header in line:
            break
    else:
        return None
    values = []
    for line in file:
        nums = line.split()
        try:
            dtype(nums[0])
        except ValueError:
            break
        for n in nums:
            values.append(dtype(n))
    return np.array(values, dtype)


def get_coordinates(fchk_path: str) -> List[np.ndarray]:
    '''Reads .fchk file and returns coordinates in Angstrom'''
    with open(fchk_path) as file:
        coordinates = get_fchk_values(file, 'Current cartesian coordinates', float)
        if coordinates is None:
            raise ParseError('Coordinates are not located')
    return list(coordinates.reshape((-1, 3)) * BOHR_TO_ANGSTROM)


def get_elements(fchk_path: str) -> List[Elements]:
    with open(fchk_path) as file:
        atoms = get_fchk_values(file, 'Atomic numbers', int)
        if atoms is None:
            raise ParseError('Atom numbers are not located')
    return [Elements(a) for a in atoms]


def get_masses(fchk_path: str) -> List[float]:
    with open(fchk_path) as file:
        masses = get_fchk_values(file, 'Real atomic weights')
        if masses is None:
            raise ParseError('Atom weights are not located')
    return masses


def get_atoms(fchk_path: str) -> Tuple[List[Elements], List[np.ndarray]]:
    return get_elements(fchk_path), get_coordinates(fchk_path)
