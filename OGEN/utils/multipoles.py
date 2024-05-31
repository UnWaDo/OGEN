import itertools
import numpy as np
from .constants import D_TO_AU_ANGSTROM


def get_dipole(coords: np.ndarray, charges: np.ndarray) -> np.ndarray:
    return (coords.T @ charges) / D_TO_AU_ANGSTROM


# Definition from 10.1039/QR9591300183
def get_quadrupole(coords: np.ndarray, charges: np.ndarray) -> np.ndarray:
    q = np.einsum('i, ij, ik -> jk', charges, coords, coords)
    theta = 0.5 * (3 * q - np.identity(3) * q.trace())

    return theta / D_TO_AU_ANGSTROM


# Definition from 10.1039/QR9591300183
def get_octopole(coords: np.ndarray, charges: np.ndarray) -> np.ndarray:
    r = np.einsum('i, ij, ik, il -> jkl', charges, coords, coords, coords)
    omega = 5 * r

    for alpha, beta, gamma in itertools.product(range(3), repeat=3):
        omega[alpha][beta][gamma] -= np.trace(
            r[alpha]) * (1 if beta == gamma else 0)

        omega[alpha][beta][gamma] -= np.trace(
            r[beta]) * (1 if gamma == alpha else 0)

        omega[alpha][beta][gamma] -= np.trace(
            r[gamma]) * (1 if alpha == beta else 0)

    return omega / 2 / D_TO_AU_ANGSTROM
