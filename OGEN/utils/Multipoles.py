import numpy as np
from .Constants import D_TO_AU_ANGSTROM


def get_dipole(coords: np.ndarray, charges: np.ndarray) -> np.ndarray:
    return (coords.T @ charges) / D_TO_AU_ANGSTROM


# Definition from 10.1039/QR9591300183
def get_quadrupole(coords: np.ndarray, charges: np.ndarray) -> np.ndarray:
    q = np.einsum('i, ij, ik -> jk', charges, coords, coords)
    theta = 0.5 * (3 * q - np.identity(3) * q.trace())

    return theta / D_TO_AU_ANGSTROM


# Definition from 10.1039/QR9591300183
# def get_octopole(coords: np.ndarray, charges: np.ndarray) -> np.ndarray:
#     r = np.einsum('i, ij, ik, il -> jkl', charges, coords, coords, coords)
#     omega = 0.5 * (5 * r - np.identity(3, 3))
#     return theta * AU_ANSGSTROM_TO_CM / D_TO_CM
