import itertools
import numpy as np
from OGEN.utils import get_dipole, get_quadrupole, get_octopole
from OGEN.utils.constants import D_TO_AU_ANGSTROM


def kronecker_delta(i: int, j: int) -> int:
    if i == j:
        return 1
    return 0


def test_dipole():
    r = np.array([
        [1, 0, 0],
        [-1, 0, 0],
    ])
    e = np.array([1, -1])

    dipole = get_dipole(r, e)
    real_dipole = np.array([2, 0, 0]) / D_TO_AU_ANGSTROM

    assert (dipole == real_dipole).all()


def test_quadrupole():
    r = np.array([
        [1, 1, 0],
        [1, -1, 0],
        [-1, 1, 0],
        [-1, -1, 0],
    ])
    e = np.array([1, -1, -1, 1])

    quadrupole = get_quadrupole(r, e)
    real_quadrupole = np.zeros((3, 3))
    for i in range(len(e)):
        for alpha, beta in itertools.product(range(3), repeat=2):
            real_quadrupole[alpha][beta] += 0.5 * e[i] * (
                3 * r[i][alpha] * r[i][beta] -
                np.linalg.norm(r[i])**2 * kronecker_delta(alpha, beta))
    real_quadrupole /= D_TO_AU_ANGSTROM

    assert (quadrupole == real_quadrupole).all()


def test_octopole():
    r = np.array([
        [1, 1, 1],
        [1, 1, -1],
        [1, -1, 1],
        [1, -1, -1],
        [-1, 1, 1],
        [-1, 1, -1],
        [-1, -1, 1],
        [-1, -1, -1],
    ])
    e = np.array([1, -1, -1, 1, -1, 1, 1, -1])

    octopole = get_octopole(r, e)
    big_r = np.zeros((3, 3, 3))
    for i in range(len(e)):
        for alpha, beta, gamma in itertools.product(range(3), repeat=3):
            big_r[alpha][beta][gamma] += (e[i] * r[i][alpha] * r[i][beta] *
                                          r[i][gamma])

    real_octopole = np.zeros((3, 3, 3))
    for alpha, beta, gamma in itertools.product(range(3), repeat=3):
        real_octopole[alpha][beta][gamma] += 0.5 * (
            5 * big_r[alpha][beta][gamma] -
            np.trace(big_r[alpha]) * kronecker_delta(beta, gamma) -
            np.trace(big_r[beta]) * kronecker_delta(gamma, alpha) -
            np.trace(big_r[gamma]) * kronecker_delta(alpha, beta))
    real_octopole /= D_TO_AU_ANGSTROM

    assert (octopole == real_octopole).all()
