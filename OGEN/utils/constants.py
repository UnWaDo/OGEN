from enum import Enum


D_TO_CM = 1e-21 / 299792458
AU_ANSGSTROM_TO_CM = 1.6e-19 * 1e-10
BOHR_TO_ANGSTROM = 0.529177249
D_TO_AU_ANGSTROM = 0.2081943


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


ATOMIC_RADII = {
    Elements.H: 0.32,
    Elements.C: 0.75,
    Elements.N: 0.71,
    Elements.O: 0.63,
    Elements.F: 0.64,
    Elements.S: 1.03,
    Elements.Cl: 0.99,
    Elements.Br: 1.14,
    Elements.I: 1.33,
}
