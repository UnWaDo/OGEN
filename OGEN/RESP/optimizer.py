import numpy as np
from scipy import optimize

from ..utils import get_dipole, get_quadrupole


D_TO_CM = 1e-21 / 299792458
AU_ANSGSTROM_TO_CM = 1.6e-19*1e-10
BOHR_TO_ANGSTROM = 0.529177249
D_TO_AU_ANGSTROM = 0.2081943


def fix_charges(charges_no_last):
    total_charge = sum(charges_no_last)
    return np.array([q for q in charges_no_last] + [-total_charge])


def unrestrained_error(charges, options, data):
    charges = fix_charges(charges)
    error = 0

    for inv, esp in zip(data['invr'], data['esp_values']):
        predicted = np.einsum('ij, j -> i', inv, charges)
        error += np.square(predicted - esp).sum()

    if len(data['dipoles']) > 0:
        error += get_dipole_error(charges[:-1], options, data)

    if len(data['quadrupoles']) > 0:
        error += get_quadrupole_error(charges[:-1], options, data)

    return error


def restrained_error(charges, options, data):
    error = unrestrained_error(charges, options, data)
    charges = fix_charges(charges)

    a, b = options['RESP_A'], options['RESP_B']
    for s, q in zip(data['symbols'], charges):
        if s in options['UNRESTRAINED']:
            continue
        error += a * ((q ** 2 + b ** 2) ** 0.5 - b)

    return error


def get_dipole_error(charges, options, data):
    charges = fix_charges(charges)

    error = 0
    for coords, dipole in zip(data['coordinates'], data['dipoles']):
        calc_dipole = get_dipole(coords, charges)

        error += np.linalg.norm(calc_dipole - dipole)

    return error * options['DIPOLES_WEIGHT']


def get_quadrupole_error(charges, options, data):
    charges = fix_charges(charges)

    error = 0
    for coords, quadrupole in zip(data['coordinates'], data['quadrupoles']):
        calc_quadrupole = get_quadrupole(coords, charges)

        error += np.linalg.norm(calc_quadrupole - quadrupole)

    return error * options['QUADRUPOLES_WEIGHT']


def bfgs_solver(options, data):
    res_unrestrained = optimize.minimize(
        unrestrained_error,
        np.zeros(data['natoms'] - 1),
        (options, data),
    )
    if not options['RESTRAINT']:
        return [fix_charges(res_unrestrained.x)], \
               [res_unrestrained.success], [res_unrestrained.message]

    res_restrained = optimize.minimize(
        restrained_error,
        res_unrestrained.x,
        (options, data),
    )
    return [fix_charges(res_unrestrained.x), fix_charges(res_restrained.x)], \
           [res_unrestrained.success, res_restrained.success], \
           [res_unrestrained.message, res_restrained.message]
