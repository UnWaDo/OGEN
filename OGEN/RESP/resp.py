import sys
import warnings
import copy
import numpy as np


from .optimizer import bfgs_solver


bohr_to_angstrom = 0.52917721092


def esp_solve(A, B):
    """Solves for point charges: A*q = B

    Parameters
    ----------
    A : ndarray
        array of matrix A
    B : ndarray
        array of matrix B

    Return
    ------
    q : ndarray
        array of charges

    """

    q = np.linalg.solve(A, B)
    # Warning for near singular matrix
    # in case np.linalg.solve does not detect singularity
    if np.linalg.cond(A) > 1/np.finfo(A.dtype).eps:
        warnings.warn("Possible fit problem; singular matrix")

    return q


def restraint(q, A_unrestrained, resp_a, resp_b, unrestrained_elems, symbols, num_conformers):
    """Adds hyperbolic restraint to matrix A

    Parameters
    ----------
    q : ndarray
        array of charges
    A_unrestrained : ndarray
        array of unrestrained A matrix
    resp_a : float
        restraint scale a
    resp_b : float
        restraint parabola tightness b
    restrained_elems : List[str]
        elements which should be left unrestrained
    symbols : ndarray
        array of element symbols
    num_conformers : int
        the number of conformers

    Returns
    -------
    a : ndarray
        restrained A array
    """

    # hyperbolic Restraint
    # [Bayly:93:10271] (Eqs. 10, 13)
    A = copy.deepcopy(A_unrestrained)
    for i in range(len(symbols)):
        # if an element is not hydrogen or if hydrogens are to be restrained
        if symbols[i] not in unrestrained_elems:
            A[i, i] = A_unrestrained[i, i] + resp_a/np.sqrt(q[i]**2 + resp_b**2) * num_conformers
            # A[i, i] = A_unrestrained[i, i] + resp_a*q[i]/np.sqrt(q[i]**2 + resp_b**2) / 2 * num_conformers

    return A


def iterate(q, A_unrestrained, B, resp_a, resp_b, unrestrained_elems, symbols, toler, maxit, num_conformers):
    """Iterates the RESP fitting procedure

    Parameters
    ----------
    q : ndarray
        array of initial charges
    A_unrestrained : ndarray
        array of unrestrained A matrix
    B : ndarray
        array of matrix B
    resp_a : float
        restraint scale a
    resp_b : float
        restraint parabola tightness b
    restrained_elems : List[str]
        elements which should be left unrestrained
    symbols : ndarray
        array of element symbols
    toler : float
        tolerance for charges in the fitting
    maxit : int
        maximum number of iterations
    num_conformers : int
        the number of conformers

    Returns
    -------
    q : ndarray
        array of the fitted charges

    """
    q_last = copy.deepcopy(q)
    niter, dif, note = 0, 2*toler, ''
    while dif > toler and niter < maxit:
        niter += 1
        A = restraint(q, A_unrestrained, resp_a, resp_b, unrestrained_elems, symbols, num_conformers)
        q = esp_solve(A, B)
        # Extract vector elements that correspond to charges
        dif = np.sqrt(np.max((q[:len(symbols)] - q_last[:len(symbols)])**2))
        q_last = copy.deepcopy(q)

    if dif > toler:
        note += ('\nCharge fitting did not converge; ' + 
               'try increasing the maximum number of iterations to ' +
               '> %i.' %maxit)
    return q[:len(symbols)], note


def intramolecular_constraints(constraint_charge, constraint_groups):
    """Extracts intramolecular constraints from user constraint input

    Parameters
    ----------
    constraint_charge : list
        list of lists of charges and atom indices list
        e.g. [[0, [1, 2]], [1, [3, 4]]]
        The sum of charges on 1 and 2 will equal 0
        The sum of charges on 3 and 4 will equal 1
    constraint_group : list
        list of lists of indices of atoms to have equal charge
        e.g. [[1, 2], [3, 4]]
        atoms 1 and 2 will have equal charge
        atoms 3 and 4 will have equal charge

    Returns
    -------
    constrained_charges : list
        list of fixed charges
    constrained_indices : list
        list of lists of indices of atoms in a constraint
        negative number before an index means
        the charge of that atom will be subtracted.

    Notes
    -----
    Atom indices starts with 1 not 0.
    Total charge constraint is added by default for the first molecule.

    """
    constrained_charges = []
    constrained_indices = []
    for i in constraint_charge:
        constrained_charges.append(i[0])
        group = []
        for k in i[1]:
            group.append(k)
        constrained_indices.append(group)

    for i in constraint_groups:
        for j in range(1, len(i)):
            group = []
            constrained_charges.append(0)
            group.append(-i[j-1])
            group.append(i[j])
            constrained_indices.append(group)

    return constrained_charges, constrained_indices


def fit(options, data):
    """Performs ESP and RESP fits.

    Parameters
    ----------
    options : list
        list of dictionaries of fitting options and internal data

    Returns
    -------
    qf : list
        list of ndarrays of fitted charges
    labelf : list
        list of strings of fitting methods i.e. ESP and RESP
    note : str
        string of notes on the fitting

    """
    qf = []
    labelf = []
    constraint_charges, constraint_indices = intramolecular_constraints(options['CONSTRAINT_CHARGE'],
                                                                        options['CONSTRAINT_GROUP'])
    natoms = data['natoms']
    ndim = natoms + 1 + len(constraint_charges) 
    A = np.zeros((ndim, ndim))
    B = np.zeros(ndim)

    # Bayly:93:10271 (Eqs. 12-14)
    for mol in range(len(data['invr'])):
        r_inverse, V = data['invr'][mol], data['esp_values'][mol]

        # Lower case a and b are the A matrix and B vector for one molecule
        # and without the addition of constraints

        # Construct a: a_jk = sum_i [(1/r_ij)*(1/r_ik)]
        a = np.einsum("ij, ik -> jk", r_inverse, r_inverse)

        # Construct b: b_j = sum_i (V_i/r_ij)
        b = np.einsum('i, ij->j', V, r_inverse)

        # Weight the moleule 
        a *= options['WEIGHT'][mol]**2
        b *= options['WEIGHT'][mol]**2

        A[:natoms, :natoms] += a
        B[:natoms] += b

    # Add total charge constraint
    A[:natoms, natoms] = 1
    A[natoms, :natoms] = 1
    B[natoms] = data['mol_charge']

    # Add constraints to matrices A and B
    for i in range(len(constraint_charges)):
        B[natoms + 1 + i] = constraint_charges[i]
        for k in constraint_indices[i]:
            if k > 0:
                A[natoms + 1 + i, k - 1] = 1
                A[k - 1, natoms + 1 + i] = 1
            else:
                A[natoms + 1 + i, -k - 1] = -1
                A[-k - 1, natoms + 1 + i] = -1

    labelf.append('ESP')
    q = esp_solve(A, B)
    qf.append(q[:natoms])
    if not options['RESTRAINT']:
        return qf, labelf, ''
    else:
        # Restrained ESP
        labelf.append('RESP')
        q, note = iterate(q, A, B, options['RESP_A'], options['RESP_B'], options['UNRESTRAINED'], data['symbols'], options['TOLER'], options['MAX_IT'], len(data['invr']))
        qf.append(q)
        return qf, labelf, note


def fit_charges(symbols, coords, sample_points, extra=[], use_bfgs=False, multipoles={}, multipoles_weights={}):
    symbols = [s.upper() for s in symbols]
    options = {
        'WEIGHT': [1],
        'RESTRAINT': True,
        'RESP_A': 0.0005,
        'RESP_B': 0.1,
        'UNRESTRAINED': ['H', 'N', 'O', 'S', 'F', 'CL', 'BR', 'I', 'X'],
        'IHFREE': True,
        'TOLER': 1e-5,
        'MAX_IT': 25,
        'CONSTRAINT_CHARGE': [],
        'CONSTRAINT_GROUP': [],
        'DIPOLES_WEIGHT': multipoles_weights.get('dipoles', 1),
        'QUADRUPOLES_WEIGHT': multipoles_weights.get('quadrupoles', 1),
    }
    data = {
        'natoms': len(symbols) + len(extra),
        'symbols': symbols + ['X'] * len(extra),
        'mol_charge': 0,
        'dipoles': [],
        'quadrupoles': [],
    }

    if len(multipoles) > 0:
        use_bfgs = True

    esps = sample_points[:,3]
    points = sample_points[:,:-1] * bohr_to_angstrom
    if len(extra):
        coords = np.concatenate((coords, extra))
    data['coordinates'] = [coords]
    data['esp_values'] = [esps]

    invr = np.zeros((len(points), len(coords)))
    for i in range(invr.shape[0]):
        for j in range(invr.shape[1]):
            invr[i, j] = 1/np.linalg.norm(points[i]-coords[j])
    data['invr'] = [invr * bohr_to_angstrom] # convert to atomic units

    dipole = multipoles.get('dipole')
    if dipole is not None:
        data['dipoles'] = [dipole]

    quadrupole = multipoles.get('quadrupole')
    if quadrupole is not None:
        data['quadrupoles'] = [quadrupole]

    if not use_bfgs:
        qf, _, _ = fit(options, data)
    else:
        qf, success, message = bfgs_solver(options, data)
        if not success[0]:
            print('Optimization error without restrain: %s' % message[0],
                  file=sys.stderr)

        elif not success[1]:
            print('Optimization error with restrain: %s' % message[1],
                  file=sys.stderr)

    x_esp = np.einsum('ij, j -> i', data['invr'][0], qf[0])
    x_resp = np.einsum('ij, j -> i', data['invr'][0], qf[1])
    return (qf[0], qf[1]), (np.square(x_esp - esps).sum(), np.square(x_resp - esps).sum())
