import numpy as np
from ..Points.la_selectors import are_colinear


def dist(a1: np.ndarray, a2: np.ndarray) -> float:
    return np.linalg.norm(a1 - a2)


def angle(a1: np.ndarray, a2: np.ndarray, a3: np.ndarray) -> float:
    v1 = a1 - a2
    v2 = a3 - a2
    return np.degrees(np.arccos(np.dot(v1, v2) / np.linalg.norm(v1) / np.linalg.norm(v2)))


def tors(a1: np.ndarray, a2: np.ndarray, a3: np.ndarray, a4: np.ndarray) -> float:
    b1 = -1 * (a2 - a1)
    b2 = a3 - a2
    b3 = a4 - a3
    b2 /= np.linalg.norm(b2)

    v = b1 - np.dot(b1, b2) * b2
    w = b3 - np.dot(b3, b2) * b2

    x = np.dot(v, w)
    y = np.dot(np.cross(b2, v), w)
    return np.degrees(np.arctan2(y, x))


def find_suitable(coord, cs, n, vars):
    if n > 3:
        raise Exception('Unimplemented n: %d' % n)
    atom_n = len(cs)
    if n > atom_n:
        raise Exception('Not enough suitable coordinates (exp %d, found %d)' % (atom_n, n))
    result = []
    for v in vars:
        if atom_n - 3 == v[0]:
            result = [num + 4 for num in v[1:n+1]]
            break
    if len(result) >= n:
        return tuple(result[:n])
    cs = [(i, c) for i, c in enumerate(cs)]
    for i, c in reversed(cs):
        if len(result) < 1 or not are_colinear([coord] + [cs[r - 1][1] for r in result] + [c]):
            result.append(i + 1)
        if len(result) >= n:
            return tuple(result[:n])
    raise Exception('Impossible to build Z-matrix')
