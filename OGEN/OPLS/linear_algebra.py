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


def find_suitable(coord, cs, n):
    if n > 3:
        raise Exception('Unimplemented n: %d' % n)
    cs = [(i, c) for i, c in enumerate(cs)]
    cs = sorted(cs[3:], key=lambda x: np.linalg.norm(coord - x[1])) + list(reversed(cs[:3]))
    if n > len(cs):
        raise Exception('Not enough suitable coordinates (exp %d, found %d)' % (len(cs), n))
    if n <= 2:
        return [i + 1 for i, _ in cs[:n]]
    cj = cs[0][1]
    for k, ck in cs[1:]:
        if not are_colinear([coord, cj, ck]):
            break
    else:
        raise Exception('Impossible to build Z-matrix')
    for l, cl in cs[1:]:
        if l == k:
            continue
        if not are_colinear([cj, ck, cl]):
            break
    else:
        raise Exception('Impossible to build Z-matrix')
    return cs[0][0] + 1, k + 1, l + 1

