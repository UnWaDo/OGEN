from typing import List, Tuple

from .la_selectors import *
from ..Topology import get_connected_indices
from ..MultiWFN.CriticalPoints import CritPoints, SpaceFunctions
from ..MultiWFN.CriticalPoints import generate_atomic_cps, generate_point_cps, generate_xyz
import numpy as np

from ..fchk_parser import Elements


TYPES = {
    'ELF': SpaceFunctions.ELF,
    'LOL': SpaceFunctions.LOL,
    'Lap': SpaceFunctions.Laplasian
}
CHALCOGENS = [Elements.O, Elements.S]
PNICTOGENS = [Elements.N]
HALOGENS = [Elements.F, Elements.Cl, Elements.Br, Elements.I]


def select_atomic_points(
    fchk_path: str,
    atoms: Tuple[List[Elements], List[np.ndarray]],
    atom_i: int,
    mode: str
) -> Tuple[List[CritPoints], List[np.ndarray]]:
    mode_rsf, mode_number, mode_fluorine = mode.split('_')
    mode_rsf = TYPES[mode_rsf]
    mode_number = int(mode_number)
    mode_fluorine = mode_fluorine == 'F'

    elements, coordinates = atoms

    if not mode_fluorine and elements[atom_i] == Elements.F:
        return ([], [])
    connected = get_connected_indices(atoms, atom_i)

    if len(connected) >= 3 and are_coplanar(
        [coordinates[atom_i]] + [coordinates[j] for j in connected]
    ):
        return ([], [])
    mode_cp_type = None
    cp_function = None
    if mode_rsf in [SpaceFunctions.ELF, SpaceFunctions.LOL]:
        if elements[atom_i] in CHALCOGENS:
            if mode_number == 1:
                mode_cp_type = [CritPoints.p3m1, CritPoints.p3p1]
                cp_function = select_one_on_reverse
            elif mode_number == 2:
                mode_cp_type = [CritPoints.p3m3]
                if len(connected) == 2:
                    cp_function = select_two_tetrahedral
                elif len(connected) == 1:
                    cp_function = select_two_on_120
            elif mode_number == 3:
                if len(connected) == 2:
                    mode_cp_type = [CritPoints.p3m1, CritPoints.p3p1]
                    cp_function = select_one_on_reverse
                elif len(connected) == 1:
                    mode_cp_type = [CritPoints.p3m3]
                    cp_function = select_two_on_120
        elif elements[atom_i] in PNICTOGENS:
            mode_cp_type = [CritPoints.p3m3]
            if len(connected) == 2:
                cp_function = select_one_on_reverse
            elif len(connected) == 3:
                cp_function = select_one_on_reverse
        elif elements[atom_i] in HALOGENS:
            mode_cp_type = [CritPoints.p3p1]
            cp_function = select_one_on_reverse
    elif mode_rsf == SpaceFunctions.Laplasian:
        if elements[atom_i] in CHALCOGENS:
            if mode_number == 1:
                cp_function = select_one_on_reverse
                if len(connected) == 2:
                    mode_cp_type = [CritPoints.p3m1]
                elif len(connected) == 1:
                    mode_cp_type = [CritPoints.p3m3]
            elif mode_number == 2:
                mode_cp_type = [CritPoints.p3p1]
                if len(connected) == 2:
                    cp_function = select_two_tetrahedral
                elif len(connected) == 1:
                    cp_function = select_two_on_120
            elif mode_number == 3:
                if len(connected) == 2:
                    mode_cp_type = [CritPoints.p3m1]
                    cp_function = select_one_on_reverse
                elif len(connected) == 1:
                    mode_cp_type = [CritPoints.p3p1]
                    cp_function = select_two_on_120
        elif elements[atom_i] in PNICTOGENS:
            if len(connected) == 2:
                mode_cp_type = [CritPoints.p3m1]
                cp_function = select_one_on_reverse
            elif len(connected) == 3:
                mode_cp_type = [CritPoints.p3p1]
                cp_function = select_one_on_reverse
        elif elements[atom_i] in HALOGENS:
            mode_cp_type = [CritPoints.p3m3]
            cp_function = select_one_on_reverse
    if mode_cp_type is None or cp_function is None:
        raise Exception('Invalid RSF mode: %s (found type? %r. found function? %r)' % (
            mode,
            mode_cp_type is not None,
            cp_function is not None
        ))
    all_cps = generate_atomic_cps(fchk_path, mode_rsf, [atom_i])
    cps = select_points_by_type(all_cps, mode_cp_type)
    cps_coordinates = cp_function(
        coordinates[atom_i], [coordinates[j] for j in connected], cps
    )
    types = [mode_cp_type[0]] * len(cps_coordinates)
    return (types, list(cps_coordinates))


def select_aromatic_points(
    fchk_path: str,
    atoms: Tuple[List[Elements], List[np.ndarray]],
    cycle: List[int],
    mode: str
) -> Tuple[List[CritPoints], List[np.ndarray]]:
    mode_rsf, mode_number, mode_fluorine = mode.split('_')
    mode_rsf = TYPES[mode_rsf]
    mode_number = int(mode_number)
    mode_fluorine = mode_fluorine == 'F'

    elements, coordinates = atoms

    mode_cp_type = CritPoints.p3p1
    center = sum([coordinates[c] for c in cycle], np.zeros(3)) / len(cycle)
    all_cps = generate_point_cps(
        fchk_path,
        mode_rsf,
        center
    )
    # cps = select_points_by_type(all_cps, mode_cp_type)
    cps_coordinates = select_one_at_center(center, all_cps[1])
    types = [mode_cp_type] * len(cps_coordinates)
    return (types, list(cps_coordinates))
