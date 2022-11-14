from typing import Dict, List, Tuple
import numpy as np
from ..fchk_parser import ATOMIC_RADII, Elements


def get_connected_indices(
    atoms: Tuple[List[Elements],List[np.ndarray]],
    atom_i: int
) -> List[int]:
    elements, coordinates = atoms
    i_coords = coordinates[atom_i]
    i_atomic_r = ATOMIC_RADII[elements[atom_i]]
    connected = []
    for atom_j, j_coords in enumerate(coordinates):
        if atom_j == atom_i:
            continue
        j_atomic_r = ATOMIC_RADII[elements[atom_j]]
        if np.linalg.norm(i_coords - j_coords) < (i_atomic_r + j_atomic_r) * 1.2:
            connected.append(atom_j)
    return connected


def remove_side_chains(connected: Dict[int, List[int]]) -> Dict[int, List[int]]:
    one_bond = [i for i in connected if len(connected[i]) == 1]
    while len(one_bond):
        for i in one_bond:
            if len(connected[i]) == 0:
                connected.pop(i)
                continue
            j = connected[i][0]
            connected[j].pop(connected[j].index(i))
            connected.pop(i)
        one_bond = [i for i in connected if len(connected[i]) == 1]
    return connected


def find_cycles(connected: Dict[int, List[int]]):
    visited = np.zeros(len(connected), int)
    map_indices = [i for i in connected]
    cycles = set()
    def find_cycle(i: int, cycle: List[int]):
        visited[map_indices.index(i)] = 1
        for j in connected[i]:
            if cycle[0] == j and len(cycle) > 2:
                cycles.add(tuple(sorted(cycle)))
            elif visited[map_indices.index(j)] == 0 and j > cycle[0]:
                cycle.append(j)
                find_cycle(j, cycle)
                cycle.pop()
        visited[map_indices.index(i)] = 0
    for i in connected:
        cycle = [i]
        find_cycle(i, cycle)
    return cycles
