import argparse
import os
from typing import List, Tuple
import networkx as nx
import numpy as np
from openbabel import pybel

ATOM = Tuple[int, int, Tuple[float, float, float]]


def atom_to_tuple(atom: pybel.ob.OBAtom) -> ATOM:
    return atom.GetId(), atom.GetAtomicNum(), (atom.GetX(), atom.GetY(),
                                               atom.GetZ())


def mol_to_graph(molecule: pybel.Molecule) -> nx.Graph:
    mol_graph = nx.Graph()
    atoms = [atom_to_tuple(a.OBAtom) for a in molecule.atoms]

    mol_graph.add_nodes_from(atoms)
    for i in range(molecule.OBMol.NumBonds()):
        b = molecule.OBMol.GetBond(i)
        mol_graph.add_edge(atom_to_tuple(b.GetBeginAtom()),
                           atom_to_tuple(b.GetEndAtom()))

    return mol_graph


def graph_to_mol(graph: nx.Graph, charge=0) -> pybel.Molecule:
    mol = pybel.ob.OBMol()
    for node in graph:
        atom = pybel.ob.OBAtom()
        atom.SetAtomicNum(node[1])
        atom.SetVector(*node[2])

        mol.AddAtom(atom)
    mol.ConnectTheDots()
    mol.SetTotalSpinMultiplicity(1)
    mol.SetTotalCharge(charge)
    mol.PerceiveBondOrders()

    return pybel.Molecule(mol)


def to_xyz(mol: pybel.Molecule, comment: str = '') -> str:
    lines = [f'{len(mol.atoms)}\n', f'{comment}\n']

    for atom in mol:
        lines.append(f'{pybel.ob.GetSymbol(atom.atomicnum):2s}'
                     f'{atom.coords[0]:15.10f}'
                     f'{atom.coords[1]:15.10f}'
                     f'{atom.coords[2]:15.10f}\n')

    return ''.join(lines)


def get_coordinates(molecule: pybel.Molecule) -> np.ndarray:
    return np.array([a.coords for a in molecule.atoms])


def find_closest(coordinates: Tuple[np.ndarray, np.ndarray]) -> np.ndarray:
    a, b = coordinates
    distances = np.linalg.norm(a[:, None, :] - b[None, :, :], axis=-1)

    i = np.argmin(distances)
    j = i % len(distances)
    i //= len(distances)
    return [coordinates[0][i], coordinates[1][j]]


def find_centroid(coordinates: np.ndarray) -> np.ndarray:
    return coordinates.mean(axis=1)


def set_molecule_coordinates(mol: pybel.Molecule, coordinates: np.ndarray) -> None:
    for a, coordinate in zip(mol.atoms, coordinates):
        a.OBAtom.SetVector(*coordinate)

def process_file(path: str, use_closest = True):
    name, ext = os.path.splitext(os.path.basename(path))
    mol = next(pybel.readfile(ext[1:], path))

    mol_graph = mol_to_graph(mol)

    parts: List[np.ndarray] = []
    molecules: List[pybel.Molecule] = []

    for i, connected in enumerate(nx.connected_components(mol_graph)):
        subgraph = nx.subgraph(mol_graph, connected)
        mol = graph_to_mol(subgraph)

        parts.append(get_coordinates(mol))
        molecules.append(mol)

    if use_closest:
        ref = find_closest(parts)
    else:
        ref = [find_centroid(c) for c in parts]

    assert len(ref) == 2, 'must use dimers'

    line = ref[1] - ref[0]
    distance = np.linalg.norm(line)

    with open('test.xyz', 'a') as file:
        molecule = molecules[0].clone

        # set_molecule_coordinates(molecules[1], parts[1] + line * (distance - 1))
        set_molecule_coordinates(molecules[1], parts[1])

        molecule.OBMol += molecules[1].OBMol

        file.write(to_xyz(molecule, f'{distance} {name[-3:]}'))

        # centroid = find_centroid(coordinates)

        # print(centroid - coordinates)
        # with open(f'{name}_{i + 1}.xyz', 'w') as f:
        #     f.write(to_xyz(mol))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('files', nargs='+')

    args = parser.parse_args()

    with open('test.xyz', 'w') as file:
        pass

    for file in args.files:
        process_file(file)
