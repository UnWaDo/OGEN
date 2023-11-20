import argparse
import os
from typing import Tuple
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


def to_xyz(mol: pybel.Molecule) -> str:
    lines = [f'{len(mol.atoms)}', '']
    for atom in mol:
        lines.append(f'{pybel.ob.GetSymbol(atom.atomicnum):2s}'
                     f'{atom.coords[0]:15.10f}'
                     f'{atom.coords[1]:15.10f}'
                     f'{atom.coords[2]:15.10f}')
    return '\n'.join(lines)


def process_file(file: str):
    name, ext = os.path.splitext(os.path.basename(file))
    mol = next(pybel.readfile(ext[1:], file))

    mol_graph = mol_to_graph(mol)

    for i, connected in enumerate(nx.connected_components(mol_graph)):
        subgraph = nx.subgraph(mol_graph, connected)
        mol = graph_to_mol(subgraph)
        with open(f'{name}_{i + 1}.xyz', 'w') as f:
            f.write(to_xyz(mol))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('files', nargs='+')

    args = parser.parse_args()
    for file in args.files:
        process_file(file)
