import subprocess
from typing import Dict, List, Tuple, Union
from rdkit.Chem.rdchem import Mol as RDMol
import rdkit.Chem as Chem
import numpy as np
from .boss_prepare import create_zmat, create_params, create_cmd
from .boss_parse import get_bonded, get_dispersion


def reorder_and_shift_atoms(mol: RDMol) -> RDMol:
    heavy_atoms_idx = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() > 1]
    hydrogens_idx = [a.GetIdx() for a in mol.GetAtoms() if a.GetAtomicNum() == 1]
    conformer = mol.GetConformer()
    center_of_mass = sum([
        conformer.GetAtomPosition(a.GetIdx()) * a.GetMass()
            for a in mol.GetAtoms()
    ], start = np.zeros(3)) / sum([a.GetMass() for a in mol.GetAtoms()])
    heavy_atoms_idx.sort(key = lambda i:
        np.linalg.norm(conformer.GetAtomPosition(i) - center_of_mass)
    )
    reordered_atoms = []
    for ai in heavy_atoms_idx:
        if ai in reordered_atoms:
            continue
        for j, aj in enumerate(reordered_atoms):
            if ai in [a.GetIdx() for a in mol.GetAtomWithIdx(aj).GetNeighbors()]:
                reordered_atoms.insert(j + 1, ai)
                break
        else:
            reordered_atoms.append(ai)
    reordered_atoms.extend(hydrogens_idx)
    reordered_mol = Chem.RenumberAtoms(mol, reordered_atoms)
    conformer = reordered_mol.GetConformer()
    for a in reordered_mol.GetAtoms():
        a.SetIntProp('order', reordered_atoms.index(a.GetIdx()))
        conformer.SetAtomPosition(
            a.GetIdx(),
            conformer.GetAtomPosition(a.GetIdx()) - center_of_mass
        )
    return reordered_mol


def calc_opls_parameters(mol: RDMol) -> Tuple[
    List[Dict[str, Union[str, float]]],
    List[Dict[str, Union[Tuple[int, int], float]]],
    List[Dict[str, Union[Tuple[int, int, int], float]]],
    List[Dict[str, Union[Tuple[int, int, int, int], float]]]
]:
    reordered_mol = reorder_and_shift_atoms(mol)
    torsions = create_zmat(reordered_mol)
    create_params('init')
    create_cmd('init')
    subprocess.run('bash bosscmd', check=True, shell=True)
    create_params('internal')
    create_cmd('internal')
    subprocess.run('bash bosscmd', check=True, shell=True)
    create_params('end')
    create_cmd('end')
    subprocess.run('bash bosscmd', check=True, shell=True)
    original_ordering = [a.GetIntProp('order') for a in reordered_mol.GetAtoms()]
    disp = get_dispersion()
    disp[:] = [disp[i] for i in original_ordering]
    bonds_energy, angles_energy, tors_energy = get_bonded(torsions)
    for b in bonds_energy + angles_energy + tors_energy:
        b['atoms'] = tuple(original_ordering.index(i) for i in b['atoms'])
    return disp, bonds_energy, angles_energy, tors_energy
