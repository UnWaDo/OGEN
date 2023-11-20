import numpy as np
import pytest
from openbabel import pybel
from OGEN.openmm_api import mol_to_openmm


@pytest.fixture
def ethanol() -> pybel.Molecule:
    mol = pybel.readstring('smi', 'CCO')
    mol.OBMol.AddHydrogens()

    return mol


def test_to_openmm(ethanol: pybel.Molecule):
    omm = mol_to_openmm(ethanol)
    residues = list(omm.residues())

    assert len(residues) == 1
    assert omm.getNumAtoms() == ethanol.OBMol.NumAtoms()
    assert omm.getNumBonds() == ethanol.OBMol.NumBonds()
