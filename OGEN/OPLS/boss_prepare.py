from typing import List, Tuple
import numpy as np
from .linear_algebra import dist, angle, tors, find_suitable
from ..Topology.topology import get_connected_indices
from ..fchk_parser import Elements
from ...ProcPDB import Molecule

HEADER = 'BOSS Z-matrix by OGEN'
GEOM_VAR = '                    Geometry Variations follow    (2I4,F12.6)'
VAR_BONDS = '                    Variable Bonds follow         (I4)'
ADD_BONDS = '                    Additional Bonds follow       (2I4)'
HARM_CONSTR = '                    Harmonic Constraints follow   (2I4,4F10.4)'
VAR_ANGLES = '                    Variable Bond Angles follow   (I4)'
ADD_ANGLES = '                    Additional Bond Angles follow (3I4)'
VAR_TORS = '                    Variable Dihedrals follow     (3I4,F12.6)'
ADD_TORS = '                    Additional Dihedrals follow   (6I4)'
DOMAIN = '                    Domain Definitions follow     (4I4)'
CONF = '                    Conformational Search (2I4,2F12.6)'
LOC_HEAT = '                    Local Heating Residues follow (I4 or I4-I4)'
FINAL = '                    Final blank line'

DEFAULT_PARNAME = 'parFile.par'


def to_lines(data: List[Tuple[str, int, int, float, int, float, int, float]]) -> List[str]:
    """
    data contains element, atomic number and three z-matrix index with positions
    """
    lines = []
    for i, d in enumerate(data):
        el, t, j, jd, k, ka, l, lt = d
        lines.append(f'{i + 1:4d} {el:3s} {t:4d} {t:4d} {j:4d}{jd:12.6f}{k:4d}{ka:12.6f}{l:4d}{lt:12.6f}')
    return lines


def create_atoms_data(
    atoms: Tuple[List[Elements], List[np.ndarray]]
) -> List[Tuple[str, int, int, float, int, float, int, float]]:
    elements, coordinates = atoms
    x = np.array([0., 0., 0.])
    y = np.array([1., 0., 0.])
    z = np.array([1., 1., 0.])

    act_coords = [x, y, z] + coordinates
    act_elements = ['X', 'X', 'X'] + [e.name for e in elements]
    act_numbers = [-1, -1, -1] + [e.value for e in elements]
    atoms_data = []
    for i, coord in enumerate(act_coords):
        cs = [c for c in act_coords[:i]]
        j, k, l = 0, 0, 0
        jd, ka, lt = 0, 0, 0
        if len(cs) >= 3:
            j, k, l = find_suitable(coord, cs, 3)
        elif len(cs) >= 2:
            j, k = find_suitable(coord, cs, 2)
        elif len(cs) >= 1:
            j, = find_suitable(coord, cs, 1)
        if j:
            jd = dist(coord, act_coords[j-1])
        if k:
            ka = angle(coord, act_coords[j-1], act_coords[k-1])
        if l:
            lt = tors(coord, act_coords[j-1], act_coords[k-1], act_coords[l-1])
        atoms_data.append((
            act_elements[i],
            act_numbers[i],
            j, jd, k, ka, l, lt
        ))
    return atoms_data


def split_variable_and_additional(
    torsions: List[Tuple[int, int, int, int, bool]]
) -> Tuple[List[Tuple[int, int, int, int, bool]], List[Tuple[int, int, int, int, bool]]]:
    varied_atoms = set()
    variable = []
    additional = []
    for torsion in torsions:
        if not torsion[-1] and torsion[0] == max(torsion[:-1]) \
            and not torsion[0] in varied_atoms:
            variable.append(torsion)
            varied_atoms.add(torsion[0])
        else:
            check = tuple(torsion[:3])
            for v in variable:
                if check == tuple(v[:3]):
                    break
            else:
                additional.append(torsion)
    return variable, additional


def create_bonds(mol: Molecule) -> List[str]:
    bonds = [[ba.num for ba in bonded] for bonded in mol.get_bonds()]
    lines = [VAR_BONDS, ADD_BONDS]
    for bd in bonds:
        lines.append(''.join([f'{b + 3:4d}' for b in bd]))
    return lines


def create_angles(mol: Molecule) -> List[str]:
    angles = [[ba.num for ba in bonded] for bonded in mol.get_angles()]
    lines = [VAR_ANGLES, ADD_ANGLES]
    for bd in angles:
        lines.append(''.join([f'{b + 3:4d}' for b in bd]))
    return lines


def create_torsions(mol: Molecule) -> List[str]:
    torsions = [[ba.num for ba in bonded[:-1]] for bonded in mol.get_torsions()]
    lines = [VAR_TORS, ADD_TORS]
    for bd in torsions:
        lines.append(''.join([f'{b + 3:4d}' for b in bd]))
    return lines


def create_zmat(atoms: Tuple[List[Elements], List[np.ndarray]], output: str = 'optzmat'):
    lines = [HEADER]
    elements, coordinates = atoms

    bonded_indices = [get_connected_indices(atoms, i) for i in range(len(elements))]
    mol = Molecule.create_from_elements(elements=[e.name for e in elements], bonds=bonded_indices)
    atoms_data = create_atoms_data(atoms)

    lines += to_lines(atoms_data)
    lines.append(GEOM_VAR)
    lines += create_bonds(mol)
    lines.append(HARM_CONSTR)
    lines += create_angles(mol)
    lines += create_torsions(mol)
    lines += [DOMAIN, CONF, LOC_HEAT, FINAL]
    with open(output, 'w+') as z:
        z.write('\n'.join(lines))


def create_params(stage: str, filename: str = DEFAULT_PARNAME):
    ztype = 'NEWZM' if stage in ['init', 'end'] else 'FULLM'
    with open(filename, 'w+') as f:
        f.write(params.format(QM='AM1SCM1A 1.14  0  0  0', ZTYPE=ztype))


def create_cmd(
    stage: str,
    zmat: str = 'optzmat',
    parname: str = DEFAULT_PARNAME,
    filename: str = 'bosscmd'
):
    mode = '211' if stage in ['init', 'end'] else '711'
    with open(filename, 'w+') as f:
        f.write(cmd.format(ZMAT=zmat, PAR=parname, MODE=mode))


params = '''\
parg49: Gas-phase optimization
 The principal solvent is
NONE 
 QM method, charge calc & scale, solute charges:
{QM}
 The solute format is
ZMAT 
 Origin of the initial solvent coordinates is
NONE 
 Cap atom, radius and force constant:
   0    0.0000    0.0000
 Optimizer, ftol, simulated annealing parameters:
BFGS    0.000100{ZTYPE}         0       0    0.0000    0.0000
 Conformational search parameters:
0           0.0000    0.0000    0.0000    0.0000    0.0000    0.0000
 Parameters for normal mode calculation:
  0  0  0    0.0000    0.0000
 NMOL, IBOX and BOXCUT for principal solvent:
   0   0  0.00
 Second solvent and number of molecules:
NONE    0
 Center atoms of solute1 & 2 and custom solvent, and ICUTAS array:
  0  0  0  0  0  0  0  0
 Atoms solute1 & 2 are rotated about:
  0  0
 Cutoff method (ICUT) and ICUTAT array:
  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 IRECNT, INDSOL, IZLONG, MAXOVL, NOXX, NOSS, NOBNDV, NOANGV, NONEBN:
  0  1  0  0  0  0  0  0  0
 No. of solvent-solvent and solute-solvent rdfs:
  0  0
 First atoms for solvent-solvent rdfs:
  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 Second atoms for solvent-solvent rdfs:
  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0
 First atoms for solute-solvent rdfs:
  0  0  0  0  0  0  0  0  0  0  0  0
 Second atoms for solute-solvent rdfs:
  0  0  0  0  0  0  0  0  0  0  0  0
 Minimum and increment for rdfs:
    0.0000    0.0000
 Min & inc for solvent-solvent pair dist.:
    0.0000    0.0000
 Min & inc for solvent-solvent energy dist.:
    0.0000    0.0000
 Min & inc for solute-solvent pair dist.:
    0.0000    0.0000
 Min & inc for solute-solvent energy dist.:
    0.0000    0.0000
 Frequency of volume & solute moves, and MAXVAR:
     0     0     0
 Freq. of coord output, NBUSE, NOCUT, NOSMTH:
     0     0     0     0
 Range for volume moves, and WKC:
    0.0000    0.0000
 Ranges for solvent translations & rotations, & radius for SASA:
    0.0000    0.0000    0.0000
 Ranges for solute1 translations and rotations:
    0.0000    0.0000
 Ranges for solute2 translations and rotations:
    0.0000    0.0000
 For local heating, solute number and temperature:
   0    0.0000
 Solvent-solvent, solute-solvent, NB cutoffs:
    0.0000    0.0000  100.0000
 Temperature, pressure, & torsion cutoff:
    0.0000    0.0000    0.0000
 Diel. constant, 1-4 Coulomb & LJ scale factors, e(Rxn Field):
    1.0000    2.0000    2.0000    1.0000
 Format for plt files, & solute for E components:
PDB     1
'''

cmd = '''\
boss="$BOSSdir/BOSS"
boxes="$BOSSdir"

export INFILE=optin
export UPFILE=optup
export SAVE=svopt
export AVERAGE=optav
export ZMATRIX="{ZMAT}"
export SLVZMAT=slvzmat
export BANGPAR=$BOSSdir/oplsaa.sb
export WATERBOX=$boxes/watbox
export ORG1BOX=$boxes/org1box
export ORG2BOX=$boxes/org2box
export SUMMARY=sum
export OUTPUT=cout
export PLTFILE=plt.pdb

cp {PAR} tmppar
cat $BOSSdir/oplsaa.par >> tmppar
export PARAMETER=tmppar

# cd run
$boss {MODE} 0 0.0 0.0 0.0
rm slvzmat svopt optav optup optin tmppar plt.pdb
cp sum optzmat
'''
