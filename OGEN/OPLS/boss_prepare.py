from itertools import combinations, product
from typing import List, Tuple
import numpy as np
from rdkit.Chem.rdchem import Mol as RDMol
from .linear_algebra import dist, angle, tors, find_suitable


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
    mol: RDMol,
    var_bonds: List[int],
    var_angles: List[int],
    var_torsions: List[int]
) -> List[Tuple[str, int, int, float, int, float, int, float]]:
    conformer = mol.GetConformer()
    x = np.array([0., 0., 0.])
    y = np.array([1., 0., 0.])
    z = np.array([1., 1., 0.])

    act_coords = [x, y, z] + [np.array(p) for p in conformer.GetPositions()]
    act_elements = ['X', 'X', 'X'] + [a.GetSymbol() for a in mol.GetAtoms()]
    act_numbers = [-1, -1, -1] + [a.GetAtomicNum() for a in mol.GetAtoms()]
    atoms_data = []
    vars = [t[:-1] for t in var_torsions] + var_angles + var_bonds
    for i, coord in enumerate(act_coords):
        cs = act_coords[:i]
        j, k, l = 0, 0, 0
        jd, ka, lt = 0, 0, 0
        if len(cs) >= 3:
            j, k, l = find_suitable(coord, cs, 3, vars)
        elif len(cs) >= 2:
            j, k = find_suitable(coord, cs, 2, vars)
        elif len(cs) >= 1:
            j, = find_suitable(coord, cs, 1, vars)
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


def create_bonds(mol: RDMol) -> Tuple[List[Tuple[int, int]], List[Tuple[int, int]]]:
    bonds = sorted([tuple(sorted([
        b.GetBeginAtomIdx(), b.GetEndAtomIdx()
    ], reverse=True)) for b in mol.GetBonds()])
    var_bonds = []
    add_bonds = []
    for b in bonds:
        for b_var in var_bonds:
            if b[0] == b_var[0]:
                break
        else:
            var_bonds.append(b)
            continue
        add_bonds.append(b)
    return var_bonds, add_bonds


def create_angles(
    mol: RDMol,
    var_bonds: List[int],
    add_bonds: List[int]
) -> Tuple[List[Tuple[int, int, int]], List[Tuple[int, int, int]]]:
    var_angles = []
    add_angles = []
    conformer = mol.GetConformer()
    for b in var_bonds:
        beginnings = [a.GetIdx() for a in mol.GetAtomWithIdx(b[-1]).GetNeighbors()
            if a.GetIdx() not in b and a.GetIdx() < b[0]
        ]
        for beg in beginnings:
            if np.isclose(180, angle(*[
                np.array(conformer.GetAtomPosition(a)) for a in list(b) + [beg]
            ])):
                add_angles.append((b[0], b[1], beg))
            elif len(var_angles) == 0 or var_angles[-1][0] != b[0]:
                var_angles.append((b[0], b[1], beg))
            else:
                add_angles.append((b[0], b[1], beg))
        endings = [a.GetIdx() for a in mol.GetAtomWithIdx(b[0]).GetNeighbors()
            if a.GetIdx() not in b and a.GetIdx() < b[-1]
        ]
        add_angles.extend([(end, b[0], b[1]) for end in endings])
    for b in add_bonds:
        beginnings = [a.GetIdx() for a in mol.GetAtomWithIdx(b[-1]).GetNeighbors()
            if a.GetIdx() not in b and a.GetIdx() < b[0]
        ]
        add_angles.extend([(b[0], b[1], beg) for beg in beginnings])
        endings = [a.GetIdx() for a in mol.GetAtomWithIdx(b[0]).GetNeighbors()
            if a.GetIdx() not in b and a.GetIdx() < b[-1]
        ]
        add_angles.extend([(end, b[0], b[1]) for end in endings])
    return var_angles, add_angles


def create_torsions(
    mol: RDMol,
    var_angles: List[int],
    add_angles: List[int]
) -> Tuple[List[Tuple[int, int, int, int, bool]], List[Tuple[int, int, int, int, bool]]]:
    var_torsions = []
    add_torsions = []
    conformer = mol.GetConformer()
    for ang in var_angles:
        beginnings = [a.GetIdx() for a in mol.GetAtomWithIdx(ang[-1]).GetNeighbors()
            if a.GetIdx() not in ang and a.GetIdx() < ang[0]
        ]
        for beg in beginnings:
            if np.isclose(180, angle(*[
                np.array(conformer.GetAtomPosition(a)) for a in list(ang[1:]) + [beg]
            ])):
                continue
            if len(var_torsions) == 0 or var_torsions[-1][0] != ang[0]:
                var_torsions.append((ang[0], ang[1], ang[2], beg, False))
            else:
                add_torsions.append((ang[0], ang[1], ang[2], beg, False))
        endings = [a.GetIdx() for a in mol.GetAtomWithIdx(ang[0]).GetNeighbors()
            if a.GetIdx() not in ang and a.GetIdx() < ang[-1]
        ]
        add_torsions.extend([(end, ang[0], ang[1], ang[2], False)
            for end in endings if not np.isclose(180, angle(*[
                np.array(conformer.GetAtomPosition(a)) for a in [end] + list(ang[:-1])
            ]))
        ])
    for ang in add_angles:
        if np.isclose(180, angle(*[np.array(conformer.GetAtomPosition(a)) for a in ang])):
            continue
        beginnings = [a.GetIdx() for a in mol.GetAtomWithIdx(ang[-1]).GetNeighbors()
            if a.GetIdx() not in ang and a.GetIdx() < ang[0]
        ]
        add_torsions.extend([(ang[0], ang[1], ang[2], beg, False)
            for beg in beginnings if not np.isclose(180, angle(*[
                np.array(conformer.GetAtomPosition(a)) for a in list(ang[1:] + [beg])
            ]))
        ])
        endings = [a.GetIdx() for a in mol.GetAtomWithIdx(ang[0]).GetNeighbors()
            if a.GetIdx() not in ang and a.GetIdx() < ang[-1]
        ]
        add_torsions.extend([(end, ang[0], ang[1], ang[2], False)
            for end in endings if not np.isclose(180, angle(*[
                np.array(conformer.GetAtomPosition(a)) for a in [end] + list(ang[:-1])
            ]))
        ])
    for a in mol.GetAtoms():
        neighbors = a.GetNeighbors()
        if len(neighbors) == 3:
            a2, a3, a4 = neighbors
            add_torsions.append((a.GetIdx(), a2.GetIdx(), a3.GetIdx(), a4.GetIdx(), True))
            add_torsions.append((a.GetIdx(), a3.GetIdx(), a4.GetIdx(), a2.GetIdx(), True))
            add_torsions.append((a.GetIdx(), a4.GetIdx(), a2.GetIdx(), a3.GetIdx(), True))
    return var_torsions, add_torsions


def create_zmat(mol: RDMol, output: str = 'optzmat') -> List[Tuple[int, int, int, int]]:
    lines = [HEADER]
    var_bonds, add_bonds = create_bonds(mol)
    var_angles, add_angles = create_angles(mol, var_bonds, add_bonds)
    var_torsions, add_torsions = create_torsions(mol, var_angles, add_angles)
    atoms_data = create_atoms_data(mol, var_bonds, var_angles, var_torsions)
    lines += to_lines(atoms_data)
    lines.extend([GEOM_VAR, VAR_BONDS])
    for bd in var_bonds:
        lines.append(f'{bd[0] + 1 + 3:4d}')
    lines.append(ADD_BONDS)
    for bd in add_bonds:
        lines.append(''.join([f'{b + 1 + 3:4d}' for b in bd]))
    lines.extend([HARM_CONSTR, VAR_ANGLES])
    for ang in var_angles:
        lines.append(f'{ang[0] + 1 + 3:4d}')
    lines.append(ADD_ANGLES)
    for ang in add_angles:
        lines.append(''.join([f'{a + 1 + 3:4d}' for a in ang]))
    lines.append(VAR_TORS)
    for tors in var_torsions:
        lines.append(f'{tors[0] + 1 + 3:4d}{-1:4d}{-1:4d}{0:12.6f}')
    lines.append(ADD_TORS)
    for tors in add_torsions:
        lines.append(''.join([f'{t + 1 + 3:4d}' for t in tors[:-1]]) + f'{-1:4d}{-1:4d}')
    lines += [DOMAIN, CONF, LOC_HEAT, FINAL]
    with open(output, 'w+') as z:
        z.write('\n'.join(lines))
    return var_torsions + add_torsions


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
