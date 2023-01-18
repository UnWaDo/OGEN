import sys
import os
import numpy as np
from openbabel import pybel


# Van der Waals radii (in angstrom) from 10.1021/jp8111556
VDW_RADII = {'H': 1.10, 'HE': 1.40,
         'LI': 1.81, 'BE': 1.53, 'B': 1.92, 'C': 1.70,
         'N': 1.55, 'O': 1.42, 'F': 1.47, 'NE': 1.54,
         'NA': 2.27, 'MG': 1.73, 'AL': 1.84, 'SI': 2.10,
         'P': 1.80, 'S': 1.80, 'CL': 1.75, 'AR': 1.88, 'BR': 1.83, 'I': 1.98}


def get_intersecting(mol1, mol2):
    for a1 in mol1:
        for a2 in mol2:
            if np.linalg.norm(a1[2] - a2[2]) <= (VDW_RADII[a1[1]] + VDW_RADII[a2[1]]) * SCALE:
                return '%s%d %s%d' % (a1[1], a1[0], a2[1], a2[0])
    return None


FILE = sys.argv[1]
if os.path.isdir(FILE):
    files = [os.path.join(FILE, f) for f in os.listdir(FILE) if f.endswith('.pdb')]
elif not FILE.endswith('.pdb'):
    exit(1)
else:
    files = [FILE]

if len(sys.argv) > 2:
    SCALE = float(sys.argv[2])
else:
    SCALE = 1

for file in files:
    name, _ = os.path.splitext(os.path.basename(file))
    mols = {}
    c = 1
    with open(file, 'r') as f:
        for line in f:
            if not line.startswith('HETATM'):
                continue
            data = line.split()
            coords = np.array([float(d) for d in data[5:8]])
            element = data[-1].upper()
            res = data[4]
            if res not in mols.keys():
                mols[res] = []
            mols[res].append((c, element, coords))
            c += 1
    res = list(mols.keys())
    for i, r1 in enumerate(res):
        for r2 in res[i+1:]:
            intersection = get_intersecting(mols[r1], mols[r2])
            if intersection is not None:
                print('%s, 1, %s' % (name, intersection))
                break
        else:
            continue
        break
    else:
        print('%s, 0' % name)