import os
import sys
import networkx as nx
import re


HETATM = re.compile(r'HETATM([\s\d]{5}) ([\s\w]{2})([\s\d]{2})(.)([\w\d]{3}) (.)([\s\d]{4})(.)   ([\s\+\-\d\.]{8})([\s\+\-\d\.]{8})([\s\+\-\d\.]{8})([\s\+\-\d\.]{6})([\s\+\-\d\.]{6})\s{10}([\s\w]{2})([\s\w\+\-\.\d]{2})')
HETATM_TYPES = {
    'serial': int, 'name_el': str, 'name_num': int,
    'altLoc': str, 'resName': str, 'chainID': str,
    'resSeq': int, 'iCode': str,
    'x': float, 'y': float, 'z': float,
    'occupancy': float, 'tempFactor': float,
    'element': str, 'charge': str
}


def cast_type(string, dt):
    if dt == str:
        return string.strip()
    try:
        return dt(string)
    except ValueError as e:
        default = {int: 0, float: 0}
        return default[dt]


def parse_hetatm(line: str):
    match = HETATM.match(line)
    if match is None:
        return None
    return {k: cast_type(match.group(i + 1), HETATM_TYPES[k]) for i, k in enumerate(HETATM_TYPES)}


def to_hetatm(hetatm_dict):
    return 'HETATM%5d %2s%-2d%1s%3s %1s%4d%1s   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s' % tuple(
        hetatm_dict[name] for name in HETATM_TYPES
    )

FILE = sys.argv[1]
if os.path.isdir(FILE):
    files = [os.path.join(FILE, f) for f in os.listdir(FILE) if f.endswith('.pdb')]
else:
    files = [FILE]

for file in files:
    with open(file, 'r') as f:
        lines = f.readlines()
    atoms = [parse_hetatm(l) for l in lines if l.startswith('HETATM')]
    conects = [[int(i) for i in l.split()[1:]] for l in lines if l.startswith('CONECT')]
    G = nx.Graph()
    G.add_nodes_from([a['serial'] for a in atoms])
    for connected in conects:
        G.add_edges_from([(connected[0], c) for c in connected[1:]])
    residues = [g for g in nx.connected_components(G)]
    if len(residues) != 2:
        print('Are you sure you using correct PDB file? The number of residues is not 2', file=sys.stderr)
    new_atoms = sorted(atoms, key=lambda a: a['serial'] not in residues[0])
    match_serials = [a['serial'] for a in new_atoms]
    for i, line in enumerate(lines):
        if line.startswith('HETATM'):
            break
    output = lines[:i]
    for new_i, a in enumerate(new_atoms):
        a['resSeq'] = 1 if a['serial'] in residues[0] else 2
        a['name_num'] = new_i
        a['name_el'] = a['name_el'].title()
        a['serial'] = new_i + 1
        output.append(f'{to_hetatm(a)}\n')
    new_conects = []
    for connected in conects:
        new_conects.append([match_serials.index(c) + 1 for c in connected])
    new_conects.sort(key=lambda x: x[0])
    for connected in new_conects:
        output.append('CONECT%s\n' % ''.join(['%5d' % c for c in connected]))
    output.extend(lines[i + len(atoms) + len(conects):])
    with open(file, 'w') as f:
        f.writelines(output)

