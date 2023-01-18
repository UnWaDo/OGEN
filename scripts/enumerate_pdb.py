import os
import re
import sys


HETATM = re.compile(r'HETATM([\s\d]{5}) ([\s\w]{2})[\s\d]{2} ([A-Z]{3})')


def enumerate_pdb(path: str, output: str):
    count = {}
    with open(output, 'w+') as out:
        with open(path, 'r') as p:
            for line in p:
                match = HETATM.match(line)
                if match is None:
                    out.write(line)
                    continue
                atom_name = match.group(2)
                if count.get(atom_name) is None:
                    count[atom_name] = 0
                count[atom_name] += 1
                out.write(HETATM.sub('HETATM% 5d % 2s%02d %3s' % (
                    int(match.group(1)),
                    atom_name,
                    count[atom_name],
                    match.group(3)
                ), line))


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Not enough arguments given\nUsage: python %s file.pdb or folder with pdbs' % sys.argv[0], file=sys.stderr)
        exit(1)
    path = sys.argv[1]
    if os.path.isdir(path):
        files = [os.path.join(path, f) for f in os.listdir(path) if f.endswith('.pdb')]
    else:
        files = [path]
    for file in files:
        out = os.path.join('form', file)
        if not os.path.exists('form'):
            os.makedirs('form')
        enumerate_pdb(file, out)
