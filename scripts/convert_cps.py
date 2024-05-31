import argparse
import os

from openbabel import openbabel
from OGEN.files.reader import read_molecule


def get_atomic_num_by_point_type(cp_type: int):
    return [2, 10, 18, 36][cp_type - 1]


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('cps', help='file containing critical points')
    parser.add_argument(
        'mol', help='file containing molecule in openbabel-readable format')
    parser.add_argument(
        '--output', help='output file name', default='output.sdf')
    args = parser.parse_args()

    cps: str = args.cps
    mol: str = args.mol
    output: str = args.output

    structure = read_molecule(mol)
    with open(cps) as f:
        for line in f:

            data = line.split()

            cp_type = int(data[0])
            coordinates = [float(c) for c in data[1:]]

            atom = openbabel.OBAtom()

            atom.SetAtomicNum(get_atomic_num_by_point_type(cp_type))
            atom.SetVector(*coordinates)

            structure.OBMol.AddAtom(atom)

    _, file_format = os.path.splitext(output)
    structure.write(file_format[1:], output, overwrite=True)
