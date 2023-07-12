import re
from typing import List, Tuple
from . import Atom, Molecule


hetatm = re.compile(r'HETATM([ \d]{5}) ([ \w]{2}[ \d]{2}).([ \w]{3}) .([ \d]{4}). {3}([\.\- \d]{8})([\.\- \d]{8})([\.\- \d]{8}).{22}([ \w]{2}).*')
conect = re.compile(r'CONECT[ \d]*')


def parse_pdb(file: str) -> Tuple[List[str], List[Molecule], List[str]]:
    header = []
    atoms = []
    footer = []
    with open(file) as f:
        l = f.readline()
        m = re.match(hetatm, l)
        while m is None:
            header.append(l)
            l = f.readline()
            m = re.match(hetatm, l)
        while m is not None:
            atoms.append(Atom(
                element = m.group(8).strip(),
                type = m.group(2).strip(),
                num = int(m.group(1)),
                res = int(m.group(4)),
                coords=[float(m.group(5)), float(m.group(6)), float(m.group(7))]
            ))
            if atoms[-1].type.isalpha():
                atoms[-1].type = '%s%d' % (atoms[-1].type, atoms[-1].num)
            l = f.readline()
            m = re.match(hetatm, l)
        m = re.match(conect, l)
        while m is None:
            l = f.readline()
            m = re.match(conect, l)
        bonds = [[i + 1] for i in range(len(atoms))]
        while m is not None:
            nums = [int(n) for n in m.group(0).split() if n.isdigit()]
            bonds[nums[0] - 1].extend(nums[1:])
            l = f.readline()
            m = re.match(conect, l)
        while l:
            footer.append(l)
            l = f.readline()
    mols = []
    for r in set([a.res for a in atoms]):
        mols.append(Molecule(
            atoms = [a for a in atoms if a.res == r],
            bonds = [b for i, b in enumerate(bonds) if atoms[i].res == r]
        ))
    return header, mols, footer
