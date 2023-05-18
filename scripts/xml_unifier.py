import os
import sys
from typing import List
from OGEN.ForceField import *


CODES = {
    'H': 1,
    'C': 3,
    'N': 5,
    'O': 7,
    'F': 11,
    'S': 13,
    'CL': 17,
    'BR': 23,
    'I': 29
}


def compound_code(elements: List[str]):
    code = 0
    found = {}
    for i, e in enumerate(elements):
        if found.get(e) is None:
            found[e] = 2
            code += CODES[e.upper()] * (i + 1)
        else:
            code += CODES[e.upper()] * found[e] * (i + 1)
            found[e] += 1
    code += 31 * len(elements)
    out = ''
    for i in range(3):
        out += chr(ord('A') + (code % 26))
        code //= 26
    return out


def rename_classes(classes, old_classes, new_classes):
    return tuple(new_classes[old_classes.index(c)] for c in classes)


if len(sys.argv) != 3:
    print('Invalid number of arguments\nUsage: python xml_unifier.py file.xml/dir [output_dir]')

xml_file = sys.argv[1]
if os.path.isdir(xml_file):
    xml_files = [os.path.join(xml_file, f)
        for f in os.listdir(xml_file) if f.endswith('.xml')]
else:
    xml_files = [xml_file]
OUTPUT_DIR = sys.argv[2]

for xml_file in xml_files:
    xml_name, ext = os.path.splitext(os.path.basename(xml_file))
    ff = ForceField.from_file(xml_file)
    name = compound_code([a.element for a in ff.atom_types if a.element is not None])
    old_types = []
    old_classes = []
    new_classes = []
    for i, at in enumerate(ff.atom_types):
        old_types.append(at.name)
        old_classes.append(at.class_name)
        new_classes.append('cl%s_%s_%d' % (
            at.element if at.element is not None else 'X',
            name,
            i
        ))
        at.name = new_classes[-1][2:]
        at.class_name = new_classes[-1]
    ff.residues[0].name = name
    for b in ff.bond_forces:
        b.classes = rename_classes(b.classes, old_classes, new_classes)
    for a in ff.angle_forces:
        a.classes = rename_classes(a.classes, old_classes, new_classes)
    for t in ff.torsion_forces:
        t.classes = rename_classes(t.classes, old_classes, new_classes)
    if not os.path.exists(OUTPUT_DIR):
        os.makedirs(OUTPUT_DIR)
    ff.to_xml(os.path.join(OUTPUT_DIR, '%s.xml' % xml_name))
