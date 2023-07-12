import re
import sys
from typing import List, Tuple, Union, TextIO
from . import Atom, Molecule


def to_pdb(
        mols: Union[Molecule, List[Molecule]],
        title: str = None,
        author: str = None,
        file: Union[str, TextIO] = sys.stdout):
    if type(mols) != list:
        mols = [mols]

    should_close = False
    if type(file) == str:
        file = open(file, 'w+')
        should_close = True

    if title is None:
        title = 'mol'
    if author is None:
        author = 'FilesProcessing'
    
    file.write('COMPND    %s\n' % title)
    file.write('AUTHOR    %s\n' % author)
    for mol in mols:
        file.write(mol.create_hetatm_section())
    for mol in mols:
        file.write(mol.create_conect_section())
    
    num = sum([len(m.atoms) for m in mols])
    file.write('MASTER     %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d %4d\n' % (
        0, 0, 0, 0, 0, 0, 0, 0, num, 0, num, 0
    ))
    file.write('END\n')
    if should_close:
        file.close()
