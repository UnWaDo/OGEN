from typing import Tuple
import numpy as np


class Atom:
    def __init__(
        self,
        element: str,
        type: str,
        num: int,
        res: int = 1,
        coords: Tuple[float, float, float] = (0, 0, 0)
):
        self.element = element.upper()
        self.type = type
        self.num = num
        self.res = res
        self.coords = np.array(coords, dtype=float)

    def __eq__(self, other):
        if self.type != other.type:
            return False
        if self.num != other.num:
            return False
        return True

    def __str__(self):
        type_el = ''.join([l for l in self.type if l.isalpha()])
        type_num = self.type[len(type_el):]
        return 'HETATM%5d %2s%-2s UNL  %4d    %8.3f%8.3f%8.3f  1.00  0.00          %2s  ' % (
            self.num,
            type_el,
            type_num,
            self.res,
            self.coords[0],
            self.coords[1],
            self.coords[2],
            self.element
        )
