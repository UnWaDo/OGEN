
import math
import os
from typing import Dict, List, Tuple, Union


def get_dispersion() -> List[Dict[str, Union[str, float]]]:
    if not os.path.exists('sum'):
        raise Exception('BOSS sum file is absent. Probably, BOSS is not installed, BOSSdir is not set or the calculation failed')
    with open('sum', 'r') as params:
        line = params.readline()
        while line and 'Final Non-Bonded Parameters' not in line:
            line = params.readline()
        for i in range(2):
            line = params.readline()
        disp = []
        data = line.split()
        while len(data):
            disp.append({
                'class': data[2],
                'charge': float(data[3]),
                'sigma': float(data[4]) / 10, # OPLS in Angstrom, OpenMM in nm
                'epsilon': float(data[5]) * 4.184, # OPLS in kcal/mol, OpenMM in kJ/mol
            })
            line = params.readline()
            data = line.split()
    return disp


def get_bonded() -> Tuple[
    List[Dict[str, Union[Tuple[int, int], float]]],
    List[Dict[str, Union[Tuple[int, int, int], float]]],
    List[Dict[str, Union[Tuple[int, int, int, int], float]]]
]:
    if not os.path.exists('cout'):
        raise Exception('BOSS out file is absent. Probably, BOSS is not installed, BOSSdir is not set or the calculation failed')
    with open('cout', 'r') as params:
        line = params.readline()
        while line and 'Bond Stretching Parameters' not in line:
            line = params.readline()
        for i in range(2):
            line = params.readline()
        bonds_energy = []
        data = line.split()
        while len(data):
            bonds_energy.append({
                'atoms': (int(data[0]) - 4, int(data[1]) - 4),
                'r_eq': float(data[2]) / 10, # OPLS in Angstrom, OpenMM in nm
                'k': float(data[3]) * 2 * 4.184 * 100 # OPLS in kcal/mol/Å^2, OpenMM in kJ/mol/nm^2 and with 1/2k
            })
            line = params.readline()
            data = line.split()
        while line and 'Angle Bending Parameters' not in line:
            line = params.readline()
        for i in range(2):
            line = params.readline()
        angles_energy = []
        data = line.split()
        while len(data):
            angles_energy.append({
                'atoms': tuple(int(d) - 4 for d in data[:3]),
                'a_eq': float(data[3]) / 180 * math.pi, # OPLS in degrees, OpenMM in radians
                'k': float(data[4]) * 2 * 4.184 # OPLS in kcal/mol/rad^2, OpenMM in kJ/mol/rad^2 and with 1/2k
            })
            line = params.readline()
            data = line.split()
    return bonds_energy, angles_energy, []
