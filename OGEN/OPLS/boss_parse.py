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


def get_bonds(file) -> List[Dict[str, Union[Tuple[int, int], float]]]:
    for i in range(2):
        line = file.readline()
    bonds_energy = []
    data = line.split()
    while len(data):
        bonds_energy.append({
            'atoms': (int(data[0]) - 4, int(data[1]) - 4),
            'r_eq': float(data[2]) / 10, # OPLS in Angstrom, OpenMM in nm
            'k': float(data[3]) * 2 * 4.184 * 100 # OPLS in kcal/mol/Å^2, OpenMM in kJ/mol/nm^2 and with 1/2k
        })
        line = file.readline()
        data = line.split()
    return bonds_energy


def get_angles(file) -> List[Dict[str, Union[Tuple[int, int, int], float]]]:
    for i in range(2):
        line = file.readline()
    angles_energy = []
    data = line.split()
    while len(data):
        angles_energy.append({
            'atoms': tuple(int(d) - 4 for d in data[:3]),
            'a_eq': float(data[3]) / 180 * math.pi, # OPLS in degrees, OpenMM in radians
            'k': float(data[4]) * 2 * 4.184 # OPLS in kcal/mol/rad^2, OpenMM in kJ/mol/rad^2 and with 1/2k
        })
        line = file.readline()
        data = line.split()
    return angles_energy


def get_torsions(file, torsions: List[Tuple[int, int, int, int]]) -> List[
    Dict[str, Union[Tuple[int, int, int, int], float, bool]]
]:
    for i in range(2):
        line = file.readline()
    tors_energy = []
    data = line.split()
    i = 0
    while len(data):
        tors_energy.append({
            'atoms': tuple(t for t in torsions[i][:-1]),
            'vs': [float(d) * 4.184 / 2 for d in data[4:8]], # OPLS in kcal/mol, OpenMM in kJ/mol and without 1/2V
            'phases': [float(d) for d in data[8:11]],
            'improper': torsions[i][-1]
        })
        i += 1
        line = file.readline()
        data = line.split()
    if i < len(torsions):
        raise Exception('Not enough torsions in cout file')
    return tors_energy

def get_bonded(torsions: List[Tuple[int, int, int, int]]) -> Tuple[
    List[Dict[str, Union[Tuple[int, int], float]]],
    List[Dict[str, Union[Tuple[int, int, int], float]]],
    List[Dict[str, Union[Tuple[int, int, int, int], float, bool]]]
]:
    if not os.path.exists('cout'):
        raise Exception('BOSS out file is absent. Probably, BOSS is not installed, BOSSdir is not set or the calculation failed')
    bonds_energy = []
    angles_energy = []
    tors_energy = []
    with open('cout', 'r') as params:
        line = params.readline()
        while line:
            if 'Dihedral' in line:
                tors_energy = get_torsions(params, torsions)
            elif 'Bond Stretching Parameters' in line:
                bonds_energy = get_bonds(params)
            elif 'Angle Bending Parameters' in line:
                angles_energy = get_angles(params)
            line = params.readline()
    return bonds_energy, angles_energy, tors_energy
