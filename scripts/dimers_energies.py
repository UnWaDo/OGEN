import os
import pandas as pd
import sys


def calc_formation_energy(dimer):
    return (energies.loc[dimer.name] \
        - energies.loc[dimer['monomer1']] \
            - energies.loc[dimer['monomer2']])['energy']


if len(sys.argv) < 3:
    print('Not enough arguments given, require reference file and energies file', file=sys.stderr)
    exit(1)

ref = pd.read_csv(sys.argv[1], index_col=0)
method, _ = os.path.splitext(os.path.basename(sys.argv[2]))
energies = pd.read_csv(sys.argv[2], index_col=0, header=None, names=['energy'])

formation_energies = ref.apply(calc_formation_energy, axis=1).rename(method)
if method in ref.columns:
    ref.update(formation_energies)
else:
    ref = pd.concat([ref, formation_energies], axis=1)
ref.to_csv(sys.argv[1])
