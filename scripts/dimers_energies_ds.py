import os
import pandas as pd
import sys


if len(sys.argv) < 3:
    print('Not enough arguments given, require reference file and energies file', file=sys.stderr)
    exit(1)

ref = pd.read_csv(sys.argv[1], index_col=0)
method, _ = os.path.splitext(os.path.basename(sys.argv[2]))
energies = pd.read_csv(sys.argv[2], index_col=0, header=None, names=[method])

formation_energies = energies.loc[[i.endswith('_diss') for i in energies.index]]
formation_energies.index = [i[:-5] for i in formation_energies.index]
if method in ref.columns:
    ref.update(formation_energies)
else:
    ref = pd.concat([ref, formation_energies], axis=1)
ref.to_csv(sys.argv[1])
