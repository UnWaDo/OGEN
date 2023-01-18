import pandas as pd
import sys


EXTRA_COLUMNS_NUM = 3

if len(sys.argv) < 2:
    print('Not enough arguments, require csv with energies', file=sys.stderr)
    exit(1)

energies_csv = sys.argv[1]
if len(sys.argv) > 2:
    output_csv = sys.argv[2]
    if not output_csv.endswith('.csv'):
        output_csv += '.csv'
else:
    output_csv = 'errors.csv'
energies = pd.read_csv(energies_csv, index_col=0)
ref = energies[energies.columns[EXTRA_COLUMNS_NUM]]

out_vdw = energies[energies.in_vdW == 0]
out_vdw_ref = out_vdw[out_vdw.columns[EXTRA_COLUMNS_NUM]]

columns = [
    'Method',
    'MAE, kcal/mol', 'Max. error, kcal/mol', 'Max. error, system',
    'MAE (vdW out), kcal/mol', 'Max. error (vdW out), kcal/mol', 'Max. error (vdW out), system',
]
errors = {c: [] for c in columns}
for c in energies.columns[EXTRA_COLUMNS_NUM + 1:]:
    tokens = c.split('_')
    diff = (energies[c] - ref).abs()
    out_vdw_diff = (out_vdw[c] - out_vdw_ref).abs()

    errors['Method'].append(c)
    errors['MAE, kcal/mol'].append(diff.mean())
    errors['Max. error, kcal/mol'].append(diff.max())
    errors['Max. error, system'].append(diff.idxmax()),
    errors['MAE (vdW out), kcal/mol'].append(out_vdw_diff.mean())
    errors['Max. error (vdW out), kcal/mol'].append(out_vdw_diff.max())
    errors['Max. error (vdW out), system'].append(out_vdw_diff.idxmax()),

errors = pd.DataFrame(errors)
errors.to_csv(output_csv)