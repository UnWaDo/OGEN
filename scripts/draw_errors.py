import pandas as pd
import matplotlib
font = {'family' : 'Times New Roman',
        'size'   : 12}
matplotlib.rc('font', **font)
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Patch
import sys, os



LOCALE = 'ru'


if len(sys.argv) < 2:
    print('Not enough arguments, require csv with errors', file=sys.stderr)
    exit(1)

errors_csv = sys.argv[1]
if len(sys.argv) > 2:
    output = sys.argv[2]
    if not output.endswith('.png'):
        output += '.png'
else:
    output, _ = os.path.splitext(os.path.basename(sys.argv[1]))
    output = '%s.png' % output

errors = pd.read_csv(errors_csv, index_col=0)
errors = errors.sort_values(by=['MAE, kcal/mol'], ignore_index=True)
colors = [
    '#1f78b4' if 'ELF' in r else
        '#33a02c' if 'LOL' in r else
            '#e31a1c' if 'Lap' in r else 
                '#ff7f00' if 'qube' in r else '#fdbf6f'
    for r in errors['Method']
]
hatch = [
    '//' if 'qube' in r else
        '' if 'opls' in r else
            '\\' if 'gfnff' in r else '-'
    for r in errors['Method']
]
edgecolor = [
    'k' if 'qube' in r else
        'white' if 'opls' in r else
            'gray' if 'gfnff' in r else 'white'
    for r in errors['Method']
]

ru = {
    'rsf': 'Функция',
    'qube_oscs': 'ВАЗ из QUBEKit',
    'no_oscs': 'Без ВАЗ',
    'qube_charges': 'Заряды QUBEKit',
    'opls_charges': 'Заряды OPLS',
    'resp_charges': 'Заряды RESP',
    'gfn_charges': 'Заряды GFN-FF',
    'qube_disp': 'Дисперсия QUBEKit',
    'opls_disp': 'Дисперсия OPLS',
    'gfn_disp': 'Дисперсия GFN-FF',
    'y': 'Средняя абсолютная ошибка, ккал/моль'
}
en = {
    'rsf': 'RSF',
    'qube_oscs': 'OSCs from QUBEKit',
    'no_oscs': 'No OSCs',
    'qube_charges': 'QUBEKit charges',
    'opls_charges': 'OPLS charges',
    'resp_charges': 'RESP charges',
    'gfn_charges': 'GFN-FF charges',
    'qube_disp': 'QUBEKit dispersion',
    'opls_disp': 'OPLS dispersion',
    'gfn_disp': 'GFN-FF dispersion',
    'y': 'Mean absolute error, kcal/mol'
}
text = {'ru': ru, 'en': en}
size = (6.5, 5) if 'X40' not in errors_csv else (6.5, 5.5)

legend_elements = [
    Patch(facecolor='#1f78b4', label='%s: ELF' % text[LOCALE]['rsf']),
    Patch(facecolor='#33a02c', label='%s: LOL' % text[LOCALE]['rsf']),
    Patch(facecolor='#e31a1c', label='%s: Lap' % text[LOCALE]['rsf']),
    Patch(facecolor='#ff7f00', label=text[LOCALE]['qube_oscs']),
    Patch(facecolor='#fdbf6f', label=text[LOCALE]['no_oscs']),
    Patch(facecolor='white', edgecolor='k', hatch='///', label=text[LOCALE]['qube_charges']),
    Patch(facecolor='white', edgecolor='k', hatch='', label=text[LOCALE]['opls_charges']),
    Patch(facecolor='white', edgecolor='k', hatch='----', label=text[LOCALE]['resp_charges']),
    Patch(facecolor='white', edgecolor='k', hatch='\\\\\\', label=text[LOCALE]['gfn_charges']),
    Line2D([0], [0], color='k', label=text[LOCALE]['qube_disp']),
    Line2D([0], [0], color='none', label=text[LOCALE]['opls_disp']),
    Line2D([0], [0], color='gray', label=text[LOCALE]['gfn_disp']),
]
fig, ax = plt.subplots()
fig.set_size_inches(*size)
bars = ax.bar(errors.index, errors['MAE, kcal/mol'], color=colors, edgecolor=edgecolor)
for i, b in enumerate(bars):
    b.set_hatch(hatch[i])
labels = [
    ('1' if '1' in r else '2' if '2' in r else '3' if '3' in r else '') + ('F' if r.endswith('F') else '')
    for r in errors['Method']
]
for i, b in enumerate(bars.patches):
    ax.text(b.get_x()+b.get_width() * 0.4, b.get_y() + b.get_height()+0.02, labels[i])
ax.set_ylabel(text[LOCALE]['y'])
ax.get_xaxis().set_visible(False)
ax.set_axisbelow(True)
# ax.legend(handles=legend_elements, ncol=2)
ax.legend(handles=legend_elements, ncol=3, loc='upper center', bbox_to_anchor=(-0.035, -0.35, 1, 0.32))
ax.grid(visible=True)
fig.savefig(output, dpi=300, bbox_inches='tight')

fig, ax = plt.subplots()
fig.set_size_inches(*size)
errors = errors.sort_values(by=['MAE (vdW out), kcal/mol'], ignore_index=True)
colors = [
    '#1f78b4' if 'ELF' in r else
        '#33a02c' if 'LOL' in r else
            '#e31a1c' if 'Lap' in r else 
                '#ff7f00' if 'qube' in r else '#fdbf6f'
    for r in errors['Method']
]
hatch = [
    '//' if 'qube' in r else
        '' if 'opls' in r else
            '\\' if 'gfnff' in r else '-'
    for r in errors['Method']
]
edgecolor = [
    'k' if 'qube' in r else
        'white' if 'opls' in r else
            'gray' if 'gfnff' in r else 'white'
    for r in errors['Method']
]
bars = ax.bar(errors.index, errors['MAE (vdW out), kcal/mol'], color=colors, edgecolor=edgecolor)
for i, b in enumerate(bars):
    b.set_hatch(hatch[i])
labels = [
    ('1' if '1' in r else '2' if '2' in r else '3' if '3' in r else '') + ('F' if r.endswith('F') else '')
    for r in errors['Method']
]
for i, b in enumerate(bars.patches):
    ax.text(b.get_x()+b.get_width() * 0.4 - 0.2 * (len(labels[i]) - 1), b.get_y() + b.get_height()+0.02, labels[i])
ax.set_ylabel(text[LOCALE]['y'])
ax.get_xaxis().set_visible(False)
ax.set_axisbelow(True)
# ax.legend(handles=legend_elements, ncol=2)
ax.legend(handles=legend_elements, ncol=3, loc='upper center', bbox_to_anchor=(-0.035, -0.35, 1, 0.32))
ax.grid(visible=True)
fig.savefig(output[:-4] + '_vdw.png', dpi=300, bbox_inches='tight')
