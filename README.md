# OGEN: force field with off-site charges generator

## Installation

Install required packages from environment.yml file.
Alternative is manual installation of `rdkit` (>= 2022.09.01), `scikit-learn` and `openbabel`. All packages are available from `conda-forge` and can be installed by running `conda install -c conda-forge package-name` in new environment. All testing was performed with Python 3.9.13 and Python 3.10.12

```
# From environment file
conda env create -f environment.yml

## OR manually
conda create -n your-env-name -yc conda-forge python=3.9.13 'rdkit>=2022.09.01' scikit-learn openbabel
```

`mamba` and `micromamba` normally works much faster than usuall `conda` or `miniconda`.

For now the only option available is to install module as source code. Module requires some additional programs to be able to run

### MultiWFN

Install from http://sobereva.com/multiwfn/download.html.
Make it accessible as `Multiwfn` from command line (by adding to `$PATH`)
or change `MWFN_EXECUTABLE` variable in `OGEN/MultiWFN/CriticalPoints.py` to match full path.

### BOSS

Required to generate full forcefield.
Bond parameters generation can be disabled by passing `--only-charges` argument to command line

## Usage

Simple way to use it is to run the module itself. It accepts the following arguments:

* path to `.fchk` file, which can be produced by Gaussian calculations;
  * real-space function to be used (ELF, LOL or Laplassian);
  * amount of OSCs on oxygen and sulfur atoms (1 or 2; mode 3 uses 2 points for C=O group and 1 elsewhere);
  * whether to add OSCs to fluorine in C−F line;
* alternatively, `--no-oscs` option can be provided, which will result in generating OPLS field with charges fitted by RESP.

If you want to use scripts, add repository folder to `PYTHONPATH`

```
export PYTHONPATH="$PYTHONPATH:/path/to/cloned/repo"
```

## Examples

```bash
# Generate FF containing only charges. Default mode is ELF3
python run.py --only-charges path/to/file.fchk

# Do not add ring OSC
python run.py --no-ring path/to/file.fchk

# Use other RSF
python run.py --rsf LOL path/to/file.fchk

# Use other OSC count
python run.py --oscs-count 1 LOL path/to/file.fchk

# Get help message with desctiption of all available flags
python run.py --help
```

## How does this work?