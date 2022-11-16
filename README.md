# OGEN: force field with off-site charges generator

## Installation

For now the only option available is to install module as source code. Module requires some additional programs to be able to run

### MultiWFN

### BOSS

## Usage

Simple way to use it is to run the module itself. It accepts the following arguments:

* path to `.fchk` file, which can be produced by Gaussian calculations;
  * real-space function to be used (ELF, LOL or Laplassian);
  * amount of OSCs on oxygen and sulfur atoms (1 or 2; mode 3 uses 2 points for C=O group and 1 elsewhere);
  * whether to add OSCs to fluorine in C−F line;
* alternatively, `--no-oscs` option can be provided, which will result in generating OPLS field with charges fitted by RESP.

## How does this work?