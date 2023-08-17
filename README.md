# pycov3

[![codecov](https://codecov.io/gh/Ulthran/pycov3/branch/main/graph/badge.svg?token=pycov3_token_here)](https://codecov.io/gh/Ulthran/pycov3)
[![CI](https://github.com/Ulthran/pycov3/actions/workflows/main.yml/badge.svg)](https://github.com/Ulthran/pycov3/actions/workflows/main.yml)

A package for generating cov3 files which are generated from sam files giving coverage information and a fasta file giving binned contigs. Cov3 files are used as input for the [DEMIC R package](https://github.com/Ulthran/DEMIC) which calculates PTR, an estimate for bacterial growth rates.

## Install it from PyPI

```bash
pip install pycov3
```

## Usage

```py
from pycov3 import BaseClass
from pycov3 import base_function

BaseClass().base_method()
base_function()
```

```bash
$ python -m pycov3
#or
$ pycov3
```

## Development

Read the [CONTRIBUTING.md](CONTRIBUTING.md) file.
