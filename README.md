# pycov3

[![CI](https://github.com/Ulthran/pycov3/actions/workflows/main.yml/badge.svg)](https://github.com/Ulthran/pycov3/actions/workflows/main.yml)
[![codecov](https://codecov.io/gh/Ulthran/pycov3/branch/main/graph/badge.svg?token=pycov3_token_here)](https://codecov.io/gh/Ulthran/pycov3)
[![Super-Linter](https://github.com/Ulthran/pycov3/actions/workflows/linter.yml/badge.svg)](https://github.com/Ulthran/pycov3/actions/workflows/linter.yml)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/a2c7aa4e154d4bae82246d5f950afa9c)](https://app.codacy.com/gh/Ulthran/pycov3/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade)

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

## The Makefile

All the utilities for the template and project are on the Makefile

```bash
❯ make
Usage: make <target>

Targets:
help:             ## Show the help.
install:          ## Install the project in dev mode.
fmt:              ## Format code using black & isort.
lint:             ## Run pep8, black, mypy linters.
test: lint        ## Run tests and generate coverage report.
watch:            ## Run tests on every change.
clean:            ## Clean unused files.
virtualenv:       ## Create a virtual environment.
release:          ## Create a new tag for release.
docs:             ## Build the documentation.
switch-to-poetry: ## Switch to poetry package manager.
init:             ## Initialize the project based on an application template.
```

## Development

Read the [CONTRIBUTING.md](CONTRIBUTING.md) file.
