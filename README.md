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

Create a SAM directory and FASTA directory, set any non-default window or coverage parameters, then create a COV3 directory and use it to generate a COV3 file for each contig set in the FASTA directory.

```py
    from pycov3.Directory import Cov3Dir, FastaDir, SamDir

    sam_d = SamDir(Path("/path/to/sams/"), False)

    window_params = {
        "window_size": None,
        "window_step": None,
        "edge_length": sam_d.calculate_edge_length(),
    }
    coverage_params = {
        "mapq_cutoff": None,
        "mapl_cutoff": None,
        "max_mismatch_ratio": None,
    }
    window_params = {k: v for k, v in window_params.items() if v is not None}
    coverage_params = {k: v for k, v in coverage_params.items() if v is not None}

    fasta_d = FastaDir(Path("/path/to/fastas/"), False)

    cov3_d = Cov3Dir(
        Path(args.out_dir),
        False,
        fasta_d.get_filenames(),
        window_params,
        coverage_params,
    )

    cov3_d.generate(sam_d, fasta_d)
```

```bash
$ python -m pycov3 -h
#or
$ pycov3 -h
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
