# pycov3

[![Tests](https://github.com/Ulthran/pycov3/actions/workflows/pr.yml/badge.svg)](https://github.com/Ulthran/pycov3/actions/workflows/pr.yml)
[![codecov](https://codecov.io/gh/Ulthran/pycov3/graph/badge.svg?token=HH27P1FDM5)](https://codecov.io/gh/Ulthran/pycov3)
[![Renovate enabled](https://img.shields.io/badge/renovate-enabled-brightgreen.svg)](https://renovatebot.com/)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/a2c7aa4e154d4bae82246d5f950afa9c)](https://app.codacy.com/gh/Ulthran/pycov3/dashboard?utm_source=gh&utm_medium=referral&utm_content=&utm_campaign=Badge_grade)
[![PyPI](https://badge.fury.io/py/pycov3.svg)](https://pypi.org/project/pycov3/)
[![Bioconda](https://anaconda.org/bioconda/pycov3/badges/downloads.svg)](https://anaconda.org/bioconda/pycov3/)
[![DockerHub](https://img.shields.io/docker/pulls/ctbushman/pycov3)](https://hub.docker.com/repository/docker/ctbushman/pycov3/)

A package for generating cov3 files which are generated from sam files giving coverage information and a fasta file giving binned contigs. Cov3 files are used as input for the [DEMIC R package](https://github.com/Ulthran/DEMIC) which calculates PTR, an estimate for bacterial growth rates.

## Installation

### PyPi

```bash
pip install pycov3
pycov3 -h
```

### Bioconda

```bash
conda create -n pycov3 -c conda-forge -c bioconda pycov3
conda activate pycov3
pycov3 -h
```

### DockerHub

```bash
docker pull ctbushman/pycov3:latest
docker run --rm --name pycov3 pycov3 pycov3 -h
```

### GitHub

```bash
git clone https://github.com/Ulthran/pycov3.git
cd pycov3/
pip install .
pycov3 -h
```

## Usage

Use `-h` to see options for running the CLI.

```bash
$ pycov3 -h
```

The FASTAs should all be in one directory with names of the format `{sample}.{bin_name}.fasta/.fa/.fna` and the SAMs should also all be in one directory with names of the format `{sample}_{bin_name}.sam`. The output COV3 files will be written to a directory with names of the format `{sample}.{bin_name}.cov3`.

You can also use the library in your own code. Create a SAM directory and FASTA directory, set any non-default window or coverage parameters, then create a COV3 directory and use it to generate a COV3 file for each contig set in the FASTA directory.

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

Alternatively, to use the bare application logic and do all the file handling yourself, you can use the `Cov3Generator` class which takes a list of generators as SAM inputs and a generator as a FASTA input.

```py
    from pycov3.Cov3Generator import Cov3Generator
    from pycov3.File import Cov3File

    cov3_generator = Cov3Generator(
        sam_generators,
        fasta_generator,
        sample,
        bin_name,
        window_params,
        **coverage_params,
    )

    cov3_dict = cov3_generator.generate_cov3()

    # Write output
    cov3_file = Cov3File(Path(/path/to/output/), "001")
    cov3_file.write_generator(cov3_generator.generate_cov3())
```

## Resource Requirements

Threads: pycov3 uses `multiprocessing` to parallelize processing of input fastas. Increasing `--thread_num` up to the number of input fastas should improve runtime, with no benefits beyond that number.

Memory: pycov3 uses generators as much as possible. The main memory users are the `Contig` objects, which each hold a contig's sequence and information for each `Window` over its length. There is also a `coverages` dictionary that could potentially grow to the size of the largest contig (but that is very unlikely). At a minimum, twice the size of the largest contig should be given per thread.

Algorithmic Complexity: Assuming enough threads are provided to have each fasta file processed separately, the time complexity is roughly `O(cwsr)`.

`c`: Number of contigs in fasta
`s`: Number of sam files
`w`: Max number of windows per contig
`r`: Max number of records per sam file

## Help

Please use the [Issues](https://github.com/Ulthran/pycov3/issues) on this repo for any problems, questions, or suggestions.