from pathlib import Path
from pycov3.File import FastaFile
from tests.unit.utils import create_sample_fasta_file


def test_fasta_file_init():
    # Test FastaFile initialization
    fasta_file_path = Path("max_bin.001.fasta")
    create_sample_fasta_file(fasta_file_path)

    fasta_file = FastaFile(fasta_file_path)
    assert fasta_file.sample == "max_bin"
    assert fasta_file.bin == "001"


def test_fasta_file_parse():
    # Test FastaFile parse method
    fasta_file_path = Path("max_bin.001.fasta")
    create_sample_fasta_file(fasta_file_path)

    fasta_file = FastaFile(fasta_file_path)
    contigs = list(fasta_file.parse())

    assert len(contigs) == 2
    assert contigs[0][0] == "Contig1"
    assert contigs[0][1] == "ATCGATCG"
    assert contigs[1][0] == "Contig2"
    assert contigs[1][1] == "GCTAGCTA"


def test_fasta_file_write():
    # Test FastaFile write method (not implemented)
    pass
