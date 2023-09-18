from pathlib import Path
from pycov3.File import Cov3File, SamFile
from tests.unit.utils import (
    create_sample_cov3_file,
    create_sample_fasta_file,
    create_sample_sam_file,
)


def test_cov3_file_init():
    # Test cov3 file initialization
    cov3_file = Cov3File(Path("max_bin.001.cov3"), "001")
    assert cov3_file.bin_name == "001"


def test_cov3_file_parse():
    # Test cov3 file parse
    cov3_fp = Path("max_bin.001.cov3")
    create_sample_cov3_file(cov3_fp)
    cov3_file = Cov3File(cov3_fp, "001")

    cov3_vals = list(cov3_file.parse())

    assert len(cov3_vals) == 4
    assert cov3_vals[0]["log_cov"] == 1.234
    assert cov3_vals[0]["GC_content"] == 0.567
    assert cov3_vals[0]["sample"] == "sample1"
    assert cov3_vals[0]["contig"] == "contig1"
    assert cov3_vals[0]["length"] == 100


def test_cov3_file_parse_sample_contig():
    # Test cov3 file parse by sample and contig
    cov3_fp = Path("max_bin.001.cov3")
    create_sample_cov3_file(cov3_fp)
    cov3_file = Cov3File(cov3_fp, "001")

    cov3_vals = list(cov3_file.parse_sample_contig())

    assert len(cov3_vals) == 3
    assert cov3_vals[0]["sample"] == "sample1"
    assert cov3_vals[0]["contig"] == "contig1"
    assert len(cov3_vals[0]["log_covs"]) == 2


def test_cov3_file_update_coverages():
    # Test cov3 file update_coverages utility
    sam_file_path = Path("Akk_001.sam")
    create_sample_sam_file(sam_file_path)
    sam_file = SamFile(sam_file_path)

    cov3_fp = Path("max_bin.001.cov3")
    create_sample_cov3_file(cov3_fp)
    cov3_file = Cov3File(cov3_fp, "001")

    sam_lines = list(sam_file.parse())
    coverages = {}

    for line in sam_lines:
        coverages = cov3_file._Cov3File__update_coverages(coverages, line, 2, 2)
        print(coverages)

    assert coverages == {1: 1}


def test_cov3_file_calculate_mapl():
    # Test cov3 file mapping length calculation
    assert Cov3File.calculate_mapl("250M") == 250
    assert Cov3File.calculate_mapl("100M50I100M") == 150
    assert Cov3File.calculate_mapl("200D50I") == 150
    assert Cov3File.calculate_mapl("250H") == 0
    assert Cov3File.calculate_mapl("*") == -1
