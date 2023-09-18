from pathlib import Path
from pycov3.File import SamFile


def create_sample_sam_file(sam_file_path):
    # Create a sample SAM file for testing
    sample_sam_content = (
        "@SQ\tSN:Contig1\tLN:100\n"
        "read1\t0\tContig1\t1\t30\t10M\t*\t0\t0\tATCGATCGAT\n"
        "read2\t0\tContig1\t11\t30\t10M\t*\t0\t0\tGCTAGCTAGC\n"
    )
    with open(sam_file_path, "w") as f:
        f.write(sample_sam_content)


def test_sam_file_init():
    # Test SamFile initialization
    sam_file_path = Path("Akk_001.sam")
    create_sample_sam_file(sam_file_path)

    sam_file = SamFile(sam_file_path)
    assert sam_file.sample == "Akk"
    assert sam_file.bin == "001"


def test_sam_file_parse():
    # Test SamFile parse method
    sam_file_path = Path("Akk_001.sam")
    create_sample_sam_file(sam_file_path)

    sam_file = SamFile(sam_file_path)
    parsed_reads = list(sam_file.parse())

    assert len(parsed_reads) == 2
    assert parsed_reads[0]["read_name"] == "read1"
    assert parsed_reads[0]["flag"] == 0
    assert parsed_reads[0]["reference_name"] == "Contig1"
    assert parsed_reads[0]["position"] == 1
    assert parsed_reads[0]["mapping_quality"] == 30
    assert parsed_reads[0]["cigar"] == "10M"
    assert parsed_reads[0]["mismatch"] == 0


def test_sam_file_parse_contig_lengths():
    # Test SamFile parse_contig_lengths method
    sam_file_path = Path("Akk_001.sam")
    create_sample_sam_file(sam_file_path)

    sam_file = SamFile(sam_file_path)
    contig_lengths = sam_file.parse_contig_lengths()

    assert len(contig_lengths) == 1
    assert contig_lengths["Contig1"] == 100


def test_sam_file_write():
    # Test SamFile write method (not implemented)
    pass
