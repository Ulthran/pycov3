def create_sample_sam_file(sam_file_path):
    # Create a sample SAM file for testing
    sample_sam_content = (
        "@SQ\tSN:Contig1\tLN:100\n"
        "read1\t0\tContig1\t1\t30\t10M\t*\t0\t0\tATCGATCGAT\n"
        "read2\t0\tContig1\t11\t30\t10M\t*\t0\t0\tGCTAGCTAGC\n"
    )
    with open(sam_file_path, "w") as f:
        f.write(sample_sam_content)


def create_sample_fasta_file(fasta_file_path):
    # Create a sample FastaFile for testing
    sample_fasta_content = ">Contig1\nATCGATCG\n>Contig2\nGCTAGCTA\n"
    with open(fasta_file_path, "w") as f:
        f.write(sample_fasta_content)
