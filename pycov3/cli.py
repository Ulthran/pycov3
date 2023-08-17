import argparse
import logging
import sys

from .Directory import FastaDir, SamDir

###
# Main steps:
# 1. Check and record parameters designated by user
# 2. Record GC content in sliding windows for contigs
# 	key subroutine:
# 	&GC_count
# 3. Distribute SAM files to threads and calculate coverages in sliding windows
# 	key subroutine:
# 	&sam_cov_parallel
# 4. Convert COV2 (coverage in samples) to COV3 (coverage in contig clusters)
# 5. Distribute COV3 files to threads and invoke R to estimate growth rates
# 	key subroutine:
# 	&cov3_estPTR_parallel
# 6. Output
###


def main(argv=None):
    p = argparse.ArgumentParser()
    p.add_argument(
        "-S", "--sam_dir", help="the directory containing the sam file(s)"
    )
    p.add_argument(
        "-F",
        "--fasta_dir",
        help="the directory containing the binned fasta file(s)",
    )
    p.add_argument(
        "-O",
        "--out_dir",
        help="the output directory",
    )
    p.add_argument(
        "-X",
        "--overwrite",
        help="overwrites any existing outputs",
        action="store_true",
    )

    p.add_argument(
        "-W",
        "--window_size",
        help="size (nt) of window for calculation of coverage (default: 5000)",
    )
    p.add_argument(
        "-D",
        "--window_step",
        help="step (nt) of window for calculation of coverage (default: 100)",
    )
    p.add_argument(
        "-M",
        "--mapq_cutoff",
        help="cutoff of mapping quality when calculating coverages (default: 5)",
    )
    p.add_argument(
        "-L",
        "--mapl_cutoff",
        help="cutoff of mapping length when calculating coverages (default: 50)",
    )
    p.add_argument(
        "-R",
        "--max_mismatch_ratio",
        help="maximum of mismatch ratio for each read as a hit (default: 0.03)",
        default=0.03,
    )
    p.add_argument(
        "-T",
        "--thread_num",
        help="set number of threads for parallel running (default: 1)",
        default=1,
    )

    p.add_argument(
        "-G",
        "--log_level",
        type=int,
        help="Sets the log level, default is info, 10 for debug (Default: 20)",
        default=20,
    )

    args = p.parse_args(argv)
    if not (args.S and args.F and args.O):
        p.print_help(sys.stderr)
        sys.exit(1)
    logging.basicConfig()
    logging.getLogger().setLevel(args.log_level)

    sam_d = SamDir(args.S, args.X)

    coverage_params = {
        "window_size": args.W,
        "window_step": args.D,
        "edge_length": sam_d.calculate_edge_length,
        "mapq_cutoff": args.M,
        "mapl_cutoff": args.L,
    }
    coverage_params = {
        k: v for k, v in coverage_params.items() if v is not None
    }

    fasta_d = FastaDir(args.F, args.X, coverage_params)
