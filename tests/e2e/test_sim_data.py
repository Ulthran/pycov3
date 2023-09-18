import logging
from pathlib import Path
from pycov3.Directory import SamDir, FastaDir, Cov3Dir

def test_sim_data():
    logging.basicConfig()
    logging.getLogger().setLevel(10)

    sim_fp = Path("sim")
    sams_fp = sim_fp / "sams"
    fastas_fp = sim_fp / "fastas"
    output_fp = sim_fp / "output"
    expected_output_fp = sim_fp / "expected_output"

    overwrite = True

    sam_d = SamDir(sams_fp, overwrite)

    window_params = {
        "window_size": 500,
        "window_step": 10,
        "edge_length": sam_d.calculate_edge_length(),
    }
    coverage_params = {
        "mapq_cutoff": None,
        "mapl_cutoff": None,
        "max_mismatch_ratio": None,
    }
    window_params = {
        k: v for k, v in window_params.items() if v is not None
    }
    coverage_params = {
        k: v for k, v in coverage_params.items() if v is not None
    }

    fasta_d = FastaDir(fastas_fp, overwrite, window_params)

    cov3_d = Cov3Dir(output_fp, overwrite, fasta_d.get_filenames(), coverage_params)

    cov3_d.generate(sam_d, fasta_d)

    import shutil
    #shutil.copyfile(output_fp / "max_bin.002.cov3", "/mnt/d/Penn/pycov3/tests/data/sim/expected_output/max_bin.002.cov3")

    cov3_1 = cov3_d.get_bin("001")
    cov3_2 = cov3_d.get_bin("002")

    for sample_contig in cov3_1.parse_sample_contig():
        print(sample_contig)

    with open(expected_output_fp / "max_bin.001.cov3") as exp_f, open(cov3_1.fp) as out_f:
        assert exp_f.readlines() == out_f.readlines()
    with open(expected_output_fp / "max_bin.002.cov3") as exp_f, open(cov3_2.fp) as out_f:
        assert exp_f.readlines() == out_f.readlines()