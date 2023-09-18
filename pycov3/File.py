import logging
import math
from abc import ABC, abstractmethod
from collections import OrderedDict
from itertools import groupby
from pathlib import Path

from .Sequence import Contig


class File(ABC):
    def __init__(self, fp: Path) -> None:
        super().__init__()
        self.fp = fp

    def exists(self) -> bool:
        if not self.fp.exists():
            logging.error(f"{self.fp} does not exist")
            return False
        if self.fp.stat().st_size == 0:
            logging.error(f"{self.fp} is empty")
            return False
        return True

    @abstractmethod
    def parse(self) -> None:
        pass

    @abstractmethod
    def write(self) -> None:
        pass


class FastaFile(File):
    def __init__(self, fp: Path) -> None:
        super().__init__(fp)

        try:
            stem = self.fp.stem.split(".")
            self.sample = stem[0]
            self.bin_name = stem[1]
        except IndexError:
            raise ValueError(
                f"FASTA filename {self.fp} not of format {{sample}}.{{bin_name}}.fasta/.fa/.fna"
            )

    def parse(self) -> list:
        with open(self.fp) as f:
            faiter = (x[1] for x in groupby(f, lambda line: line[0] == ">"))

            for header in faiter:
                # drop the ">"
                header_str = header.__next__()[1:].strip().split(" ")[0]

                # join all sequence lines to one.
                seq_str = "".join(s.strip() for s in faiter.__next__())

                yield (header_str, seq_str)

    def write(self):
        pass


class SamFile(File):
    def __init__(self, fp: Path) -> None:
        super().__init__(fp)
        stem = self.fp.stem.split("_")
        self.sample = stem[0]
        self.bin_name = stem[1].split(".")[0]

    def parse(self) -> list:
        with open(self.fp, "r") as f:
            for line in f:
                if line.startswith("@"):
                    continue  # Skip header lines
                fields = line.split("\t")
                read_name = fields[0]
                flag = int(fields[1])
                reference_name = fields[2]
                position = int(fields[3])
                mapping_quality = int(fields[4])
                cigar = fields[5]
                mismatch = 0
                try:
                    if reference_name != "*":
                        # Find XM field for mismatches
                        for i in [13, 14, 15, 12, 11]:
                            if fields[i].split(":") == "XM":
                                mismatch = int(fields[i].split(":")[2])
                                break
                except IndexError:
                    pass

                parsed_read = {
                    "read_name": read_name,
                    "flag": flag,
                    "reference_name": reference_name,
                    "position": position,
                    "mapping_quality": mapping_quality,
                    "cigar": cigar,
                    "mismatch": mismatch,
                }
                yield parsed_read

    def parse_contig_lengths(self) -> list:
        lengths = {}
        with open(self.fp, "r") as sam_file:
            for line in sam_file:
                if line.startswith("@"):
                    if line.startswith("@SQ"):
                        contig_name = line.split("\t")[1][3:]
                        contig_length = int(line.split("\t")[2][3:])
                        lengths[contig_name] = contig_length
                else:
                    break

        return lengths

    def write(self):
        pass


class Cov3File(File):
    def __init__(
        self,
        fp: Path,
        bin_name: str,
        mapq_cutoff: int = 5,
        mapl_cutoff: int = 50,
        max_mismatch_ratio: float = 0.03,
    ) -> None:
        super().__init__(fp)
        self.bin_name = bin_name

        self.mapq_cutoff = mapq_cutoff
        self.mapl_cutoff = mapl_cutoff
        self.max_mismatch_ratio = max_mismatch_ratio

        if not (0 <= self.mapq_cutoff <= 30):
            raise ValueError(
                f"MapQ cutoff of {self.mapq_cutoff} is not between 0 and 30"
            )
        if not (30 <= self.mapl_cutoff <= 80):
            raise ValueError(
                f"MapL cutoff of {self.mapl_cutoff} is not between 30 and 80"
            )
        if not (0.01 <= self.max_mismatch_ratio <= 0.3):
            raise ValueError(
                f"Max mismatch ratio of {self.max_mismatch_ratio} is not between 0.01 and 0.30"
            )

        self.min_cov_window = 0.1
        self.min_window_count = 5

    def parse(self):
        with open(self.fp) as f:
            for line in f.readlines():
                fields = line.split(",")
                yield {
                    "log_cov": float(fields[0]),
                    "GC_content": float(fields[1]),
                    "sample": fields[2],
                    "contig": fields[3],
                    "length": int(fields[4]),
                }

    def parse_sample_contig(self):
        with open(self.fp) as f:
            data_dict = {}
            for line in f.readlines():
                fields = line.split(",")
                parsed_line = {
                    "log_cov": float(fields[0]),
                    "GC_content": float(fields[1]),
                    "sample": fields[2],
                    "contig": fields[3],
                    "length": int(fields[4]),
                }

                sample = parsed_line["sample"]
                contig = parsed_line["contig"]

                if (sample, contig) not in data_dict:
                    data_dict[(sample, contig)] = {
                        "sample": sample,
                        "contig": contig,
                        "contig_length": parsed_line["length"],
                        "log_covs": [],
                        "GC_contents": [],
                    }

                data_dict[(sample, contig)]["log_covs"].append(parsed_line["log_cov"])
                data_dict[(sample, contig)]["GC_contents"].append(
                    parsed_line["GC_content"]
                )

            for values in data_dict.values():
                yield values

    def write(self, sams: list, fasta: FastaFile, window_params: dict):
        sam_generators = {sam.fp.stem: sam.parse() for sam in sams}
        next_lines = OrderedDict(
            sorted({name: next(sg, {}) for name, sg in sam_generators.items()}.items())
        )

        with open(self.fp, "w") as f_out:
            for contig_name, seq in fasta.parse():
                contig = Contig(
                    contig_name, seq, fasta.sample, fasta.bin_name, **window_params
                )
                logging.debug(f"Current contig: {contig_name}")

                for name, line in next_lines.items():
                    mut_line = line
                    coverages = {}
                    while True:
                        if not mut_line:
                            break  # Generator is exhausted
                        if mut_line["reference_name"] == "*":
                            next_lines[name] = {}
                            break  # SAM file has no more mapped reads
                        if mut_line["reference_name"] != contig_name:
                            break  # This contig is unmapped by this SAM file
                        if contig.windows:
                            coverages = self.__update_coverages(
                                coverages,
                                mut_line,
                                contig.edge_length,
                                contig.window_step,
                            )

                        mut_line = next(sam_generators[name], {})

                    next_lines[
                        name
                    ] = mut_line  # Instead of updating with every iteration of the while loop
                    if coverages:
                        for info in self.__log_cov_info(
                            contig,
                            coverages,
                            contig.edge_length,
                            contig.window_size,
                            contig.window_step,
                        ):
                            f_out.write(
                                f"{info},{name},{contig_name},{contig.seq_len}\n"
                            )

        logging.debug(f"Finished writing {self.fp}")

    def __update_coverages(
        self, coverages: dict, line: dict, edge_length: int, window_step: int
    ) -> dict:
        mapl = self.calculate_mapl(line["cigar"])
        if (
            line["mapping_quality"] >= self.mapq_cutoff
            and mapl >= self.mapl_cutoff
            and line["mismatch"] <= self.max_mismatch_ratio * mapl
        ):
            start_step = int((line["position"] - 1 - edge_length) / window_step)
            end_step = int((line["position"] - 1 + mapl - edge_length) / window_step)

            if start_step not in coverages:
                coverages[start_step] = 0
            coverages[start_step] += window_step - (
                (line["position"] - 1 - edge_length) % window_step
            )
            if end_step not in coverages:
                coverages[end_step] = 0
            coverages[end_step] += (
                line["position"] - 1 + mapl - edge_length
            ) % window_step

            for step in range(start_step + 1, end_step):
                if step not in coverages:
                    coverages[step] = 0
                coverages[step] += window_step

        return coverages

    def __log_cov_info(
        self,
        contig: Contig,
        coverages: dict,
        edge_length: int,
        window_size: int,
        window_step: int,
    ) -> list:
        first_i = contig.windows[0].start
        last_i = contig.windows[-1].end
        first_step = int((first_i - 1 - edge_length) / window_step)
        last_step = int((last_i - 1 - edge_length) / window_step)

        cov_step = []
        cov_window_sum = 0
        qualified_info = []  # Information to output
        n = 0  # Window index

        for step in range(first_step, last_step):
            if step in coverages.keys():
                cov_step.append(coverages[step])
                cov_window_sum += coverages[step]
            else:
                cov_step.append(0)

            if len(cov_step) == window_size / window_step:
                avg_cov_window = cov_window_sum / window_size
                window = contig.windows[n]
                gc_content = window.gc_content
                n += 1
                cov_window_sum -= cov_step.pop(0)
                if avg_cov_window > self.min_cov_window:
                    log_cov = round(math.log(avg_cov_window) / math.log(2), 4)
                    qualified_info.append(f"{log_cov},{gc_content}")

        if n >= self.min_window_count and len(qualified_info) == n:
            return qualified_info
        return []

    @staticmethod
    def calculate_mapl(cigar: str) -> int:
        operations = []
        current_length = ""

        if cigar == "*":
            return -1

        for char in cigar:
            if char.isdigit():
                current_length += char
            else:
                operations.append((int(current_length), char))
                current_length = ""

        return sum([n for n, c in operations if c == "M" or c == "D"]) - sum(
            [n for n, c in operations if c == "I"]
        )
