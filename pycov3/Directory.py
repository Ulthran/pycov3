import logging
from abc import ABC, abstractmethod
from pathlib import Path

from .File import FastaFile, SamFile

class Directory(ABC):
    def __init__(self, fp: Path, overwrite: bool) -> None:
        super().__init__()
        self.fp = fp.resolve()
        self.overwrite = overwrite

        if not fp.exists():
            logging.info(f"{self.fp} does not exist, creating it now")
            self.fp.mkdir(parents=True, exist_ok=True)
        if not fp.is_dir():
            logging.error(f"{self.fp} is not a directory")
            raise ValueError
        if any(self.fp.iterdir()) and not self.overwrite:
            logging.error(f"{self.fp} is a non-empty directory, please either point output to an empty or non-existent directory or run with the overwrite flag")
            raise ValueError
        
    
    
class FastaDir(Directory):
    def __init__(self, fp: Path, overwrite: bool, coverage_params: dict) -> None:
        super().__init__(fp, overwrite)

        self.fastas = [FastaFile(x, coverage_params) for x in self.fp.iterdir() if x.endswith((".fasta", ".fa", ".fna"))]
        if not self.fastas:
            logging.error(f"No files found ending in .fasta, .fa, or .fna in {self.fp}")
            raise ValueError
        
class SamDir(Directory):
    def __init__(self, fp: Path, overwrite: bool) -> None:
        super().__init__(fp, overwrite)

        self.sams = [SamFile(x) for x in self.fp.iterdir() if x.endswith(".sam")]
        if not self.sams:
            logging.error(f"No files found ending in .sam in {self.fp}")
            raise ValueError
    
    def calculate_edge_length(self) -> int:
        return 0