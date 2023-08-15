from abc import ABC, abstractmethod
from itertools import groupby
from pathlib import Path

class File(ABC):
    def __init__(self, fp: Path) -> None:
        self.fp = fp
    
    @abstractmethod
    def parse(self) -> None:
        pass

    @abstractmethod
    def write(self) -> None:
        pass


class FastaFile(File):
    def __init__(self, fp: Path) -> None:
        super().__init__(fp)

    def parse(self) -> list:
        with open(self.fp) as f:
            faiter = (x[1] for x in groupby(f, lambda line: line[0] == ">"))

            for header in faiter:
                # drop the ">"
                header_str = header.__next__()[1:].strip()

                # join all sequence lines to one.
                seq_str = "".join(s.strip() for s in faiter.__next__())

                yield (header_str, seq_str)

    def write(self):
        pass