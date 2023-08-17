import logging
from abc import ABC, abstractmethod

class Sequence(ABC):
    def __init__(self, seq: str) -> None:
        super().__init__()
        self.seq = seq
        if not self.seq:
            logging.error("No sequence or empty sequence given")
            raise ValueError
    
    def get(self) -> str:
        return self.seq
        
class Contig(Sequence):
    def __init__(self, seq: str, edge_length: int, window_size: int = 5000, window_step: int = 100) -> None:
        super().__init__(seq)
        self.edge_length = edge_length
        self.window_size = window_size
        self.window_step = window_step

        if not (10 <= self.window_step <= 1000):
            logging.error(f"Window step of {self.window_step} is not between 10 and 1,000, please choose a value in this range")
            raise ValueError
        if not (500 <= self.window_size <= 10000):
            logging.error(f"Window size of {self.window_size} is not between 500 and 10,000, please choose a value in this range")
            raise ValueError
        if self.window_size < self.window_step * 2:
            logging.error(f"Window size must be at least twice the window step value")
            raise ValueError
        if self.window_size % self.window_step != 0:
            logging.error(f"Window step must evenly divide window size")
            raise ValueError
    
    def calculate_GC_content(self):
        gc_count = self.seq.count('G') + self.seq.count('C')
        total_count = len(self.seq)
        gc_content = (gc_count / total_count) * 100
        self.gc_content = gc_content
    
    def calculate_GC_skew(self):
        gc_skew = [0]
        for i in range(1, len(self.seq)):
            gc_skew.append(gc_skew[i - 1] + (self.seq[i] == 'G') - (self.seq[i] == 'C'))
        self.gc_skew = gc_skew
