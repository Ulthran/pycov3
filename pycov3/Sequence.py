import logging

class Window():
    def __init__(self, start: int, end: int, sub_seq: str) -> None:
        self.start = start
        self.end = end
        self.len = self.end - self.start

        if self.start < 0:
            logging.error(f"Window start must be larger than 0 ({self.start} > 0)")
            raise ValueError
        if self.len <= 0:
            logging.error(f"Window end must be larger than start ({self.end > self.start})")
            raise ValueError
        
        self.gc_content = self.calculate_GC_content(sub_seq)
        self.gc_skew = self.calculate_GC_skew(sub_seq)

    """
    Returns subsequence represented by window (inclusive on start and exclusive on end)
    """
    def get_window(self, seq: str):
        if self.end > len(seq):
            logging.error(f"Window end must be smaller than length of sequence ({self.end} < {len(seq)})")
            raise ValueError
        return seq[self.start:self.end]
    
    @staticmethod
    def generate_windows(seq_len: int, edge_length: int, window_size: int, window_step: int) -> list:
        return [(x, x + window_size) for x in range(edge_length, seq_len - edge_length - window_size + 1, window_step)]
    
    @staticmethod
    def calculate_GC_content(seq: str) -> float:
        gc_count = seq.count("G") + seq.count("C")
        total_count = len(seq)
        gc_content = gc_count / total_count
        return gc_content

    @staticmethod
    def calculate_GC_skew(seq: str): # CHECK THIS
        gc_skew = [0]
        for i in range(1, len(seq)):
            gc_skew.append(
                gc_skew[i - 1] + (seq[i] == "G") - (seq[i] == "C")
            )
        return gc_skew


class Contig():
    def __init__(
        self,
        name: str,
        seq: str,
        sample: str,
        bin: str,
        edge_length: int,
        window_size: int = 5000,
        window_step: int = 100,
    ) -> None:
        self.name = name
        self.seq_len = len(seq)
        self.sample = sample
        self.bin = bin

        self.edge_length = edge_length
        self.window_size = window_size
        self.window_step = window_step

        if not (10 <= self.window_step <= 1000):
            logging.error(f"Window step of {self.window_step} is not between 10 and 1,000")
            raise ValueError
        if not (500 <= self.window_size <= 10000):
            logging.error(f"Window size of {self.window_size} is not between 500 and 10,000")
            raise ValueError
        if self.window_size < self.window_step * 2:
            logging.error(
                f"Window size must be at least twice the window step value"
            )
            raise ValueError
        if self.window_size % self.window_step != 0:
            logging.error(f"Window step must evenly divide window size")
            raise ValueError
        
        if self.seq_len >= self.window_size + 2 * self.edge_length:
            self.windows = [Window(s, e, seq[s:e]) for s, e in Window.generate_windows(self.seq_len, self.edge_length, self.window_size, self.window_step)]
            logging.debug(f"Num windows: {len(self.windows)}")
            import sys
            #logging.debug(f"Mem of windows: {sum([sys.getsizeof(i) for i in self.windows])}")
        else:
            #logging.debug(f"Contig {self.name} is too short and will be ignored")
            self.windows = []
