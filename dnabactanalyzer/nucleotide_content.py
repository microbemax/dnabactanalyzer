from __future__ import annotations 
from typing import Tuple, List
from scipy import signal

class NucleotideContent:
    """
    Analyzes and visualizes nucleotide content in DNA sequences in wig files.

    This class reads DNA sequences from a FASTA file, calculates GC and AT content, and
    generates WIG files for visualization. It's designed to handle large sequences and apply
    smoothing for clearer analysis of nucleotide distribution.

    Parameters:
    fasta_path (str): Path to the FASTA file with the DNA sequence.
    smooth (int, optional): Smoothing factor for content calculation, defaults to 101.
                             Should be an odd number for best results.

    Usage:
    Create an instance with the path to your FASTA file. Then call `GC_content_wig` or
    `AT_content_wig` with an output path to generate WIG files for GC or AT content.
    """
    def __init__(self, fasta_path: str, smooth: int = 101):
        self.fasta_path: str = fasta_path
        self.smooth: int = smooth
        self.accession: str
        self.sequence: str
        self.accession, self.sequence = self.read_fasta()
        if self.smooth % 2 == 0:
            self.smooth = self.smooth + 1 # make sure the value is odd

    def read_fasta(self) -> Tuple[str, str]:
        with open(self.fasta_path, 'r') as infile:
            lines = infile.read().splitlines()
            return lines[0][1:], lines[1]

    def calculate_content(self, nucleotides: List[str]) -> List[float]:
        content: List[float] = []
        window: int = 100
        for i in range(len(self.sequence)):
            ss: str = self.sequence[i:window+i+1]
            count: float = sum(ss.count(n) for n in nucleotides) / float(window)
            content.append(count)
        return signal.savgol_filter(content, self.smooth, 3)

    def write_wig(self, content: List[float], wig_output: str, content_type: str) -> None:
        header1: str = f'track\ttype=wiggle_0\tname={content_type}_content\tgraphType=points\tvisibility=full\tcolor=168,130,88\n'
        header2: str = f'fixedStep\tchrom={self.accession}\tstart=1\tstep=1\tspan=1\n'
        with open(wig_output, 'w') as f:
            f.write(header1)
            f.write(header2)
            for i in content:
                f.write(f'{i:.1f}\n')

    def GC_content_wig(self, wig_output: str) -> None:
        gc_content: List[float] = self.calculate_content(['G', 'C'])
        self.write_wig(gc_content, wig_output, 'GC')

    def AT_content_wig(self, wig_output: str) -> None:
        at_content: List[float] = self.calculate_content(['A', 'T'])
        self.write_wig(at_content, wig_output, 'AT')
