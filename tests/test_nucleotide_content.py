import pytest
from dnabactanalyzer.nucleotide_content import NucleotideContent
from scipy import signal
import os
import numpy as np
path = 'test_data/test.fasta'
output_file = 'test_data/test_out.wig'
def test_read_fasta():
    # Test the fasta file reading functionality
    nc = NucleotideContent(path)
    assert nc.accession is not None
    assert nc.sequence is not None

def test_calculate_content_GC():
    # Test GC content calculation
    nc = NucleotideContent(path)
    gc_content = nc.calculate_content(['G', 'C'])
    assert gc_content is not None
    assert isinstance(gc_content, np.ndarray)
    assert all(isinstance(item, float) for item in gc_content)

def test_calculate_content_AT():
    # Test AT content calculation
    nc = NucleotideContent(path)
    at_content = nc.calculate_content(['A', 'T'])
    assert at_content is not None
    assert isinstance(at_content, np.ndarray)
    assert all(isinstance(item, float) for item in at_content)

def test_smooth_parameter():
    # Test the impact of the smooth parameter
    nc = NucleotideContent(path, smooth=100)
    content = nc.calculate_content(['G', 'C'])
    assert content is not None
    assert len(content) == len(nc.sequence)

def test_write_wig():
    # Test writing to WIG file
    nc = NucleotideContent(path)
    content = nc.calculate_content(['G', 'C'])
    nc.write_wig(content, output_file, 'GC')
    assert os.path.exists(output_file)
    os.remove(output_file)  # Clean up after test

def test_GC_content_wig():
    # Test GC_content_wig functionality
    nc = NucleotideContent(path)
    nc.GC_content_wig(output_file)
    assert os.path.exists(output_file)
    os.remove(output_file)  # Clean up after test

def test_AT_content_wig():
    # Test AT_content_wig functionality
    nc = NucleotideContent(path)
    nc.AT_content_wig(output_file)
    assert os.path.exists(output_file)
    os.remove(output_file)  # Clean up after test
