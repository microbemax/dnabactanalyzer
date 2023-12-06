import pytest
from dnabactanalyzer.genbank_file_generator import GenBankProcessor
import pandas as pd
import os

@pytest.fixture
def processor():
    test_gbk_file = 'test_data/test.gbk'
    return GenBankProcessor(test_gbk_file)

def test_write_full_seq(processor):
    output_path = 'test_data/test_output.fasta'
    processor.write_full_seq(output_path)
    assert os.path.exists(output_path)
    os.remove(output_path)  # Cleanup

def test_get_full_info_df_genome(processor):
    df = processor.get_full_info_df(record_type='genome')
    assert isinstance(df, pd.DataFrame)
    assert all(isinstance(col, str) for col in df.columns)  # Testing column types

def test_get_full_info_df_plasmid(processor):
    df = processor.get_full_info_df(record_type='plasmid')
    assert isinstance(df, pd.DataFrame)
    assert all(isinstance(col, str) for col in df.columns)  # Testing column types
