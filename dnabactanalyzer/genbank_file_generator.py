from __future__ import annotations 
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature
import pandas as pd

class GenBankProcessor:
    """
    A class to process GenBank files and extract useful information.

    Attributes:
    -----------
    gbk_file_path : str
        Path to the GenBank file.

    Methods:
    --------
    _get_features(record: SeqRecord) -> pd.DataFrame:
        Extracts and returns various genomic features from a SeqRecord object.

    write_full_seq(out_path: str) -> None:
        Writes the full nucleotide sequence from the GenBank file to a specified output path in FASTA format.

    get_full_info_df(record_type: str = 'genome') -> pd.DataFrame:
        Parses the GenBank file and compiles a DataFrame of features from either genome or plasmid records.
    """

    def __init__(self, gbk_file_path: str):
        self.gbk_file_path: str = gbk_file_path

    def _get_features(self, record: SeqRecord) -> pd.DataFrame:
        data = {'locus_tag': [], 'old_locus_tag': [], 'gene': [], 'function': [],
                'start': [], 'stop': [], 'strand': [], 'nts': [], 'aa': []}

        for feat in record.features:
            if feat.type in ['CDS', 'rRNA', 'tRNA']:
                for key in data.keys():
                    try:
                        if key in ['start', 'stop', 'strand']:
                            value = getattr(feat.location, key)
                        elif key in ['nts', 'aa']:
                            seq = feat.location.extract(record.seq)
                            value = str(seq) if key == 'nts' else str(seq.translate())
                        else:
                            value = feat.qualifiers.get(key, ['na'])[0]
                    except (AttributeError, IndexError):
                        value = 'na'
                    data[key].append(value)

        gene_df = pd.DataFrame(data)
        #gene_df['gene_length'] = gene_df['stop'] - gene_df['start']
        return gene_df

    def write_full_seq(self, out_path: str) -> None:
        record = SeqIO.read(self.gbk_file_path, 'gb')
        with open(out_path, 'w') as f:
            f.write(f'>{record.name}\n{str(record.seq)}')

    def write_gene_fna(self, out_path: str) -> None:
        pass

    def write_gene_faa(self, out_path: str) -> None:
        pass


    def get_full_info_df(self, record_type: str = 'genome') -> pd.DataFrame:
        dfs = []
        for record in SeqIO.parse(self.gbk_file_path, 'genbank'):
            if (record_type == 'genome' and not record.id.startswith('plasmid')) \
               or (record_type == 'plasmid' and record.id.startswith('plasmid')):
                dfs.append(self._get_features(record))

        return pd.concat(dfs) if dfs else pd.DataFrame()

if __name__ == '__main__':
# Usage Example

    gbk = '/Users/mahmoudal-bassam/Documents/genomes/Streptomyces_venezuelae_NRRL_B-65442.gb'
    outfile = '/Users/mahmoudal-bassam/Documents/test.fasta'
    processor = GenBankProcessor(gbk)
    genome_df = processor.get_full_info_df(record_type='genome')