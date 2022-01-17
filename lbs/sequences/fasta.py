import os
from typing import List


def seq_list_to_fasta(sequence_list: List[str], path: str = 'file.fasta'):
    '''
    convert list of sequences to fasta formatted file
    params:
        sequnence_list List[str] list of sequence strings
        path str output file, should end with .fasta
    return
        None
    '''
    assert path.endswith('.fasta'), '`path` must have `.fasta` extension'
    with open(path, 'wt') as file:
        for i, sequence in enumerate(sequence_list):
            sequence_block = f'>SEQUENCE_{i}\n{sequence}\n'
            file.write(sequence_block)
            