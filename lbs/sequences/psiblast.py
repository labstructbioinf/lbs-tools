import os
import tempfile
import subprocess
from shutil import copy

import pandas as pd

def run_psiblast(path_psiblast, fasta_file, max_target_seqs=2, dbtype='prot'):
    '''
    run psiblast on fasta file return dataframe with:
        `'qid','sid','ident', 'cov','evalue'` columns
    params:
        path_psiblast: str psiblast directory
        fasta_file: str input fasta file
        max_target_seqs: int 
        dbtype: str = prot
    return:
        blast_data: pd.DataFrame
    '''
    assert os.path.isfile(fasta_file), 'invalid `fasta_file` no such file'
    assert os.path.isdir(path_psiblast), '`path_psiblast` must be valid directory'
    cwd = os.getcwd()
    with tempfile.TemporaryDirectory(dir=cwd, suffix='_psiblast') as tmp_dir:
        tmp_fasta = os.path.abspath(tmp_dir + '/data.fasta')
        tmp_df = os.path.abspath(tmp_dir + '/blast_output.csv')
        # copy fasta file to temp_dir
        copy(fasta_file, tmp_fasta)
        # define commands
        cmd_makeblastdb = [
            f'{path_psiblast}/makeblastdb',
                '-dbtype', f'{dbtype}',
                '-in',  f'{tmp_fasta}'
        ]
        cmd_psiblast = [
            f'{path_psiblast}/psiblast',
                '-query', f'{tmp_fasta}',
                '-db',  f'{tmp_fasta}', 
                '-outfmt', '"6 qseqid sseqid pident qcovs evalue qcovs"',
                '-evalue', '1e-1',
                '-num_threads', '4',
                '-max_target_seqs', f'{max_target_seqs}',
                '-out', f'{tmp_df}'
        ]
        # join arguments - passing list of args as in documnetation dont work properly
        cmd_makeblastdb = ' '.join(cmd_makeblastdb)
        cmd_psiblast = ' '.join(cmd_psiblast)
        handle = subprocess.Popen(cmd_makeblastdb, shell=True,
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        try:
            output, errors = handle.communicate(timeout=30)
            print(output, end='\r')
        except TimeoutExpired:
            raise BaseException(errors)
        finally:
            if errors:
                raise BaseException(errors)
        handle = subprocess.Popen(cmd_psiblast, shell=True, 
                                  stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        try:
            output, errors = handle.communicate(timeout=30)
            print(output, end='\r')
        except TimeoutExpired:
            raise BaseException(errors)
        data = pd.read_csv(tmp_df, sep='\t', names=['qid','sid','ident', 'cov','evalue'])
    return data