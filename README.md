# lbs-tools

## sequence
### MMseq2 clustering
```python
from lbs.sequences import MMSeqsClusterer

mmseq = MMSeqsClusterer(mmseqs_loc='path_to_mmseq2_bin', tmp_dir='/tmp/')
```
### list of sequences to fasta file
```python
from lbs.sequences import seq_list_to_fasta

seq_list_to_fasta(sequence_list: List[str], path: str = 'file.fasta')
```
### run psiblast
```python
from lbs.sequences import run_psiblast

run_psiblast(path_psiblast: str, #path to psiblast directory
             fasta_file: str, #input .fasta file
             max_target_seqs: int =2, dbtype: str ='prot')
```