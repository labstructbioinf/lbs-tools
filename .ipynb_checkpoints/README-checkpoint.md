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

## md

```python
from lbs.md import Params, OpenMM

params = Params() # MD params object
path_pdb = 'example.pdb'

mm = OpenMM(params)
mm.prepare_components() # initialize force field and langevin integrator
pdb = mm.prepare_pdb(path_pdb) # fix pdb issues & add solvent

df = mm.run(pdb, 'result.pdb') # run simulation store structure in 'result.pdb'
                               # df - energy etc. over time
```