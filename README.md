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

### coloring sequence

```python
from lbs.utils.color_sequence import DivColorScaling
from IPython.display import HTML, display
clr = DivColorScaling()
# change color map
clr.cmap = plt.cm.Greens

sequence = "LBS is the best!"
# use float 0-1 range or int 0-255
seq_importance = [0, 125, 130, 140, 150, 200, 0, 50, 250, 100, 100, 100, 60, 60, 60, 60, 60]
# iterate over sequence and importance
html_string = ''
for letter, importance in zip(sequence, seq_importance):
    html_string += clr.html_colored_letter(letter, importance)
# it may not be viewed correctly :)
HTML(html_string)
```
<body>
<span style= "color:rgb(247, 252, 245)">L</span><span style= "color:rgb(119, 197, 120)">B</span><span style= "color:rgb(112, 194, 116)">S</span><span style= "color:rgb(96, 186, 108)"> </span><span style= "color:rgb(80, 178, 100)">i</span><span style= "color:rgb(25, 130, 62)">s</span><span style= "color:rgb(247, 252, 245)"> </span><span style= "color:rgb(211, 238, 205)">t</span><span style= "color:rgb(0, 74, 29)">h</span><span style= "color:rgb(154, 214, 149)">e</span><span style= "color:rgb(154, 214, 149)"> </span><span style= "color:rgb(154, 214, 149)">b</span><span style= "color:rgb(202, 234, 195)">e</span><span style= "color:rgb(202, 234, 195)">s</span><span style= "color:rgb(202, 234, 195)">t</span><span style= "color:rgb(202, 234, 195)">!</span>
</body>



### scripts

#### calculate protein sequence embeddings
Script to create embeddings from sequences via prot_t5_xl_half_uniref50 stored in dataframe
by default `seq` column in used as embedder input. Records are stored
as list maintaining dataframe order. 
Example use
```bash
python scripts/embeddings.py -i data.csv -o data.emb
```
In python load via:
```python 
import torch
torch.load(..)
```
or
```python 
import pickle
with open(.., 'rb') as f:
    embs = pickle.load(f)
```
