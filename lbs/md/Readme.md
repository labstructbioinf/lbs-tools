## openmm example

#### minimal
```python
from lbs.md import Params, OpenMM

path_pdb = ... #pdb file location
params = Params()

mm = OpenMM(params)
mm.prepare_components()
pdb = mm.prepare_pdb(path_pdb)
df = mm.run(pdb, 'temp.pdb')
```


#### adaptive box size & pdb repair
```python
from lbs.md import Params, OpenMM, calculate_box_size, Patcher

path_pdb = ... #pdb file location
params = Params()
mm = OpenMM(params)
mm.prepare_components()
# calculate adaptive box wall size
params.boxSize = calculate_box_size(path_pdb, factor=1.5)
# patcher requires PDBFixer installed
patch = Patcher(mm.forcefield)
pdb, rm_residues = patch.read_and_repair(path_pdb)
# rm_residues are list of objects to remove from orginal file
# such as DNA, lignads etc. some custom residues may fall in it
# watch out!
pdb = mm.prepare_pdb(pdb, rm_residues)
df = mm.run(pdb, 'temp.pdb')
```