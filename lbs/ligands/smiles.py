from rdkit import Chem, DataStructs 
import pandas as pd
import itertools

class ligand:

	def __init__(self, name, smiles):
		self.name = name
		self.smiles = smiles
		self.mol = Chem.MolFromSmiles(self.smiles)
		self.fp = Chem.RDKFingerprint(self.mol)

	def __repr__(self):
		return(self.name)
	
	def compare(self, l2):
		return DataStructs.TanimotoSimilarity(self.fp,l2.fp)

if __name__ == "__main__":

	# SMILES codes obtained from http://ligand-expo.rcsb.org/ld-download.html
	smilesdf = pd.read_csv('./data/Components-smiles-cactvs.smi.tar.gz', names=['SMILES', 'ligand', 'name'], sep='\t')
	smilesdf = smilesdf.dropna()
	
	assert smilesdf.ligand.is_unique
	smilesdf.set_index('ligand', inplace=True)

	ligands = [ligand(l_name, smilesdf.loc[l_name].SMILES) for l_name in ['NAD', 'NDP', 'SAM']]
	
	for comb in itertools.combinations(ligands, 2):
		l1, l2 = comb
		print(comb, l1.compare(l2))