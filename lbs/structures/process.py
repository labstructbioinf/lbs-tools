from biopandas.pdb import PandasPdb
import numpy as np

def merge_chains(pdb1, pdb2):

	"""
	takes pdb1: AF2 model with two chains A and B and stores
	stores pdb2: the same model but with a single chain A and renumbered residues 
	"""

	model = PandasPdb().read_pdb(pdb1)
	# keep only ATOMs
	model._df = {'ATOM':model.df['ATOM']}
	model_df = model.df['ATOM']

	# there should be 2 unique chains
	assert len(model_df.chain_id.unique())==2, pdb1

	# model_df = model_df[model_df.chain_id.isin(['A', 'B'])]

	# Calculate res_count
	res_count = model_df[model_df.chain_id == 'A']['residue_number'].max()

	# Iterate through rows where chain_id is 'B'
	for res_id, g in model_df[model_df.chain_id == 'B'].groupby('residue_number'):
		model_df.loc[g.index, 'residue_number'] = res_id + res_count

	ca_residue_numbers = model_df[model_df.atom_name == 'CA']['residue_number']

	unique_values, counts = np.unique(ca_residue_numbers, return_counts=True)
	duplicates = unique_values[counts > 1]

	assert ca_residue_numbers.is_unique, (pdb1, res_count, duplicates)
	assert ca_residue_numbers.is_monotonic_increasing, pdb1

	model_df.chain_id = 'A'

	model.to_pdb(path=pdb2, 
			records=None, 
			gz=False, 
			append_newline=True)