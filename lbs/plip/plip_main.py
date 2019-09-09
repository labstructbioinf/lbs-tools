from plip.modules.preparation import PDBComplex
import plip_output_parser as plip_parse
import pandas as pd
import pickle

def convert_interaction_for_df(pdb, inter):
	# in: interaction data
	# out: list ready for df

	# select only data that will not go into columns
	row_data   = ['restype_l', 'reschain_l', 'resnr_l', 'inter_type', 'restype', 'reschain', 'resnr']
	inter_data = { k:v for k,v in inter.items() if k not in row_data }

	row = [pdb] + [ inter[d] for d in row_data ] + [inter_data]

	return row

def run_plip(pdb, pdb_db_path):
	# input: one pdb index, (list of sequences to check - maybe process from all data later)
	# main: load pdb from path
	# out: status of job, save results in file during function run
	# out-format: dict of ligands: keys are ligands ids (name:chain:resnum) and values -
	# lists of dicts, each dict one interaction
	# out-format-alt: dataframe for pdb with ligands as rows and columns:
	# restype_l, reschain_l, resnr_l, inter_type, restype, reschain, resnr, inter_data(maybe without already present columns)

	pdb_filepath        = pdb_db_path + pdb[1:3] + '/' + pdb + '.pdb1'
	ligand_interactions = []

	# PLIP analysis
	mol = PDBComplex()
	# check if file present
	try:
		mol.load_pdb(pdb_filepath)
	except:
		print('=> ' + pdb + ': File not found')
		return False
	ligand_list_from_mol = [ x.hetid for x in mol.ligands ]
	print('=> ' + pdb + ': ' + str(len(ligand_list_from_mol)) + ' ligands to process')
	mol.analyze()

	for bsid in [":".join([x.hetid, x.chain, str(x.position)]) for x in mol.ligands]:
		# for every ligand encountered in pdb

		interactions = mol.interaction_sets[bsid] # Contains all interaction data
		for inter in interactions.all_itypes:
			inter_parsed = plip_parse.parse_inter(inter)
			ligand_interactions.append(convert_interaction_for_df(pdb, inter_parsed))

	# create df for pdb
	df = pd.DataFrame(ligand_interactions, columns=['pdb', 'restype_l', 'reschain_l', 'resnr_l', 'inter_type',
													'restype', 'reschain', 'resnr', 'inter_data'])

	pickle.dump(df, open(pdb + '.p', 'wb')) # dump df
	return True
