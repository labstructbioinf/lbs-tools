from plip.modules.preparation import PDBComplex
from . import plip_output_parser as plip_parse
import pandas as pd
import pickle, tempfile, gzip, os
from Bio.PDB import *
from Bio.PDB.PDBExceptions import PDBConstructionException

def convert_interaction_for_df(pdb, inter):
	# in: interaction data
	# out: list ready for df

	# select only data that will not go into columns
	row_data   = ['restype_l', 'reschain_l', 'resnr_l', 'inter_type', 'restype', 'reschain', 'resnr']
	inter_data = { k:v for k,v in inter.items() if k not in row_data }

	row = [pdb] + [ inter[d] for d in row_data ] + [inter_data]

	return row
	
	
def get_unzipped_tempfile(gz_fn):
    f = tempfile.NamedTemporaryFile(dir='/tmp/', mode='wt', delete=False)
    fh = open(f.name, 'w')
    fh.write(gzip.open(gz_fn, mode='rt').read())
    fh.close()
    return f.name

def run_plip(pdb, pdb_filepath, dump=True, dumpdir='', resdf=False, max_res_count=4000):
	# input: one pdb index, (list of sequences to check - maybe process from all data later)
	# main: load pdb from path
	# out: status of job, save results in file during function run
	# out-format: dict of ligands: keys are ligands ids (name:chain:resnum) and values -
	# lists of dicts, each dict one interaction
	# out-format-alt: dataframe for pdb with ligands as rows and columns:
	# restype_l, reschain_l, resnr_l, inter_type, restype, reschain, resnr, inter_data(maybe without already present columns)

	if dumpdir!='' and dump:
		assert dumpdir[-1]=='/', "dumpdir must end with '/'"

	picklename = dumpdir + pdb + '.p'

	if dump:
		if os.path.isfile(picklename):
			print('=> ' + pdb + ': in cache')
			if resdf:
				return pd.read_pickle(picklename)
			else:
				return True

	# unzip if pdb is gzipped
	if pdb_filepath[-3:] == '.gz':
		try:
			pdb_filepath=get_unzipped_tempfile(pdb_filepath)
		except FileNotFoundError:
			print('=> ' + pdb + ': Gz file not found')
			return False			
		gz=True
	else:
		gz=False	
		
	# first try to parse the structure with Bio.PDB
	#parser = PDBParser(PERMISSIVE=0)
	#try:
	#structure = parser.get_structure('temp', pdb_filepath)		
	#except PDBConstructionException:
	#	print('=> ' + pdb + ': Problem with parsing')
	#	return False
		
	# check the number of residues
	#res_count = sum([len(list(model.get_residues())) for model in structure])
	
	#if res_count>=max_res_count:
	#	print('=> ' + pdb + ': The input molecule is too big (%d residues)' % res_count)
	#	return False		

	ligand_interactions = []
	
	# PLIP analysis
	mol = PDBComplex()
	try:
		mol.load_pdb(pdb_filepath)
	except:
		print('=> ' + pdb + ': File not found')
		return False
	ligand_list_from_mol = [ x.hetid for x in mol.ligands ]
	print('=> ' + pdb + ': ' + str(len(ligand_list_from_mol)) + ' ligands to process')
	
	if ligand_list_from_mol==0:
		return False
	
	try:
		mol.analyze()
	except TypeError:
		print('=> ' + pdb + ': Problem with parsing')
		return False
		
	for bsid in [":".join([x.hetid, x.chain, str(x.position)]) for x in mol.ligands]:
		# for every ligand encountered in pdb

		interactions = mol.interaction_sets[bsid] # Contains all interaction data
		for inter in interactions.all_itypes:
			inter_parsed = plip_parse.parse_inter(inter)
			ligand_interactions.append(convert_interaction_for_df(pdb, inter_parsed))

	# create df for pdb
	df = pd.DataFrame(ligand_interactions, columns=['pdb', 'restype_l', 'reschain_l', 'resnr_l', 'inter_type',
													'restype', 'reschain', 'resnr', 'inter_data'])
	if dump:
		pickle.dump(df, open(picklename, 'wb')) # dump df
		
	if gz:
		# del tempfile
		os.remove(pdb_filepath)
	
	if resdf:
		return df	
	else:
		return True