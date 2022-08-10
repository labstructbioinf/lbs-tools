import gzip, subprocess
from Bio.PDB.DSSP import _make_dssp_dict
import pandas as pd


def parse_dssp_output(dssp_fn: str, use_gzip: bool=False):

	if use_gzip:
		f = gzip.open(dssp_fn, 'rt')
	else:
		f = open(dssp_fn, 'r')

	lines = [line.rstrip() for line in f.readlines()[28:]]

	f.close()
	dssp = {int(line[0:5].strip()): {'pdb_num': line[5:11].strip(), 'pdb_chain': line[11:12].strip(), 
									 'pdb_resn': line[13].strip(), 'pdb_ss': line[16:17]} for line in lines}
	dssp = pd.DataFrame.from_dict(dssp, orient='index')
	return dssp

def run_dssp(pdb_path, dssp_bin='/opt/apps/dssp-3.1.4/bin/mkdssp'):
	dssp_path = f'{pdb_path.rstrip(".pdb")}.dssp'

	cmd = f'{dssp_bin} {pdb_path} > {dssp_path}'
	p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

	stdout, stderr = p.communicate()
	stdout, stderr = stdout.decode("utf-8"), stderr.decode("utf-8")

	assert stderr=="", 	f'Failed to run {dssp_bin} for {pdb_path}:\n\n\n{stderr}'	

	dssp_data = parse_dssp_output(dssp_path)
	return dssp_data 

def get_dssp_seq(fn: str) -> str:
	"""
	Extracts sequence from DSSP output file.
	:param fn: input filename
	:return: sequence in PDB structure
	# TODO: Connectivity info??
	"""
	f = open(fn, 'r')
	out_dict, keys = _make_dssp_dict(f)
	seq = [(out_dict[key][0]) for key in keys]
	f.close()
	return ''.join(seq)





