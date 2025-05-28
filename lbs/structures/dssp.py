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

def run_dssp(pdb_path, dssp_path='', dssp_bin='/opt/apps/dssp-3.1.4/bin/mkdssp'):
    if dssp_path == '':
        dssp_path = f'{pdb_path.rstrip(".pdb")}.dssp'

    cmd = f'{dssp_bin} {pdb_path} > {dssp_path}'
    p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)

    stdout, stderr = p.communicate()
    stdout, stderr = stdout.decode("utf-8"), stderr.decode("utf-8")

    assert stderr == "", f'Failed to run {dssp_bin} for {pdb_path}:\n\n\n{stderr}'

    dssp_data = parse_dssp_output(dssp_path)
    return dssp_data

# def run_dssp(pdb_path, dssp_bin='/opt/apps/dssp-3.1.4/bin/mkdssp'):
#     # Create a temporary file to store the uncompressed PDB content
#     with tempfile.NamedTemporaryFile(delete=False, mode='w') as temp_pdb_file:
#         temp_pdb_path = temp_pdb_file.name
#         
#         # Open the .gz PDB file and decompress it into the temporary file
#         with gzip.open(pdb_path, 'rt') as gz_file:
#             temp_pdb_file.write(gz_file.read())
#     
#     try:
#         # Generate a temporary file to store the DSSP output
#         with tempfile.NamedTemporaryFile(delete=False, mode='w') as temp_dssp_file:
#             dssp_path = temp_dssp_file.name
# 
#         # Run the DSSP command with the temporary file for output
#         cmd = f'{dssp_bin} {temp_pdb_path} > {dssp_path}'
#         p = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
# 
#         stdout, stderr = p.communicate()
#         stdout, stderr = stdout.decode("utf-8"), stderr.decode("utf-8")
# 
#         assert stderr == "", f'Failed to run {dssp_bin} for {pdb_path}:\n\n\n{stderr}'
# 
#         # Parse the DSSP data
#         dssp_data = parse_dssp_output(dssp_path)
# 
#     finally:
#         # Clean up the temporary files after use
#         if os.path.exists(temp_pdb_path):
#             os.remove(temp_pdb_path)
#         if os.path.exists(dssp_path):
#             os.remove(dssp_path)
# 
#     return dssp_data

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





