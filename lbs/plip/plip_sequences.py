from Bio import SeqIO, PDB
import re

from Bio import pairwise2
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from Bio.SubsMat import MatrixInfo as matlist
import numpy as np


def get_sequence_biopython(pdbpath, pdb_id):

	''' get sequence from pdb file using Biopython '''
	''' input: path to pdb file, pdb id with chain, example input: '/home/pdb/1ztm.pdb', '1ztm_A' '''
	''' output: sequence (string) of desired chain '''

	sequence = ''
	handle   = open(pdbpath, 'r')

	for record in SeqIO.parse(handle, "pdb-atom"):
		# go through sequence records in pdb and select only that with desired chain id
		if record.id.split(':')[-1] == pdb_id.split('_')[-1]:
			sequence = record.seq

	handle.close()

	return sequence

def get_first_residue_id_dssp(pdbname, pdbpath, pdb_id):

	''' get id of first residue in pdb file (pdb numbering) '''
	''' use dssp data where pdb ids are keys of dssp data dictionary '''
	''' input: pdbname, path to pdb file, pdb id with chain, example input: 'pdb1ztm', '/home/pdb/1ztm.pdb', '1ztm_A' '''
	''' output: integer id of first residue in pdb, example output: 44 '''

	p          = PDB.PDBParser()
	structure  = p.get_structure(pdbname, pdbpath)
	model      = structure[0]
	dssp       = PDB.DSSP(model, pdbpath)
	chain_id   = pdb_id.split('_')[-1] # fetch chain id from pdb_id
	chain_data = []

	for res in list(dssp.keys()):
		# consider only residues of desired chain
		if res[0] == chain_id:
			chain_data.append(res)

	if chain_data:
		return chain_data[0][1][1]
	else:
		return False

def run_mapper(pdb, pdb_db_path, max_gaps=0.2):

	# load sequence from Rossmann data
	seq_ross = pdb[3]

	# find structure and load seq from structure
	pdb_id        = pdb[0]
	pdbpath       = pdb_db_path + pdb_id[1:3] + '/' + pdb_id[:4] + '.pdb1'
	try:
		seq_pdb       = get_sequence_biopython(pdbpath, pdb_id)
		first_res_pdb = get_first_residue_id_dssp('pdb' + pdb_id[:4], pdbpath, pdb_id)
		if not first_res_pdb:
			return(0, pdb[0], 'No chain data')
	except:
		return(0, pdb[0], 'No file or dssp error')

	# print(seq_ross)
	# print(seq_pdb)

	# SeqIO.write(SeqRecord(Seq(seq1, IUPAC.protein), 'seq1', 's1', 's1'), 'seq1.fas', 'fasta')
	# SeqIO.write(SeqRecord(seq2, 'seq2', 's2', 's2'), 'seq2.fas', 'fasta')

	# alignments = pairwise2.align.globalxx(seq_ros[3], seq)
	# print(pairwise2.format_alignment(*alignments[0]))

	# print()
	aln = sorted(pairwise2.align.localdd(seq_pdb, seq_ross, matlist.blosum62, -1000, -1000, -1, 0), key=lambda i:i[2])[0]
	print(aln)


	#print(aln[3], aln[4])


	# find start and end positions of aligned match
	i = re.search('[^-]', aln[1]).start()
	# print(i)
	j = re.search('[^-]', aln[1][::-1]).start()
	# print(j)
	j = len(aln[1])-j
	# print(j)

	# count gaps in aln when go

	def check_if_rest_of_seq_X(aln, pos, pdbseq2rosseq_values):

		# count X in segment of pdb
		X_count = 1
		last_pos = pos+1
		for res in aln[0][pos+1:]:
			if res == 'X':
				X_count += 1
				last_pos += 1
			else:
				break

		# check if mismatch present and if so count all values up to now
		# if there is enough X for all previous residues send mismatch signal
		if 'mismatch' in pdbseq2rosseq_values:
			aln_before = len(pdbseq2rosseq_values)
			if X_count >= aln_before:
				return ('mismatch', last_pos)

		# if problem not before X then maybe after
		# count number of residues left in rosseq
		rosseq_left = len(aln[1][pos+1:].replace('-', ''))

		if X_count >= rosseq_left:
			return ('rest_X', False)
		else:
			return (False, False)

	pdb_idx = first_res_pdb
	#print('start idx:', pdb_idx)
	pdbseq2rosseq  = {}
	j_corrected_to_X = j
	i_corrected_to_X = i
	for pos in range(0, j):
		# if res in pdb == '-': do not add to dict, no change pdb_idx
		if aln[0][pos] == '-':
			pass
		# if res in seq == '-': do not add to dict, increase pdb_idx
		elif aln[1][pos] == '-':
			pdb_idx += 1
			if aln[0][pos] == 'X':
				X_status = check_if_rest_of_seq_X(aln, pos, pdbseq2rosseq.values())
				# if there is X and rest of sequence is all X in pdb: break
				if X_status[0] == 'rest_X':
					j_corrected_to_X = pos
					break
				# if there is X but rest of seq is ok then check what happens before: if mismatch - reset pdbseq2rosseq
				elif X_status[0] == 'mismatch':
						i_corrected_to_X = X_status[1]
						pdbseq2rosseq  = {}

		else:
			#print(pdb_idx)
			#print(aln[1][pos], aln[0][pos])

			if aln[1][pos] != aln[0][pos]:
				pdbseq2rosseq[pdb_idx] = 'mismatch'
			else:
				pdbseq2rosseq[pdb_idx] = aln[0][pos]
			pdb_idx += 1
	#print(pdbseq2rosseq)

	# reject sequences with too many gaps
	print(i_corrected_to_X, j_corrected_to_X)
	if aln[1][i_corrected_to_X:j_corrected_to_X].count('-') / len(aln[1][i_corrected_to_X:j_corrected_to_X]) < max_gaps:
		return (pdb[0], pdbseq2rosseq)
	else:
		return(0, pdb[0], 'Too many gaps')

	# x=np.asarray(list(aln[0][i:j]))
	# y=np.asarray(list(aln[1][i:j]))
	# pos = (y!='-') & (y!='X')
    #
	# # print(pos)
    #
	# # assure that all aligned residues pairs are equal
	# if all(x[pos] == y[pos]):
    #
	# 	# reject sequences with too many gaps
	# 	if 1.*aln[1][i:j].count('-')/len(aln[1][i:j]) < max_gaps:
	# 		temp = i
	# 		seq = aln[0][i:j]
	# 	else:
    #
	# 		#if pdb=='2wic_A':
	# 		#    print(aln[1][i:j])
	# 		#    print(aln[0][i:j])
	# 		#    sys.exit(-1)
	# 		reason="too many gaps in the aln"
	# 		return (0, pdb[0])

	# from Bio.Blast.Applications import NcbiblastpCommandline
	# print(NcbiblastnCommandline(version=1)())
	# blast_cline = NcbiblastpCommandline(query='seq1.fas', subject='seq2.fas', evalue=1, outfmt=6)()[0]
	# print(blast_cline)
