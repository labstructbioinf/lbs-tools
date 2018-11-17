import time
from Bio import Entrez, SeqIO
import numpy as np

Entrez.email = "s.dunin-horkawicz@cent.uw.edu.pl"

def getseqs(accessions, chunksize=1000, format='fasta'):
	"""
	Download protein sequences from NCBI
	
	
	Arguments:
	
		accessions: list of NCBI accession numbers
		chunksize: number of proteins to be downloaded at once
		format: "fasta" or "gb"
	
	Returns:
		a list of protein sequences
	
	
	"""

	assert format in ['fasta', 'gb']

	accessions = np.array(accessions)

	l = np.array_split(accessions, max(len(accessions)/chunksize, 1))

	result = []
		
	for gis_sub in l:

		search_results = Entrez.read(Entrez.epost("protein", id=",".join(gis_sub)))
		webenv = search_results["WebEnv"]
		query_key = search_results["QueryKey"] 
		fetch_handle = Entrez.efetch(db="protein", rettype=format, webenv=webenv, query_key=query_key)	
	
		if format=='fasta':
			res = list(SeqIO.FastaIO.FastaIterator(fetch_handle))
		elif format=='gb':
			res = list(SeqIO.InsdcIO.GenBankIterator(fetch_handle))
	
		result.extend(res)
	
		# let's not make NCBI angry
		time.sleep(1)

	return result		

if __name__ == '__main__':
	pass
	
