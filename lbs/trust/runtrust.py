import sys
import os.path
from multiprocessing import Pool

from Bio import SeqIO


def run_TRUST(fasta, binary="~/apps/trust"):

	binary = binary.rstrip("/")

	gi = fasta.description

	if not os.path.isfile(str(gi) + '.txt'):
		print(gi)
		tfname = '%s_temp.fas' % gi
		SeqIO.write(fasta, open(tfname, 'w'), "fasta")
		cmd = f'java -Xmx2700m -cp {binary} nl.vu.cs.align.SelfSimilarity -matrix {binary}/BLOSUM62 -fasta {tfname} -noseg -gapo 8 -gapx 2 > {gi}.txt' 
		os.system(cmd)
		os.remove(tfname)
	else:
		pass
		

if __name__ == "__main__":
	if len(sys.argv) != 2:
		print("Usage: python script.py <input_fasta_file>")
		sys.exit(1)

	infile = sys.argv[1]

	MyPool = Pool(4)
	sequences = SeqIO.parse(infile, 'fasta')
	_ = list(MyPool.map(run_TRUST, sequences))
	MyPool.close()
	MyPool.join()
