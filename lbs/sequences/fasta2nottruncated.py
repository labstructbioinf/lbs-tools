
from Bio import SeqIO
import sys

def removeTruncated(filename, outfilename, lrange, lfraction, rrange, rfraction, format='fasta'):

	def gapfraction(s):
		return round(1. * s.count('-') / len(s), 2)
		
	of = open(outfilename, 'w')
	
	removed, ok = 0, 0

	for s in SeqIO.parse(filename, format):
		l = s.seq[:lrange]
		r = s.seq[-rrange:]
		
		lg = gapfraction(l)
		rg = gapfraction(r)
		
		if lg <= lfraction and rg <= rfraction:
			of.write(s.format(format))
			ok+=1
			status=''
		else:
			removed+=1
			status='REMOVED'
			
		print ("{} {:>5} ....... {} {:>5} {}".format(l, lg, r, rg, status))

	of.close()
	
	print ('{}/{} sequences removed'.format(removed, removed+ok))

if __name__ == "__main__":
	removeTruncated(sys.argv[1], sys.argv[2], 20, 0.7, 10, 0.7)
	
	
	