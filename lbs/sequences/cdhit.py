import os, tempfile, subprocess
from Bio import SeqIO

# example

'''
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

sequences = [SeqRecord(Seq(row.sseq), id=row.sseqid) for _, row in res_df.iterrows()]

clust = cdhit.cdhit(sequences, cdhitbin='/opt/apps/cd-hit/cd-hit')
clusters = clust.run(identity=0.95)
'''


class cdhit:
	def __init__(self, sequences, cdhitbin="", cpus=1):
		assert os.path.isfile(cdhitbin), 'please provide full path to the cd-hit binary'
		self.cdhitbin = cdhitbin
		self.sequences = sequences
		self.cpus=cpus
	
	def run(self, identity=0.9, coverage=0):
		assert identity>=0.4 and identity<=1.0
		
		# define word length
		if identity>0.7:
			word=5
		elif identity>0.6:
			word=4
		elif identity>0.5:
			word=3
		else:
			word=2
		
		print ('Word lendth set to', word)	

		temp_idents = []

		# prepare input fasta file
		inf = tempfile.NamedTemporaryFile(dir='/tmp/', mode='wt', delete=False)
		maxdesclen = 0
		for s in self.sequences:
			ident = s.id
			assert ident.find('...')==-1 and ident.find('>')==-1, 'invalid sequence name'
			assert ident not in temp_idents, 'duplicate seq ids'
			temp_idents.append(ident)
			maxdesclen = max([maxdesclen, len(ident)])
			inf.write(s.format('fasta'))
		inf.close()
		
		print ('Maximal description length', maxdesclen)
	
		# prepare empty output file
		outf = tempfile.NamedTemporaryFile(dir='/tmp/', mode='wt', delete=False)
		outf.close()
	
		cmd = "%s -i %s -o %s -d %s -T %s -c %s -n %s -A %s" % (self.cdhitbin, inf.name, outf.name, maxdesclen+1, self.cpus,
														  identity, word, coverage)
	
		status = subprocess.call(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
		
		assert status==0, 'non-clean exit from cdhit. check in %s and out %s files' % (inf.name, outf.name)
		
		os.remove(inf.name)
		
		# parse fasta results
		self.results_fasta = list(SeqIO.parse(outf.name, 'fasta'))
		os.remove(outf.name)
		
		# parse cluster assignments
		results_clusters={}
		for l in open(outf.name + '.clstr'):
			if l[0]=='>':
				cluster = int(l[8:-1])
			else:
				ident = l[l.find('>')+1:l.find('...')]
				if not cluster in results_clusters:
					results_clusters[cluster]=[ident]
				else:
					results_clusters[cluster].append(ident)
				
		os.remove(outf.name + '.clstr')	
		
		return results_clusters		

	
		
		
	
	