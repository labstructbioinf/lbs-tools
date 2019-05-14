import numpy as np
import pandas as pd

#Then I used your scripts to compute the rmsds but for the RMSD->score scaling function, I used a sigmoid. I got the base form of a sigmoid from Wikipedia
#1./(1 + numpy.exp(-value))

#I extended this basic function by two parameters like this:
#1./(1 + numpy.exp(-factor*(rmsd_value - shift)))

#I honestly speaking do not know a nice way of finding 'good' (the definition of 'good' is part of the issue) 
#values for the parameters factor and shift (rmsd_value is one input RMSD). The shift adjusts which x-axis value
# (RMSD in our case) corresponds to a function value of 0.5. 'factor' can be used to adjust the steepness of the sigmoid. 
#I played around with a couple of values and found factor=4 and shift=15 to work OK-ish in this case.


def conversion_function1(value, factor=1.5, shift=2): # RMSD-like (lower the better)
    return 1./(1 + math.exp(- factor *(value - shit)))

def conversion_function2(value): #  (higher the better)
    return math.exp(-1*(value/1 ))

def conversion_empty_function(value):
	return value

class clanswriter():


	header = """<param>
maxmove=0.1
pval=1
usescval=false
complexatt=true
cooling=1.0
currcool=1.0
attfactor=10.0
attvalpow=1
repfactor=10.0
repvalpow=1
dampening=1.0
minattract=1.0
cluster2d=false
blastpath=X
formatdbpath=X
showinfo=false
zoom=1.0
dotsize=2
ovalsize=10
groupsize=4
usefoldchange=false
avgfoldchange=false
colorcutoffs=0.0;0.1;0.2;0.3;0.4;0.5;0.6;0.7;0.8;0.9;
colorarr=(230;230;230):(207;207;207):(184;184;184):(161;161;161):(138;138;138):(115;115;115):(92;92;92):(69;69;69):(46;46;46):(23;23;23):
</param>
<rotmtx>
1.0;0.0;0.0;
0.0;1.0;0.0;
0.0;0.0;1.0;
</rotmtx>
<seq>
"""

	def __init__(self, df, pickmin=True, cf=conversion_empty_function):
		"""
		df - a pandas df with id1, id2, value data
		"""
		self.df = df
		self.keys = list(set(df.id1.tolist() + df.id2.tolist()))
		# for bi-directional hits pick the one with smallest/biggest value
		self.df[['id1', 'id2']] = np.sort(self.df[['id1', 'id2']], axis=1)
		if pickmin:
			self.links = self.df.loc[self.df.groupby(['id1','id2'])["evalue"].idxmin()]
		else:
			self.links = self.df.loc[self.df.groupby(['id1','id2'])["evalue"].idxmax()]
		self.links['evalue'] = self.links['evalue'].apply(cf)
		self.links=self.links.values
			
	def write(self, outfile):
		f = open(outfile, 'w')
		f.write('sequences=%s\n' % len(self.keys))
		f.write(self.header)
		key2pos = {}
		for pos, k in enumerate(self.keys):
			key2pos[str(k)]=pos
			f.write(">%s\n" % k)
			# FIXME!
			f.write("X\n")
		f.write("</seq>\n<hsp>\n")
		for l in self.links:
			f.write("%s %s:%s\n" % (key2pos[l[0]], key2pos[l[1]], l[2]))
		f.write("</hsp>")
		f.close()
		
if __name__ == "__main__":
	
	# Let's prepare some data using BLAST
	# makeblastdb -in PF00672_seed.txt -dbtype prot
	# blastp -db PF00672_seed.txt -query PF00672_seed.txt -outfmt "6 qseqid sseqid evalue" -out PF00672_seed.BLAST.csv -num_threads 4

	df = pd.read_csv('/Users/sdh/Documents/projekty/tab2clans/PF00672_seed.BLAST.csv', sep='\t', names=['id1', 'id2', 'evalue'])
	df = df[df.evalue <= 1e-3]
	cw = clanswriter(df)
	cw.write('temp.clans')

