import numpy as np
import pandas as pd
import sys, math
import seaborn as sns

#Then I used your scripts to compute the rmsds but for the RMSD->score scaling function, I used a sigmoid. I got the base form of a sigmoid from Wikipedia
#1./(1 + numpy.exp(-value))

#I extended this basic function by two parameters like this:
#1./(1 + numpy.exp(-factor*(rmsd_value - shift)))

#I honestly speaking do not know a nice way of finding 'good' (the definition of 'good' is part of the issue) 
#values for the parameters factor and shift (rmsd_value is one input RMSD). The shift adjusts which x-axis value
# (RMSD in our case) corresponds to a function value of 0.5. 'factor' can be used to adjust the steepness of the sigmoid. 
#I played around with a couple of values and found factor=4 and shift=15 to work OK-ish in this case.


def conversion_function_RMSD(value, factor=1.5, shift=0.5): # RMSD-like (lower the better)
    return 1./(1 + math.exp(- factor *(value - shift)))

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

	def __init__(self, df_links, df_desc, pickmin=True, cf=conversion_empty_function):
		"""
		Args:
			df_links: A Pandas DatFrame with id1, id2, value data
			df_desc:  A Pandas DataFrame with id (index!), sequence, group data
			pickmin:  In the case of redundant links (id1->id2 and id2->id1 present) the one
					  with lowest (pickmin=True) or highest (pickmin=Flase) score will be 
					  taken
			cf:		  Function to apply to the scores *after* picking the best one (see above)
			
		"""
		
		self.df_desc = df_desc
		
		self.df_links = df_links
		self.df_links[['id1', 'id2']] = self.df_links[['id1', 'id2']].astype(str)
		
		# for bi-directional hits pick the one with smallest/biggest value
		self.df_links[['id1', 'id2']] = np.sort(self.df_links[['id1', 'id2']], axis=1)
		
		if pickmin:
			self.links = self.df_links.loc[self.df_links.groupby(['id1','id2'])["score"].idxmin()]
		else:
			self.links = self.df_links.loc[self.df_links.groupby(['id1','id2'])["score"].idxmax()]
			
		self.keys = list(set(self.links.id1.tolist() + self.links.id2.tolist()))
		
		self.links['score'] = self.links['score'].apply(cf)		
		self.links=self.links[['id1', 'id2', 'score']].values.astype(str)
		
		self.key2pos = {}
		for pos, k in enumerate(self.keys):
			self.key2pos[str(k)]=pos
			
	def write(self, outfile, groups=None):
		"""
		Writes stored edges as CLANS file
		
		"""
		f = open(outfile, 'w')
		f.write('sequences=%s\n' % len(self.keys))
		f.write(self.header)

		for k in self.keys:
			tmp = self.df_desc.loc[int(k)]
			f.write(">%s [%s]\n" % (tmp['name'], tmp['group']))
			f.write("%s\n" % tmp.seq)

		f.write("</seq>\n<hsp>\n")
		for l in self.links:
			f.write("%s %s:%s\n" % (self.key2pos[l[0]], self.key2pos[l[1]], l[2]))
		f.write("</hsp>")
		
		# Add groups
		nr_groups = self.df_desc.group.value_counts().shape[0]		
		palette = sns.color_palette(None, nr_groups)

		for g, color in zip(self.df_desc.groupby(by='group'), palette):
		
			g_ids = [str(self.key2pos[str(i)]) for i in g[1].index]
		
			f.write("""\n<seqgroups>
name=%s
type=%s
size=6
hide=0
color=%s;%s;%s
numbers=%s;
</seqgroups>""" % (g[0], 
				   0, 
				   int(color[0]*255), 
				   int(color[1]*255),
				   int(color[2]*255), 
				   ";".join(g_ids)
				   ))
		
		f.close()
		
if __name__ == "__main__":
	
	pass
	

