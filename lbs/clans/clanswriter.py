import sys, math
import numpy as np
import pandas as pd
from sklearn import preprocessing
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
cluster2d=true
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
		
		assert df_desc.index.is_unique
		
		self.df_desc = df_desc.copy()
		self.df_desc.index = df_desc.index.map(str)
		
		self.df_links = df_links.copy()
		self.df_links[['id1', 'id2']] = self.df_links[['id1', 'id2']].astype(str)
		
		# for bi-directional hits pick the one with smallest/biggest value
		self.df_links[['id1', 'id2']] = np.sort(self.df_links[['id1', 'id2']], axis=1)
		
		if pickmin:
			self.links = self.df_links.loc[self.df_links.groupby(['id1','id2'])["score"].idxmin()]
		else:
			self.links = self.df_links.loc[self.df_links.groupby(['id1','id2'])["score"].idxmax()]
			
		#self.keys = list(set(self.links.id1.tolist() + self.links.id2.tolist()))
		#assert len(set(self.keys) - set(self.df_desc.index.tolist()))==0
		
		self.keys = self.df_desc.index.tolist()
		assert set(self.links.id1.tolist() + self.links.id2.tolist()) <= set(self.keys)
		
		self.links['score'] = self.links['score'].apply(cf)		
		self.links=self.links[['id1', 'id2', 'score']].values.astype(str)
		
		self.key2pos = {}
		for pos, k in enumerate(self.keys):
			self.key2pos[k]=pos
			

			
	def write(self, outfile, min_group_size=1, color_palette='hot', sort_key=str, sort_reverse=False):
		"""
		Writes stored edges as CLANS file
		
		
		sort_key and sort_reverse: 
		
		"""
		f = open(outfile, 'w')
		f.write('sequences=%s\n' % len(self.keys))
		f.write(self.header)

		for k in self.keys:
			tmp = self.df_desc.loc[k]
			f.write(">%s {%s} [%s]\n" % (tmp['name'], tmp['full_name'], tmp['group']))
			f.write("%s\n" % tmp.seq)

		f.write("</seq>\n<hsp>\n")
		for l in self.links:
			f.write("%s %s:%s\n" % (self.key2pos[l[0]], self.key2pos[l[1]], l[2]))
		f.write("</hsp>")
		
		# Add groups
		groups = self.df_desc.group.unique().tolist()
		sorted_groups = sorted(groups, key=sort_key, reverse=sort_reverse)
		
		print('groups order', sorted_groups)
		
		#nr_groups = self.df_desc.group.value_counts().shape[0]		
		palette = sns.color_palette(color_palette, len(groups))

		groupid2seqgroup = {}

		#for g, color in zip(self.df_desc.groupby(by='group'), palette):
		
		for g_pos, g_name in enumerate(sorted_groups):
		
			color = palette[g_pos]	
			g_ids = [str(self.key2pos[i]) for i in self.df_desc[self.df_desc.group==g_name].index]
			assert len(g_ids)>0
			if len(g_ids) < min_group_size: continue
		
			groupid2seqgroup[g_name] = """\n<seqgroups>
name=%s
type=%s
size=6
hide=0
color=%s;%s;%s
numbers=%s;
</seqgroups>""" % (g_name, 
				   0, 
				   int(color[0]*255), 
				   int(color[1]*255),
				   int(color[2]*255), 
				   ";".join(g_ids)
				   )
				   		
		for g in sorted_groups:
			f.write(groupid2seqgroup[g])
		
		f.close()
		
def gen_df_desc(df_links, cores, group_label="", full_name_label=""):
    
    # create desc df and fill with sequences
    ids = list(set(df_links.id1.tolist() + df_links.id2.tolist()))
    df_desc = pd.DataFrame(ids, columns=['id'])
    df_desc.set_index(df_desc.id, inplace=True)
    df_desc['seq'] = df_desc.index.map(lambda x:cores.iloc[x]['sequence'])

    
    def f(x):
        #d = preds.loc[cores.iloc[x].ecod][['pdb_chain', 'x_name', 'h_name', 't_name', 'f_name', 'source']] # , 'maxscore_bin'
        
        d = cores.iloc[x]
        
        # Define groups based on ECOD
        return pd.Series({"name":f'{d.name}_{d.pdb_chain}', 
                          "group": d[group_label],
                          "full_name":d[full_name_label]})
    
        # Define groups based on bined scores
        #return pd.Series({"name":d.pdb_chain, "group":d.maxscore_bin})

    df_desc=df_desc.merge(df_desc.id.apply(lambda x:f(x)), left_index=True, right_index=True)
    df_desc=df_desc.drop(['id'], axis=1)
    return df_desc
    
    
		
def matrix2clans(matrix, data, percentile=75, group_label="", full_name_label=""):

	"""
	
	Args:
		matrix:     similarity matrix (0 to 1)
		percentile: percentile of the best scores to consider
	
	"""

	tril = np.tril_indices(matrix.shape[0])
	tril_array = np.column_stack(tril)

	CUT_OFF = np.percentile(matrix[tril], percentile)
	print('CUT_OFF', CUT_OFF)
	tril_array_best = tril_array[matrix[tril]  >= CUT_OFF,:]
	tril_array_best_tup = (tril_array_best[:,0], tril_array_best[:,1])
		
	min_max_scaler = preprocessing.MinMaxScaler()

	# Add and scale scores
	df_links = pd.DataFrame(matrix[tril_array_best_tup], columns=['score_raw'])
	df_links['score'] = min_max_scaler.fit_transform(df_links['score_raw'].to_numpy().reshape(-1, 1))
	# Add hits 
	df_links = pd.concat([df_links, pd.DataFrame(tril_array_best, columns=['id1', 'id2'])], axis=1)
	# Remove self-hits
	df_links=df_links[df_links.id1!=df_links.id2]

	df_desc = gen_df_desc(df_links, data, group_label=group_label, full_name_label=full_name_label)

	return df_links, df_desc 






















		
		
		
if __name__ == "__main__":
	
	pass
	

