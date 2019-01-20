import pandas as pd
import numpy as np

class taxdb:
	def __init__(self, ncbidir):
		if ncbidir[-1]=='/': ncbidir=ncbidir[:-1]
		self.ncbidir = ncbidir
		
		"""
	tax_id								-- node id in GenBank taxonomy database
 	parent tax_id						-- parent node id in GenBank taxonomy database
 	rank								-- rank of this node (superkingdom, kingdom, ...) 
 	embl code							-- locus-name prefix; not unique
 	division id							-- see division.dmp file
 	inherited div flag  (1 or 0)		-- 1 if node inherits division from parent
 	genetic code id						-- see gencode.dmp file
 	inherited GC  flag  (1 or 0)		-- 1 if node inherits genetic code from parent
 	mitochondrial genetic code id		-- see gencode.dmp file
 	inherited MGC flag  (1 or 0)		-- 1 if node inherits mitochondrial gencode from parent
 	GenBank hidden flag (1 or 0)        -- 1 if name is suppressed in GenBank entry lineage
 	hidden subtree root flag (1 or 0)   -- 1 if this subtree has no sequence data yet
 	comments							-- free-text comments and citations
		"""
		
		self.nodes = pd.read_csv(ncbidir+'/nodes.dmp', sep='|', index_col=False, 
																dtype={'taxid': np.int, 'parent_taxid':np.int} ,
															     names=['taxid',
																		'parent_taxid',
																		'rank',
																		'embl_code',
																		'divid',
																		'divflag', 
																		'gencode',
																		'gencode_flag',
																		'mitgencode',
																		'mitgencode_flag',
																		'gb_flag',
																		'subtree_flag',
																		'comments'])
																		
		for fix in ['embl_code', 'rank', 'comments']:															
			self.nodes[fix] = self.nodes[fix].apply(lambda x:x.replace('\t', '').replace('\f', ''))
	
	def keepgenomes(self, genome_taxids):
		"""
		from self.nodes remove species (rank) that are in genome_taxids list
		"""
		
		self.nodes.drop(self.nodes.index[(self.nodes['rank'] == 'species') & (~self.nodes['taxid'].isin(genome_taxids))], inplace=True)

	def getchildren(self, taxid):
		"""
		get all leafs (nodes without children) that are descendants of the given taxid
		"""

		children = self.nodes[self.nodes.parent_taxid == taxid].taxid.tolist()

		if len(children)==0:
			#if self.nodes[self.nodes.taxid==taxid]['rank'].iloc[0] == 'species':
			yield taxid
		else:
			for child in children:
				yield from self.getchildren(child)
		
		
		
	