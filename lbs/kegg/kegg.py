import os, sys
sys.path.append('..')
from utils import memoized
from inspect import getmembers, isfunction
from Bio.KEGG import REST
from itertools import chain
import numpy as np

def fix(i):
    assert isinstance(i, list)
    if i[-1] == ['']: 
        return i[:-1]
    else:
        return i

@memoized.cache
def get_organisms():
	organisms = REST.kegg_list('organism').read()
	organisms = [i.split('\t') for i in organisms.split('\n')]
	return organisms

@memoized.cache
def get_pathways(organism):
	pathways = REST.kegg_list('pathway', organism).read()
	pathways = [i.split('\t')[0] for i in pathways.split('\n')]
	return pathways

@memoized.cache
def get_kos(pathway):
	kos = REST.kegg_link('ko', pathway).read()
	kos = [i.split('\t') for i in kos.split('\n')]
	return kos

@memoized.cache  
def get_genes(ko):
	genes = REST.kegg_link('genes', ko).read()
	genes = [i.split('\t') for i in genes.split('\n')]
	return genes
	
@memoized.cache	
def get_seq(gene):
    seq = REST.kegg_get(gene, 'aaseq').read()
    # seq = list(seq.split('\n'))
    # seq = ''.join(seq[1:-1])
    return seq

@memoized.cache
def get_chunks_of_genes(df):
	all_genes_from_df = []
	for org, data in df.iterrows():
		temp = list(chain(*data))
		temp = [org + ":" + i for i in temp]
		all_genes_from_df.extend(temp)
	all_genes_from_df = np.array_split(all_genes_from_df, max(len(all_genes_from_df)/8, 1))
	# Converting np.arrays to lists
	for i in range(len(all_genes_from_df)):
		all_genes_from_df[i] = list(all_genes_from_df[i])
	return all_genes_from_df

@memoized.cache
def genes_to_sequences(list_of_genes):
    genes_and_sequences = {}
    seqs = get_seq(list_of_genes).split('>')
    seqs = [seq.split('(A)')[1].replace('\n', '') for seq in seqs[1:]]
    genes_and_sequences = dict(zip(list_of_genes, seqs))
    return genes_and_sequences


def configure(cachedir='/tmp', verbose=False):
	for f in dir(sys.modules[__name__]):

		if f[:4] == 'get_':
			exec(f+".cachedir = cachedir")
			exec(f+".verb = verbose")
			print(f'function {f} has been configured')
