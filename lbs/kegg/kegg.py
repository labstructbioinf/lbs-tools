import os, sys
sys.path.append('..')
from utils import memoized
from inspect import getmembers, isfunction
from Bio.KEGG import REST

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
def get_pathways():
	pathways = REST.kegg_list('pathway').read()
	pathways = [i.split('\t') for i in pathways.split('\n')]
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
    seq = list(seq.split('\n'))
    seq = ''.join(seq[1:-1])
    return seq



def configure(cachedir='/tmp', verbose=False):
	for f in dir(sys.modules[__name__]):

		if f[:4] == 'get_':
			exec(f+".cachedir = cachedir")
			exec(f+".verb = verbose")
			print(f'function {f} has been configured')
