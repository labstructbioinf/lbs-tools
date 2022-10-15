
from csb.bio.io.hhpred import HHOutputParser
import os
import pandas as pd

def parse_allvsall(hhr_dir, df_desc, CLANS_EVAL_CUT=1e-3):

	links = []

	for query_index in df_desc.index.tolist():
	
		
	
		hhr_file = f'{hhr_dir}/{df_desc.loc[query_index].name}.hhr'
		
		#print(hhr_file)

		if not os.path.isfile(hhr_file):
			raise IOError(f'Missing hhr file {hhr_file}')
	
		hits = HHOutputParser(alignments=False).parse_file(hhr_file)

		for hit in hits:
	
			#print(hit.id)
	
			if hit.evalue > CLANS_EVAL_CUT: continue
				
			hit_indexes = df_desc[df_desc.full_name == hit.id].index.tolist()

			assert hit_indexes!=[]
	
			for hit_index in hit_indexes:
				if query_index != hit_index:
					links.append((hit.evalue, int(query_index), int(hit_index)))

	df_links = pd.DataFrame(links, columns = ['score', 'id1', 'id2'])

	return df_links
