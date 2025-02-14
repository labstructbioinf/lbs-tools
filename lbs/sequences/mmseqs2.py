import tempfile
import os
import pandas as pd

class MMSeqsClusterer:
	def __init__(self, mmseqs_loc='/opt/apps/MMseqs2/bin/mmseqs', tmp_dir='/tmp/', debug=False):
		"""
		Wrapper for clustering protein sequences (fragments) with MMSeqs2 (Steinegger et al. 2017)
		:param mmseqs_loc: location of the MMSeqs2 binary
		:param tmp_dir: location of the temporary directory for intermediate files
		"""
		self.mmseqs_loc = mmseqs_loc
		self.tmp_dir = tmp_dir
		self.debug = debug
	
	def cluster(self, df, min_identity=0.25, coverage=0.5, cov_mode=0, cluster_mode=0):
		"""
		Cluster sequences in the pandas DataFrame ('sequence' column) with MMSeqs2
		Manual at https://github.com/soedinglab/mmseqs2/wiki
		
		:param df: DataFrame storing protein sequences ('sequence' column)
		:param min_identity: minimum identity [0, 1]
		:param coverage: minimum alignment coverage [0, 1]
		:param cov_mode: coverage mode [0-3] (0 coverage mode should be used to cluster full length protein sequences)
		:param cluster_mode: cluster mode [0-2]
		:return: updated DataFrame containing 'clust_id' column denoting the cluster each sequence was assigned to
				 (note that the cluster name denotes its representative sequence's id)
		"""

		fas_fn = self.tmp_dir + next(tempfile._get_candidate_names())
		out_prefix = self.tmp_dir + next(tempfile._get_candidate_names())
	
		with open(fas_fn, 'w') as f:
			for id, data in df.iterrows():
				f.write('>{}\n{}\n'.format(id, data['sequence']))

		if self.debug:
			aux = "-v 1"
		else:
			aux = "-v 0 >/dev/null"

		cmd = f'{self.mmseqs_loc} easy-cluster {fas_fn} {out_prefix} tmp --min-seq-id {min_identity} -c {coverage} --cov-mode {cov_mode} --cluster-mode {cluster_mode} {aux}'
																			 
		status = os.system(cmd)
		if status != 0:
			raise RuntimeError('MMSeqs2 clustering failed! Check input parameters.')
			
		# Parse MMSeqs2 output
		tmp_df = pd.read_csv('{}_cluster.tsv'.format(out_prefix), header=None, sep='\t')
		tmp_df.columns = ['clust_rep', 'clust_entities']
		clust_data = tmp_df.groupby(by='clust_rep').agg({'clust_entities': list}).to_dict(orient='index')
	
		# Remove temporary files
		os.system('rm {}'.format(fas_fn))
		os.system('rm {}*'.format(out_prefix))
		
		cluster_indexes = {entry: representative for representative, entries in clust_data.items() for entry in
						   entries['clust_entities']}
		df = df.copy()
		df['clust_id'] = df.index.map(lambda x: cluster_indexes[x])
		return df
