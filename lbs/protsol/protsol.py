### Methods:

# Local:
# https://protein-sol.manchester.ac.uk/ 
# https://loschmidt.chemi.muni.cz/soluprot/
# https://github.com/mahan-fcb/DSResSol
# https://github.com/jcchan23/GraphSol (many dependencies) 

# Only server:
# https://www-cohsoftware.ch.cam.ac.uk/


import tempfile, os, glob, shutil
import pandas as pd

class predictor():

	def __init__(self, bindir, debug=False):
		self.bindir = bindir
		assert os.path.isdir(self.bindir), 'binary path not found %s' % bindir
		self.debug = debug
	
	def predict(self):
		assert 'sequence' in self.seq_df.columns, "seq_df must contain a `sequence` column"
		
	def write_fasta(self, file):
		f = open(file, 'wt')
		for idx, row in self.seq_df.iterrows():
			f.write(f'>{idx}\n{row.sequence}\n')
		f.close()
		
class soluProt(predictor):

	def __init__(self, bindir, condash, debug=False):
		self.condash = condash
		assert os.path.isfile(condash)
		super().__init__(bindir, debug=debug)
			
	def predict(self, seq_df):
		self.seq_df = seq_df
		super().predict()
		
		infile = tempfile.NamedTemporaryFile(delete=False)
		infile.close()
		
		outfile = tempfile.NamedTemporaryFile(delete=False)
		outfile.close()
		
		self.write_fasta(infile.name)
	
		cmd = f'/bin/bash -c "source {self.condash} && conda activate soluprot && python {os.path.join(self.bindir, "soluprot.py")} --i_fa {infile.name} --o_csv {outfile.name} --tmp_dir /tmp && conda deactivate"'
		os.system(cmd)
		
		res_df = pd.read_csv(outfile.name, skiprows=1, names = ['c', 'id', 'score'])
		res_df.set_index('id', inplace=True)
		res_df.drop(columns=['c'], inplace=True)
		
		os.remove(infile.name)
		os.remove(outfile.name)
		
		return res_df
		
class protSol(predictor):

	def predict(self, seq_df):
		self.seq_df = seq_df
		super().predict()
		
		res = []
		
		with tempfile.TemporaryDirectory() as tmpdirname:
			if self.debug: print(f'{tmpdirname} created')
			for ext in ['pl', 'sh', 'txt']:
				for file in glob.glob(os.path.join(self.bindir, f'*.{ext}')):
					shutil.copy(file, tmpdirname)
		
			cwd = os.getcwd()
			os.chdir(tmpdirname)
			
			self.write_fasta('input.fasta')
			os.system('./multiple_prediction_wrapper_export.sh input.fasta')
				
			id_set = set([str(i) for i in self.seq_df.index])
			for l in open('seq_prediction.txt').readlines():
				if l.startswith('SEQUENCE PREDICTIONS'):				
					row = l.split(',')
					ident, score = row[1][1:], float(row[3])
					assert ident in id_set, f'{ident} not found in df index'
					res.append((ident, score))			
			os.chdir(cwd)
			
			
		res_df = pd.DataFrame(res, columns=['id', 'score'])
		res_df.set_index('id', inplace=True)
		return res_df
	
				
			
if __name__ == "__main__":

	df=pd.DataFrame(['MRIAVVDGQGGGIGKALVEKLKQELPDIHVIALGTNALATSLMLKAGANEGASGENAIVFNAPRVQIIAGGIGIIGANSMLGELTPLMAKAISESPARKVLIPFNRCHLHVVGVKPVSLPELVEEAVSCIKQLCGEV',
				     'MRITVIDGQGGGIGKAITEKLRKNLPKETEIIALGTNALATALMLKAGANDGASGENAVVRSVAETDVIIGSLAIIVPHAMMGELTPAMAEAIAMSKAPKILLPLNRCGIDVVGVYAEPLPHLIDHILKKVKDLMEK'],
				columns=['sequence'])

	p = protSol(bindir='/home/users/sdunin/apps/protsol', debug=True)
	print(p.predict(df))
	
	p = soluProt(bindir='/home/users/sdunin/apps/soluprot-1.0.1.0', debug=True,
				 condash='/home/users/sdunin/miniconda3/etc/profile.d/conda.sh')
	print(p.predict(df))
	
	