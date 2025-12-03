import os
import tempfile
import subprocess
import pandas as pd
from pathlib import Path
import shutil


class MMSeqsClusterer:
	def __init__(self, mmseqs_loc='/opt/apps/MMseqs2/bin/mmseqs', tmp_dir='/tmp', debug=False, cpus=1):
		self.mmseqs_loc = mmseqs_loc
		self.tmp_dir = tmp_dir
		self.debug = debug
		self.cpus = cpus

	def cluster(self, df, min_identity=0.25, coverage=0.9, cov_mode=0, cluster_mode=0,
				sensitivity=None, kmer_per_seq=None, max_seqs=None, short_seq_threshold=15):
		"""
		Cluster sequences in a DataFrame using MMSeqs2, with automatic support for short sequences.
		"""
		assert 'sequence' in df.columns, "Input DataFrame must contain a 'sequence' column"

		shortest = df['sequence'].map(len).min()
		is_short_mode = shortest < short_seq_threshold

		if is_short_mode and self.debug:
			print(f"[MMSEQS] Shortest sequence has {shortest} aa. Switching to short-sequence mode.")

		sensitivity = sensitivity if sensitivity is not None else (6.0 if is_short_mode else 4.0)
		if is_short_mode:
			kmer_per_seq = kmer_per_seq if kmer_per_seq is not None else 21
			max_seqs = max_seqs if max_seqs is not None else 50

		tmp_root = tempfile.mkdtemp(dir=self.tmp_dir)
		fas_fn = os.path.join(tmp_root, "input.fasta")
		out_prefix = os.path.join(tmp_root, "out")

		with open(fas_fn, 'w') as f:
			for idx, row in df.iterrows():
				f.write(f">{idx}\n{row['sequence']}\n")

		cmd = [
			self.mmseqs_loc, "easy-cluster", fas_fn, out_prefix, os.path.join(tmp_root, "tmp"),
			"--min-seq-id", str(min_identity),
			"-c", str(coverage),
			"--cov-mode", str(cov_mode),
			"--cluster-mode", str(cluster_mode),
			"-s", str(sensitivity),
			"--threads", str(self.cpus),
			"-v", "1"
		]

		if kmer_per_seq is not None:
			cmd.extend(["--kmer-per-seq", str(kmer_per_seq)])
		if max_seqs is not None:
			cmd.extend(["--max-seqs", str(max_seqs)])

		if self.debug:
			print(f"[MMSEQS DEBUG] Running command: {' '.join(cmd)}")
			print(f"[MMSEQS DEBUG] Temp dir: {tmp_root}")

		try:
			subprocess.run(
				cmd,
				check=True,
				stdout=subprocess.PIPE,
				stderr=subprocess.PIPE,
				text=True
			)
		except subprocess.CalledProcessError as e:
			print("MMSeqs2 clustering failed.")
			print("Command:", ' '.join(cmd))
			print("Return code:", e.returncode)
			print("STDOUT:")
			print(e.stdout)
			print("STDERR:")
			print(e.stderr)
			raise RuntimeError("MMSeqs2 clustering failed. See output above for details.")

		cluster_tsv = f"{out_prefix}_cluster.tsv"
		if not os.path.exists(cluster_tsv):
			raise RuntimeError(f"Expected output file not found: {cluster_tsv}")

		tmp_df = pd.read_csv(cluster_tsv, sep='\t', header=None, names=["clust_rep", "clust_entity"])
		tmp_df["clust_rep"] = tmp_df["clust_rep"].astype(str)
		tmp_df["clust_entity"] = tmp_df["clust_entity"].astype(str)
		cluster_map = tmp_df.set_index("clust_entity")["clust_rep"].to_dict()

		df = df.copy()
		df['clust_id'] = df.index.map(lambda x: cluster_map.get(str(x), None))

		if self.debug:
			print(f"[MMSEQS DEBUG] Clustering complete. Found {len(tmp_df['clust_rep'].unique())} clusters.")
			print(f"[MMSEQS DEBUG] Temp files retained at: {tmp_root}")
		else:
			# solidne sprzątanie całego katalogu tymczasowego, łącznie z podkatalogami
			try:
				shutil.rmtree(tmp_root)
			except Exception:
				pass

		return df
