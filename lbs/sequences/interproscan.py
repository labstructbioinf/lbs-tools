import time
import requests
import pandas as pd

class InterProScan:
	"""
	A class to interact with the InterProScan REST API.
	"""

	BASE_URL = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5"

	def __init__(self, email, analyses=None):
		"""
		Initialize the InterProScan instance.

		:param email: User email required for job submission.
		:param analyses: List of analysis tools to use. Defaults to ["SMART", "TMHMM"].
		"""
		self.email = email
		self.analyses = analyses if analyses else ["SMART", "TMHMM"]
		self.sequences = ""
		self.all_results = []
		self.pid2domains = {}

	def submit_job(self, sequences):
		"""
		Submit a job to the InterProScan API.

		:param sequences: The protein sequences in FASTA format.
		:return: Job ID if successful, None otherwise.
		"""
		payload = {
			'email': self.email,
			'title': 'Batch InterProScan Job',
			'sequence': sequences,
			'appl': ','.join(self.analyses),
			'goterms': 'false',
			'pathways': 'false'
		}
		response = requests.post(f"{self.BASE_URL}/run/", data=payload)
		if response.status_code == 200:
			return response.text.strip()
		print(f"Error submitting job: {response.text}")
		return None

	def check_status(self, job_id):
		"""
		Check the status of a submitted job.

		:param job_id: The job ID to check.
		:return: True if job is finished, False otherwise.
		"""
		while True:
			response = requests.get(f"{self.BASE_URL}/status/{job_id}")
			if response.status_code == 200:
				status = response.text.strip()
				if status == 'FINISHED':
					return True
				elif status in ['RUNNING', 'PENDING', 'QUEUED']:
					print(f"Job {job_id} is {status}...")
					time.sleep(10)
				else:
					print(f"Job {job_id} failed with status: {status}")
					return False
			else:
				print(f"Error checking status: {response.text}")
				return False

	def get_results(self, job_id):
		"""
		Retrieve results for a completed job.

		:param job_id: The job ID.
		:return: JSON response containing the results, or None if failed.
		"""
		response = requests.get(f"{self.BASE_URL}/result/{job_id}/json")
		if response.status_code == 200:
			return response.json()
		print(f"Error retrieving results: {response.text}")
		return None

	def read_fasta(self, file_path):
		"""
		Read sequences from a FASTA file.

		:param file_path: Path to the FASTA file.
		"""
		with open(file_path, 'r') as file:
			self.sequences = file.read()

	def read_dict(self, seq_dict):
		"""
		Read sequences from a dictionary.

		:param seq_dict: Dictionary where keys are sequence identifiers and values are sequences.
		"""
		self.sequences = ''.join([f'>{ident}\n{seq}\n' for ident, seq in seq_dict.items()])

	def scan(self, batch_size=50):
		"""
		Submit sequences in batches to InterProScan and collect results.

		:param batch_size: Number of sequences per batch.
		"""
		sequence_batches = self.sequences.strip().split('>')[1:]

		for i in range(0, len(sequence_batches), batch_size):
			batch = '>' + '>'.join(sequence_batches[i:i + batch_size])
			print(f"Submitting batch {i // batch_size + 1}")
			job_id = self.submit_job(batch)
			if job_id:
				print(f"Job submitted successfully. Job ID: {job_id}")
				if self.check_status(job_id):
					results = self.get_results(job_id)
					if results:
						print(f"Results for batch {i // batch_size + 1} are ready!")
						self.all_results.append(results)
				else:
					print(f"Job {job_id} did not complete successfully.")
			else:
				print("Failed to submit job.")
			time.sleep(5)

	def parse_result(self, results):
		"""
		Parse the InterProScan results and return a pandas DataFrame.

		:param results: JSON results from InterProScan.
		:return: DataFrame containing domain annotations.
		"""
		data = []
		matches = results.get('matches', [])

		for match in matches:
			signature = match.get('signature', {})
			locations = match.get('locations', [])

			for location in locations:
				data.append({
					"signature_accession": signature.get('accession', ''),
					"signature_name": signature.get('name', ''),
					"signature_description": signature.get('description', ''),
					"signature_library": signature.get('signatureLibraryRelease', {}).get('library', ''),
					"signature_version": signature.get('signatureLibraryRelease', {}).get('version', ''),
					"model_ac": match.get('model-ac', ''),
					"start": location.get('start', ''),
					"end": location.get('end', ''),
					"evalue": location.get('evalue', ''),
					"score": location.get('score', ''),
					"hmmStart": location.get('hmmStart', ''),
					"hmmEnd": location.get('hmmEnd', ''),
					"hmmLength": location.get('hmmLength', '')
				})
				
		# properly handle the situation then matches are empty (yet we still want a df to be returned)
		return pd.DataFrame(data, columns=[
			"signature_accession", "signature_name", "signature_description", 
			"signature_library", "signature_version", "model_ac", 
			"start", "end", "evalue", "score", "hmmStart", "hmmEnd", "hmmLength"
		]).sort_values(by='start')


	def parse_all_results(self):
		"""
		Parse all collected results and store them in a dictionary.
		"""
		for batch in self.all_results:
			for result in batch['results']:
				pid = result['xref'][0]['id']
				if pid in self.pid2domains:
					raise ValueError(f"Duplicate protein ID detected: {pid}")
				self.pid2domains[pid] = self.parse_result(result)
