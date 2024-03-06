import sys
import re
import numpy
import random

def make_id(a, b):
	return str(max(a, b)) + "-" + str(min(a, b))

class ClansException(Exception):
	def __init__(self, comment):
		self.comment = comment
	
	def __str__(self):
		return repr(self.comment)

class ClansParseException(ClansException):
	pass

class ClansSeq:
	def __init__(self, pos, gi, sequence, defline):
		self.pos = pos
		self.gi = gi
		self.sequence = sequence.replace("\n", "")
		self.title = defline
		self.in_cluster = False
		self.cluster = None
		self.connections = []

	def cluster_name(self):
		if self.cluster:
			return self.cluster.name
		else:
			return "UNKNOWN"

	def sorted_connections(self):
		return sorted(self.connections, key=lambda x: x.value)

	def __str__(self):
		return ">%s\n%s" % (self.title, self.sequence)

class ClansCluster:
	def __init__(self, name, ctype, size, hide, color, cseqs, nrseq):
		for seq in cseqs:
			seq.cluster = self
	
		self.name = name
		self.type = ctype
		self.size = size
		self.nrseq = nrseq
		self.hide = hide
		self.color = color
		self.seqs = cseqs

	def add_seq(self, new_seq):
		self.seqs.append(new_seq)
		self.nrseq += 1

class Connection:
	def __init__(self, id1, id2, value):
		self.value = value
		self.id1 = id1
		self.id2 = id2

class ClansRun:
	def __init__(self, file_name, only_warn=False, read_hsp=True, hsp_cut=10, acc_parser=lambda x: x, warn=True, debug = False):

		run_file = open(file_name, 'rt')


		def read_and_store(store):
			line = run_file.readline()
			store += line 
			return line, store

		self.pos2gi = {}
		self.gi2pos = {}
		self.pos2seq = {}
		self.connections = {}
		self.clusters = {}
		counter = 0
		line, self.header = read_and_store("")

		while line != "<seq>\n":
			line, self.header = read_and_store(self.header)	
		
		line = run_file.readline()	
		self.header = self.header[:self.header[:-1].rfind('\n')]

		while line and line != "</seq>\n":
			if line[0] == ">":
				defline = line[1:].strip('\n') # without ">"
				gi = acc_parser(defline)
			
				if gi not in self.gi2pos:
					self.gi2pos[gi] = str(counter)		
				else:
					self.gi2pos[gi] += ";" + str(counter)
			
				self.pos2gi[str(counter)] = gi
				counter += 1
			else:
				self.pos2seq[counter-1] = ClansSeq(counter-1, gi, line, defline[:-1])
			line = run_file.readline()

		while line != "<hsp>\n":
			line = run_file.readline()

		line = run_file.readline()

		if read_hsp:
			while line and line != "</hsp>\n":
				line = line[:-2]
			
				if line.count(' ') > 1:
					line = line[:line.find(' ', line.find(' ')+1)]
			
				f, t = line.split(' ')
				t, v = t.split(':')
				f, t, v = int(f), int(t), float(v)

				if v <= hsp_cut:
					f_seq = self.pos2seq[f]
					t_seq = self.pos2seq[t]
					c = Connection(f_seq, t_seq, v)
					ident = make_id(f, t)

					if ident in self.connections:
						if self.connections[ident].value != c.value:
							pass
					else:
						self.connections[ident] = c
						f_seq.connections.append(c)
						t_seq.connections.append(c)
				line = run_file.readline()

		run_file.seek(0)
		hit_re = re.compile("<seqgroups>(.*?)</seqgroups>", re.DOTALL)
		rec = re.compile("name=(.*)\ntype=(.*)\nsize=(.*)\nhide=(.*)\ncolor=(.*)\nnumbers=(.*)\n")
		hits = hit_re.findall(run_file.read())

		if debug:
			print("ClansRun: %s CLANS seqgroup(s) were/was found" % len(hits))
		
		if len(hits) > 0: 
			clu = rec.findall("".join(hits))
		else:
			clu = [] 

		if debug:
			print("ClansRun:", len(clu), "clusters were read")
	
		counter = 0
		dup_control = {}	
	
		for c in clu:
			cluster_name = c[0]

			while cluster_name in self.clusters:
				cluster_name += "_" + str(counter)
				counter += 1

			seqs_in_cluster = []

			for numb in c[5].split(';'):
				if numb != '':	
					self.pos2seq[int(numb)].in_cluster = True	
					seqs_in_cluster.append(self.pos2seq[int(numb)])

			self.clusters[cluster_name] = ClansCluster(cluster_name, c[1], c[2], c[3], c[4], seqs_in_cluster, len(seqs_in_cluster))

			for s in self.clusters[cluster_name].seqs:
				s.cluster = self.clusters[cluster_name]

			for numb in c[5].split(';'):
				if numb != '':		
					if numb in dup_control:
						if warn:
							print("seq nr %s already in cluster %s, cannot add to culster %s" % (numb, dup_control[numb], cluster_name))
						if not only_warn:
							sys.exit(-1)
					else:
						dup_control[numb] = cluster_name

		if debug:
			print("ClansRun:", len(self.clusters), "final clusters were formed")
		
		run_file.close()
