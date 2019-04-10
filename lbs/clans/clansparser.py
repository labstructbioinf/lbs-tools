# CLANS runfile classes

import sys, re, numpy, random

debug = False

def makeid(a,b):
        return str(max((a), (b)))+"-"+str(min((a), (b)))

class clansException(Exception):
	def __init__(self, comment):
		self.comment = comment
	def __str__(self):
		return repr(self.comment)


class clansParseException(clansException): pass

# to powinna byc pochodna Fasta.Record
class clansseq:
	def __init__(self, pos, gi, sequence, defline):
		self.pos = pos
		self.gi = gi
		self.sequence = sequence.replace("\n", "")
		self.title = defline
		self.incluster = False
		self.cluster = None 
		self.connections = []
		
	def clustername(self):
		if self.cluster != None:
			return self.cluster.name
		else:
			return "UNKNOWN"

	def sortedconnections(self):
		return sorted(self.connections, key=lambda x: x.value)

	def __str__(self):
		return ">%s\n%s" % (self.title, self.sequence)

class clanscluster:
	def __init__(self, name, ctype, size, hide, color, cseqs, nrseq):
	
		for seq in cseqs:
			seq.cluster = self
			
	
		self.name = name
		self.type = ctype
		self.size = size # dot (other shape) size!
		self.nrseq = nrseq
		self.hide = hide
		self.color = color
		self.seqs = cseqs
		
	def addseq(self, newseq):
		self.seqs.append(newseq)
		self.nrseq += 1

class connection:
	def __init__(self, id1, id2, value):
		self.value = value
		self.id1 = id1
		self.id2 = id2


class clansrun:
	def __init__(self, runfile, onlywarn=False, readHSP=True, HSPcut=10, accparser=lambda x: x, warn=True):




		def readandstore(store):
			line = runfile.readline()
			store+=line 
			return line, store

		self.pos2gi = {}
		self.gi2pos = {}
		self.pos2seq = {}
		self.connections = {}
		self.clusters = {}	
		
		
		counter = 0

		line, self.header = readandstore("")

		while not line=="<seq>\n":
			line, self.header = readandstore(self.header)	
			
		line = runfile.readline()	

		self.header = self.header[:self.header[:-1].rfind('\n')]

		
	
		while line and not line=="</seq>\n":
			if line[0]==">":
				defline = line[1:] # without ">"
				
				gi = accparser(defline)
								
				if not self.gi2pos.has_key(gi):
					self.gi2pos[gi] = str(counter)		
				else:
					self.gi2pos[gi] = self.gi2pos[gi] + ";" + str(counter)
				#print gi2pos[gi]
				self.pos2gi[str(counter)] = gi
				counter = counter + 1
			
			else:
				#print line
				self.pos2seq[counter-1]=clansseq(counter-1, gi, line, defline[:-1])

			
			line = runfile.readline()


		while not line=="<hsp>\n":

			line = runfile.readline()

		line = runfile.readline()

		if readHSP:

			while line and not line=="</hsp>\n":
			
				line = line[:-2]
			
				if line.count(' ')>1:
					#print line
					line = line[:line.find(' ', line.find(' ')+1)]
					#print "'%s'" % line

				#print line
			

				f,t = line.split(' ')
				t,v = t.split(':')
	
				f,t,v = int(f),int(t), float(v)


				if v <= HSPcut:

					fseq = self.pos2seq[f]
					tseq = self.pos2seq[t]

					c = connection(fseq, tseq, v)

					ident = makeid(f,t)

					if self.connections.has_key(ident):
						if self.connections[ident].value != c.value:
							#pass					
							print ("error!")
							sys.exit(-1)

					else:
						self.connections[ident] = c
						fseq.connections.append(c)
						tseq.connections.append(c)

			
				line = runfile.readline()


		runfile.seek(0)


#0 name#1 type#2 size#3 hide# 4 color# 5 numbers=3371;2318;2122;2077;1891;1845;906;719;690;591;311;

		
		hitre = re.compile("<seqgroups>(.*?)</seqgroups>", re.DOTALL)
		
		rec = re.compile("name=(.*)\ntype=(.*)\nsize=(.*)\nhide=(.*)\ncolor=(.*)\nnumbers=(.*)\n")

		hits = hitre.findall(runfile.read())
		if debug:
			print ("clansrun: %s CLANS seqgroup(s) were/was found" % len(hits))
		
		
		if len(hits)>0: 
			clu = rec.findall("".join(hits))
		else:
			clu = [] # no clusters defined
			

		if debug:
			print ("clansrun:", len(clu), "clusters were read")

			
		
		counter = 0
		dupcontrol = {}	
		for c in clu:
		
			clustername = c[0]
			#print clustername
			while self.clusters.has_key(clustername):
				#print "error in CLANS run file, duplicate cluster name %s" % (c[0])
				clustername = clustername+"_"+str(counter)
				counter += 1
				#sys.exit(-1)

			seqsincluster = []


			# zbieranie sekwencji do klastra

			# linkowanie klaster -> sek
			
			for numb in c[5].split(';'):
				if not numb=='':	
					self.pos2seq[int(numb)].incluster = True	
					seqsincluster.append(self.pos2seq[int(numb)])

			self.clusters[clustername]=clanscluster(clustername,c[1],c[2],c[3],c[4],seqsincluster, len(seqsincluster))
			
			# linkowanie sekwencji sekwencja->klaster
			
			for s in self.clusters[clustername].seqs:
				s.cluster = self.clusters[clustername]

			for numb in c[5].split(';'):
				if not numb=='':		
					if dupcontrol.has_key(numb):
						if warn:
							print ("seq nr %s already in cluster %s, cannot add to culster %s" % (numb, dupcontrol[numb], clustername))
						if not onlywarn:
							sys.exit(-1)
					else:
						dupcontrol[numb]=clustername

		if debug:
			print ("clansrun:", len(self.clusters), "final clusters were formed")
	# http://www.macs.hw.ac.uk/~pdw/topology/ScaleFree.html
	def degree(self):
		counts = []
		for s in self.pos2seq.values():
			c = len(s.connections)
			assert c != 0
			counts.append(c)

		counts, bins = numpy.histogram(counts, range=(1,161), bins=40, normed=False)

		for b, c in zip(bins, counts):
			#print "%s" % (round(c/float(len(self.seqs.keys())),4))
			print ("%s\t%s" % (b,c))

	#def save(self):
	#	print self.header
	#	print "<seq>"
	#	for s in self.pos2seq.values():
	#		print ">"+s.title
	#		print s.sequence


	


	
if __name__ == '__main__':
	test = clansrun(open(sys.argv[1]), onlywarn=True)

	#print len(test.connections.keys())

	for c in test.connections.values():
		print ("test", c.value)

	#test.save()


				

