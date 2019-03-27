import tempfile, sys, math, numpy, subprocess, os, itertools
import os.path
import numpy as np

from sympy import nsimplify
from Bio.PDB.PDBParser import PDBParser

ahelix_p = 3.62705
ahelix_d = 1.51
ahelix_R1 = 2.26

#4hpq
pihelix_p = 4.3985
# Wychodzi 1.1689 (helanal), spore stdev - 0.3423
pihelix_d = 1.1689
	
def inTheCore_new(ph1, w1, rep_len):
	#ph1+=5
	
	#w1 = 102.8 # 7/2
	#w1 = 98.18181818181819 # 11/3
	#rep_len = 11
	
	#if ph1<0:
	#	ph1 = 360 + ph1
	
	def adj(ang):
		ang += 180
	
		if ang>360:
			c = int(ang/360)
			return ang - (360*c) - 180 
		else:
			return ang - 180
	
	fullrepeat = [adj(ph1+(w1*i)) for i in range(rep_len)]
		
	#print fullrepeat
	
	return sorted(fullrepeat, key=lambda x:abs(x))[0]
	
	#sys.exit(-1)
	
	
def inTheCore(ph1, rep_len):
	
	#w1 = 360.0 / (1.0 * rep_len / rep_turn)
	w1 = 360.0 / rep_len 
	

	l=[ph1]

	p=ph1
	while True:
		p=p+w1
		if p>=180: break
		l.append(p)

	p=ph1
	while True:
		p=p-w1
		if p<=-180: break
		l.append(p)
	
	assert len(l)==rep_len
	temp = sorted(l, key=lambda x:abs(x))
	#print temp
	return temp[0]
	
def inTheCore_pandas(row): 
	return inTheCore_new(row['ph1'], row['w1'], row['rep_len'])
	#return inTheCore(row['ph1'], row['rep_len'])

def calcPER(P, tolerance = 0.01):
	f = nsimplify(P, tolerance=tolerance)
	return f.p, f.q

def w0_P_to_p(w0, P):
	return P / (1 + P * w0 / 360)

def w0_to_P(w0, residues_per_turn):
	return residues_per_turn / ( 1 - (residues_per_turn * w0 / 360) )
	
def P_to_w0(P, residues_per_turn):
	return  (- residues_per_turn / P + 1) / residues_per_turn * 360

def calcphi(heptad, phi1):

	heporder = 'dgcfbea'
	heptadlen   = 360/7.0

	distto_f = 3 - heporder.index(heptad)   
	phi=phi1-heptadlen*distto_f  

	return phi+9 # -26 -- 26 range
	
	

def make_bb(f, outname, chainlen, chains  = ['A', 'B', 'C', 'D', 'E']):
	
    fpdb = open(outname + '.pdb', 'w')
    
    atom_id = 1
    res_id  = 1
     
    for coords in f:

		coords = numpy.round(coords, 3)
		fpdb.write('ATOM  ' + '{0: >5}'.format(str(atom_id)) + '{0: >4}'.format('CA') + ' {0: >4}'.format('GLY') + ' ' + chains[(res_id-1) / chainlen] + 
				   '{0: >4}'.format(str(res_id)) + '    ' + '{0: >8}'.format(coords[0]) + '{0: >8}'.format(coords[1]) + '{0: >8}'.format(coords[2]) + 
				   '{0: >6}'.format('1.00') + '{0: >6}'.format('0.00') + '\n')
		res_id  += 1
		atom_id += 4

    fpdb.close()



def CCCPgenerate(outname, cN, chL, R0, R1, w0, w1, A, ph1, cr, dph0, Zoff, Zoff_type):
	"""
	dla zadanych parmetrow uruchamia CCCP + odbudowa backbone BBQ
	"""

	tf = tempfile.NamedTemporaryFile(delete=False)
	tf.write("generateCrickBB(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, struct(\'%s\', 1))" % (cN, chL, R0, R1, w0, w1, A, ph1, cr, dph0, Zoff, Zoff_type))
	tf.close()
	
	cmd = "octave-cli -q %s" % (tf.name)
	
	curdir = os.getcwd()
	os.chdir('/Users/sdh/apps/cccp')
	
	p = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
	p.wait()
	stdout, stderr = p.communicate()	
	
	if stderr!='': 
		print('cos poszlo nie tak...')
		sys.exit(-1)
	
	os.chdir(curdir)
	
	coords = numpy.array([i.split() for i in stdout.split('\n')[2:-2]], dtype=numpy.dtype(float) )
		
	make_bb(coords, outname, chL)
	
	cmd = 'java -classpath /Users/sdh/apps/BBQ/:/Users/sdh/apps/BBQ/jbcl.jar BBQ -d=q_50_xyz.dat -r=./' + str(outname) + '.pdb > ./' + str(outname) + '.bbq.pdb'
	p = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
	p.wait()
	
	f = open('temp.seq', 'w')
	print >>f, "A"*(cN*chL) 
	f.close()
	
	cmd = '/Users/sdh/apps/scwrl4/Scwrl4 -s temp.seq -i %s.bbq.pdb -o %s.bbq.scwrl.pdb' % (outname, outname)
	p = subprocess.Popen(cmd, shell=True, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
	p.wait()
	
	os.system('rm temp.seq hot.grp %s.pdb %s.bbq.pdb' % (outname, outname))
	
	return True
	
	# usuwanie skrajnych pozycji
	"""
	
	p = PDBParser(PERMISSIVE=1)

	structure_id = outname
	filename = outname + '.bbq.pdb'
	s = p.get_structure(structure_id, filename)
	
	for model in s:
		for chain in model:
			print chain
			for id in [chL*(i)+1 for i in range(cN)]:
				print id
				#chain.detach_child((' ', id, ' '))
	"""
				

def CCCPgenerate_periodicity(P, name, cN, chL, R0, ph1, cr, dph0, Zoff, Zoff_type, residues_per_turn, d, R1, justprint=False):
	"""
	P - periodicity of a bundle
	"""
	
	if P==residues_per_turn:
		print "P==p"
		sys.exit(-1)
	
	name = str(name)
	
	#cN = 2 # number of chains
	#chL = 56 # number of residues
	#R0 = 5.0 # superhelix radius
	#R1 = 2.26 # helix radius
	
	#ph1 = -9
	#cr = "[1]"
	#dph0 = "[180]"
	#Zoff = "[0]"
	
	
	# wartosci dla a-helisy
	#residues_per_turn = 3.62705 # residues / turn in a straight helix
	#d = 1.51 # axial rise per residue in a straight helix
	
	# in the science h-bond paper w1 is fixed (supp. data page 7)
	w1 = 360.0 / P # as in CCCP paper (not "p"!); page 1087
	
	#w1 = 70
	
	# option 1 
	#w0 = (- residues_per_turn / P + 1) / residues_per_turn * 360 # as in periodicity calculations in YadA paper
	#A = math.degrees(math.asin((R0 * math.radians(w0)) / d)) # pitch angle as in CCCP
	
	# option 2 (Lupas & Gruber 2005)
	# gives identical results
	
	delta_t = 2*math.pi * ( 1/residues_per_turn - 1/P )
	PITCH = (2*math.pi/delta_t) * math.sqrt(math.pow(d,2) - math.pow(R0*delta_t,2))
	
	A = math.degrees(math.atan(2*math.pi*R0 / PITCH))
	w0 = math.degrees((d * math.sin(math.radians(A))) / R0)
	
	w0 = numpy.around(w0, decimals=3)
	A = numpy.around(A, decimals=3)
	w1 = numpy.around(w1, decimals=3)
	PITCH = numpy.around(PITCH, decimals=3)

	if not justprint:
		CCCPgenerate(name, cN, chL, R0, R1, w0, w1, A, ph1, cr, dph0, Zoff, Zoff_type)
	
	return ' '.join([str(i[0])+"_"+i[1] for i in zip([name, P, d, R0, w0, A, R1, w1, ph1, dph0, Zoff, cr, PITCH, cN], ['name', 'P', 'd', 'R0', 'w0', 'A', 'R1', 'w1', 'ph1', 'dph0', 'Zoff', 'cr', 'PITCH', 'chains'])])


class FitParser():

	short_key = {'starting heptad position':'heptad',\
					  'R1 (A)':'R1',\
					  'R0 (A)':'R0',\
					  'rise per residue (A)': 'd',\
					  'alpha (rad)':'A',\
					  'w0 (rad/res)':'w0',\
					  'w1 (rad/res)':'w1',\
					  'error':'fit_rmsd',\
					  'pitch (A)':'PITCH',\
					  'ph1 (rad)':'ph1',\
					  'ph0 (rad)':'ph0',\
					  'absolute ap zoff (A)':'ab_ap_zoff',\
					  'dph0_p_ap (rad)':'dph0_p_ap',\
					  'dph0_2 (rad)':'dph0_2',\
					  'dph0_3 (rad)':'dph0_3',\
					  'dph0_4 (rad)':'dph0_4',\
					  'absolute zoff_4 (A)':'abs_zoff_4',\
					  'absolute zoff_2 (A)':'abs_zoff_2',\
					  'absolute zoff_3 (A)':'abs_zoff_3',\
					  'message':'message',\
					  'Z_aa for chain 2 (A)':'zaa2',\
					  'Z_aa for chain 3 (A)':'zaa3',\
					  'Z_aa for chain 4 (A)':'zaa4',\
					  }

	def __init__(self, mode=1):
		pass

	def parse(self, filename, todegs=False, prefix=""):

	

		# GENERAL (mode=0)

		#5 R0 (A) = 7.359923
		#6 R1 (A) = 2.259965
		#7 w0 (rad/res) = -0.042699
		#8 w1 (rad/res) = 1.795185
		#9 alpha (rad) = -0.209592
		#10 ph1 (rad) = -0.645646
		#11 dph0_2 (rad) = 1.570807
		#12 dph0_3 (rad) = -3.141592
		#13 dph0_4 (rad) = -1.570786
		#14 starting heptad position = b
		#15 pitch (A) = 217.397019
		#16 rise per residue (A) = 1.510448
		#17 register Z-offset for chain 2 (A) = 1.785075
		#18 Z_aa for chain 2 (A) = 1.694279
		#19 register Z-offset for chain 3 (A) = -0.000002
		#20 Z_aa for chain 3 (A) = -0.000002
		#21 register Z-offset for chain 4 (A) = 1.785069
		#22 Z_aa for chain 4 (A) = 1.694273
		#23 absolute zoff_2 (A) = 0.502193
		#24 absolute zoff_3 (A) = 0.000005
		#25 absolute zoff_4 (A) = 0.502185
		#26 error = 0.000563


		# SYMMETRIC (mode=1)

		#5 R0 (A) = 6.531083
		#6 R1 (A) = 2.265139
		#7 w0 (rad/res) = -0.044128
		#8 w1 (rad/res) = 1.794917
		#9 alpha (rad) = -0.191279
		#10 ph1 (rad) = -1.340822
		#11 dph0_p_ap (rad) = -1.851448
		#12 starting heptad position = b
		#13 pitch (A) = 211.912218
		#14 rise per residue (A) = 1.515943
		#15 register Z-offset for chain 2 (A) = 2.369028
		#16 Z_aa for chain 2 (A) = 2.365918
		#17 register Z-offset for chain 3 (A) = 2.369028
		#18 Z_aa for chain 3 (A) = 2.365918
		#19 register Z-offset for chain 4 (A) = 0.000000
		#20 Z_aa for chain 4 (A) = 0.000000
		#21 absolute ap zoff (A) = -3.507711
		#22 error = 0.700876


		f=open(filename)
		l=f.readline()
		while l.find('structure summary:')<>-1:
			l=f.readline()
		
		x = dict([x[:-1].split(' = ') for x in f.readlines()+[l]])
		
		for old_key in x.keys():
			data = x.pop(old_key)
			new_key = prefix+self.short_key[old_key]
			try:
				data = float(data)
				if todegs and new_key in ['w0', 'w1', 'A', 'ph1']:
					data = math.degrees(data)
			except:
				pass # string
				
			x[new_key] = data
		
		x['parfile'] = os.path.basename(filename)
		
		self.res = x
		return self.res

	def show(self):
	
		for param in self.CCCP_params:
			val = res[param]
			try:
				val = float(val)
				if param in ['w0', 'w1', 'alpha', 'dph0_2', 'dph0_3', 'dph0_4', 'ph1']:
					val = math.degrees(val)
				val = round(val,2)
			
			except:
				pass
			print '%s' % val ,
		print

	def renamewebCCCP(self, filename, phisteps, zsteps):
		assert phisteps==zsteps
		number = int(filename.split('.')[0])-1
	
		axial_rotation, axial_shift = number / zsteps, number % phisteps
	
		newname = "cccpweb_%s_%s.pdb" % (axial_shift, axial_rotation)
	
		return newname
	


def tetramer(ap=False):

	chains = 4
	delta_phi0 = "[90,-180,-90]" 
	pre_Zoff = "[{0},0,{0}]"

	R0_range = numpy.linspace(7.36, 7.36, 1) # CCCP default 7.36 

	if not ap:
		orient = "[1,1,1]" # parallel
	else:
		orient = "[0,1,0]" # antiparallel
		
	return chains, delta_phi0, pre_Zoff, orient, R0_range
	
if __name__ == "__main__":	


	#Zoff_type = "registerzoff"
	Zoff_type = "zoffaa"
	
	P = 3.5
	pmax = 16
	#Zoff_raw = 0
	length = 7*3
	
	f = open('bb_gen.log', 'w')
	
	for chains, delta_phi0, pre_Zoff, orient, R0_range in [tetramer(ap=True)]:
		R0 = R0_range[0]
		assert len(R0_range)==1
	
		for Zoff_raw in numpy.linspace(1, 4, 7):
	
			Zoff = pre_Zoff.format(Zoff_raw)

			pos=0
			for phi1 in numpy.linspace(-pmax, pmax, (pmax/2)+1):
				phi1 = phi1 - 6.5		

				name = "-".join([str(chains), str(Zoff_raw), str(pos)])

				res = CCCPgenerate_periodicity(P, name, chains, length, R0, "[%s]" % ','.join(itertools.repeat(str(phi1), chains)),\
													   orient, delta_phi0, Zoff, Zoff_type, 3.62705, ahelix_d, ahelix_R1, justprint=False)
				print >>f, res
				print res
				pos+=1



	#print CCCPgenerate_periodicity(3.631, 'prosta', 2, 36, 5, -9, '[1]', '[180]', '[0]', 'zoffaa', 3.63, 1.51, 2.26)
	#print CCCPgenerate_periodicity(3.5, 'krzywa', 2, 36, 5, -9, '[1]', '[180]', '[0]', 'zoffaa', 3.63, 1.51, 2.26)


	#print calcPER(3.55, tolerance=0.5)
	#sys.exit(-1)
	#print w0_P_to_p(-2.45, 3.5)
	#sys.exit(-1)

	#print inTheCore_new(44.331749962, 360/11.0, 11)
	
	#pass


	#print round(w0_to_P((np.degrees(float(sys.argv[1]))), 3.62125352),3)
	#print calcPER(4.45, tolerance=0.025)
	