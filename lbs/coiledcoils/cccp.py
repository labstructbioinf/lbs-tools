# This file uses PEP8 conventions.
# Variables that also appear in the CCCP Generator, but under a different name, are listed below.

""" This script uses the CCCP generator developed and described by G. Grigoryan and W. F. DeGrado in
"Probing Designability via a Generalized Model of Helical Bundle Geometry", J. Mol. Biol., 405(4): 1079-1100 (2011).
See also http://grigoryanlab.org/cccp/.

Description of variables:

	Variables that also appear in the CCCP generator:

	num_chains 	 - the number of chains.
				   In the CCCP generator, it is alternately called cN and chains (k is used as a chain index)
	chain_length - the length of each chain (in amino acid residues). Must be an integer.
	r0   - the superhelical radius, called r0 in the CCCP Generator (R0 in the paper). Given in Angstrom.
	r1   - the helical radius, called r1 in the CCCP Generator (R1 in the paper). Given in Angstrom.
	w0   - the superhelical frequency, in degrees per residue.
	w1   - the helical frequency, in degrees per residue.
	a    - the pitch angle (called A in the CCCP generator, alpha_p in the paper), in degrees.
	ph1  - the starting helical phase, in degrees. Can be passed as a scalar if its the same for all helices.
		   Also referred to as the helical phase offset.
	dph0 - the superhelical phase offset, passed as a vector of length K = num_chains. 
		   For any chain k, dpho[k - 1] indicates the phase offset from the first chain, in degrees.
	c_dir  - The orientation of all subsequent chains relative to the first chain.
			 c_dir is a vector of length K = num_chains.
			 c_dir has (alternating) entries 1 and 0, 1 indicating parallel and 0 indicating antiparallel orientation.
			 The first entry must be 1 (the first chain is parallel to itself).
			 Example: for any chain k > 1, if c_dir(k - 1) = 1, its orientation is parallel to the first chain.
	z_off  - Alternately called zoff and Zoff in the CCCP generator, refers to the absolute vertical offset of a chain k
			 relative to the first chain.
			 z_off is a vector of length K = num_chains.
			 For any chain k, z_off[k - 1] indicates the offset in the z direction relative to the first chain.
			 Since absolute values might not be useful for fitting, the parameter z_off_type can be used to modify z_off.
	z_type - Also referred to as varargin in the CCCP generator, used to modify z_off. Passed as a string.
			 Description taken / paraphrased from the CCCP generator:
			 zoffaa - z_off referes to the vertical offset between 'a' positions on opposing chains
			 registerzoff - z_off will refer to points on the Crick curves that point directly into the interface
			 apNNzoff - the Z offset is interpreted as between the N termini of chains 1 and k if k runs parallel to 1,
			 and the N terminus of chain 1 and the C terminus of chain k if k is anti-parallel to 1.

	Other variables defined for this script:
	
	rep_len - repeat length

	Where variables are renamed or reassigned in the functions and class defined below, clear indications are given 
	in the function, class or method docstrings. Where no indication is given, or where the variable is simply
	introduced, its meaning should match the one defined above.

The script also calls the BBQ (Backbone Building from Quadrulaterals) developed by Dominik Gront, Sebastian Kmiecik 
and Andrzej Kolinski and described in "Backbone building from quadrilaterals: A fast and accurate algorithm for protein 
backbone reconstruction from alpha carbon coordinates", Journal of Computational Chemistry, March 2007. 
See also: https://doi.org/10.1002/jcc.20624 
"""

import tempfile
import sys # use sys.exit(-1) to debug
import math
import numpy
import subprocess
import os
import os.path
import itertools

from sympy import nsimplify
from Bio.PDB.PDBParser import PDBParser

ahelix_p = 3.62705
ahelix_d = 1.51
ahelix_R1 = 2.26

#4hpq
pihelix_p = 4.3985
# Wychodzi 1.1689 (helanal), spore stdev - 0.3423
pihelix_d = 1.1689
	

def in_the_core(ph1, w1, rep_len):
	""" This function calculates all residue angles in the given repeat length and outputs the smallest one.

	It takes three inputs:
	ph1 - starting helical phase (helical phase offset)
	w1  - helical frequency in degrees per residue, examples: w1 = 102.8 for 7/2; w1 = 98.18181818181819 for 11/3;
	rep_len - repeat length, passed as a float.
	"""

	def adj(ang: float) -> float:
		""" This function checks whether an angle exceeds 180 degrees, and converts it to a negative angle if it does.
		The input angle must be in degrees.

		Examples: If the angle is 270 degrees, it will be converted to -90 degrees. 
				  If the angle is 90 degrees, it will not be changed. 
				  If the angle is 450 degrees, it will be converted to 90 degrees. 
		"""

		ang += 180
	
		if ang > 360:
			c = int(ang/360) # Calculate by how many full turns the angle exceeds 360 degrees
			return ang - (360*c) - 180 
		else:
			return ang - 180
	
	full_repeat = [adj(ph1 + (w1*i)) for i in range(rep_len)] # Create vector with residue angles in degrees

	return sorted(full_repeat, key = lambda x:abs(x))[0] # Return smallest residue angle in degrees

	
def in_the_core_pandas(row):
	""" This function calls in_the_core (defined above), using values specified in a pandas data frame as inputs. """
	return in_the_core(row['ph1'], row['w1'], row['rep_len'])

def calcPER(P, tolerance = 0.01):
	# TODO: What does this function do?
	f0 = nsimplify(P, tolerance = tolerance)
	return f0.p, f0.q

def w0_P_to_p(w0, P):
	return P / (1 + P*w0 / 360)

def w0_to_P(w0, residues_per_turn):
	return residues_per_turn / (1 - (residues_per_turn*w0 / 360))
	
def P_to_w0(P, residues_per_turn):
	return  (- residues_per_turn / P + 1) / residues_per_turn*360

def calcphi(heptad, phi1):

	heporder = 'dgcfbea'
	heptadlen   = 360/7.0

	distto_f = 3 - heporder.index(heptad)   
	phi = phi1 - heptadlen*distto_f

	return phi+9 # -26 -- 26 range

def make_bb(f, outname, chain_length, chain_names  = ['A', 'B', 'C', 'D', 'E']):
	""" This function takes a file with coordinates XYZ as its input and writes a PDB file. 
	The input file is generated using the CCCP Generator ("generateCrickBB.m").  
	It is called in the cccp_generate wrapper below. 

	This function takes four inputs:

	f - the file containing the XYZ coordinates
	outname - the name of the outputted pdb file
	chain_length - the length of the individual chains
	chain_names  - a list of letters to assign to individual chains
	"""
	
	fpdb = open(outname + '.pdb', 'w') # Initialize PDB file to be written

	atom_id = 1 # Start with first CA as first atom
	res_id  = 1 # Start with first residue

	for coords in f:

		coords = numpy.round(coords, 3) # round coordinates to three digits
		fpdb.write('ATOM  ' + '{0: >5}'.format(str(atom_id)) + '{0: >4}'.format('CA') + ' {0: >4}'.format('GLY') + ' ' + chain_names[int((res_id-1) / chain_length)] + 
				   '{0: >4}'.format(str(res_id)) + '    ' + '{0: >8}'.format(coords[0]) + '{0: >8}'.format(coords[1]) + '{0: >8}'.format(coords[2]) + 
				   '{0: >6}'.format('1.00') + '{0: >6}'.format('0.00') + '\n') # Write PDB, use placeholders for consistent formatting
		res_id  += 1 # Move to the next residue
		atom_id += 4 # Move to the next CA

	fpdb.close()

def cccp_generate(outname, num_chains, chain_length, r0, r1, w0, w1, a, ph1, c_dir, dph0, z_off, z_type):
	""" This function takes the specified Crick parameters as its input and acts as a wrapper that:
	
	1. Calls the octave script (generateCrickBB.m, see module docstring for citation and doi)
	2. Writes the output coordinates (backbone trace) to a PDB file
	3. Call BBQ (Backbone Building from Quadrulaterals) to convert the backbone trace to poly-Alanin

	With the exception of the name of the output files (outname), all parameters are defined in the module docstring (see above).
	"""

	# Write parameters to temporary file to pass them to Octave
	tf = tempfile.NamedTemporaryFile(delete=False)
	tf.write("generateCrickBB(%s, %s, %s, %s, %s, %s, %s, %s, %s, %s, %s, struct(\'%s\', 1))" % (num_chains, chain_length, r0, r1, w0, w1, a, ph1, c_dir, dph0, z_off, z_type))
	tf.close()
	
	# Call Octave, get generated backbone trace
	cmd = "octave-cli -q %s" % tf.name # ... and pass parameters to generate backbone trace
	
	current_working_directory = os.getcwd() # re-define command for better readability (access current working directory)
	os.chdir('/Users/sdh/apps/cccp') # Change directory to CCCP app location
	
	p = subprocess.Popen(cmd, shell = True, stderr = subprocess.PIPE, stdout = subprocess.PIPE)
	p.wait()
	stdout, stderr = p.communicate() # Get the output (the generated backbone trace)	
	
	if stderr != '':
		print('something went wrong...')
		sys.exit(-1)
	
	os.chdir(current_working_directory) 
	
	# Process the backbone trace XYZ
	coords = numpy.array([i.split() for i in stdout.split('\n')[2:-2]], dtype = numpy.dtype(float) ) # Write to numpy array
	make_bb(coords, outname, chain_length) # Create PDB file
	
	# Convert backbone trace to poly-Alanin using BBQ
	cmd = 'java -classpath /Users/sdh/apps/BBQ/:/Users/sdh/apps/BBQ/jbcl.jar BBQ -d=q_50_xyz.dat -r=./' + str(outname) + '.pdb > ./' + str(outname) + '.bbq.pdb'
	p = subprocess.Popen(cmd, shell = True, stderr = subprocess.PIPE, stdout = subprocess.PIPE)
	p.wait()
	
	f = open('temp.seq', 'w') # Write Alanin one letter code to amino acid sequence file
	print >>f, "A"*(num_chains*chain_length) 
	f.close()
	
	cmd = '/Users/sdh/apps/scwrl4/Scwrl4 -s temp.seq -i %s.bbq.pdb -o %s.bbq.scwrl.pdb' % (outname, outname)
	p = subprocess.Popen(cmd, shell = True, stderr  =subprocess.PIPE, stdout = subprocess.PIPE)
	p.wait()
	
	os.system('rm temp.seq hot.grp %s.pdb %s.bbq.pdb' % (outname, outname))
	
	return True
	
	# removing extreme positions
	#"""
	
	#p = PDBParser(PERMISSIVE=1)

	#structure_id = outname
	#filename = outname + '.bbq.pdb'
	#s = p.get_structure(structure_id, filename)
	
	#for model in s:
		#for chain in model:
			#print chain
			#for id in [chain_length*(i)+1 for i in range(num_chains)]:
				#print id
				#chain.detach_child((' ', id, ' '))
	#"""
				

def cccp_generate_periodicity(P, name, num_chains, chain_length, r0, ph1, c_dir, dph0, z_off, z_type, residues_per_turn, d, r1, justprint=False):
	"""
	P - periodicity of a bundle
	"""
	
	if P == residues_per_turn:
		print("P == p")
		sys.exit(-1)
	
	name = str(name)
	
	# values ​​for a-helix
	#residues_per_turn = 3.62705 # residues / turn in a straight helix
	#d = 1.51 # axial rise per residue in a straight helix
	
	# in the science h-bond paper w1 is fixed (supp. data page 7)
	w1 = 360.0 / P # as in CCCP paper (not "p"!); page 1087
	
	# option 1 
	#w0 = (- residues_per_turn / P + 1) / residues_per_turn * 360 # as in periodicity calculations in YadA paper
	#a = math.degrees(math.asin((r0 * math.radians(w0)) / d)) # pitch angle as in CCCP
	
	# option 2 (Lupas & Gruber 2005)
	# gives identical results
	
	delta_t = 2*math.pi * ( 1/residues_per_turn - 1/P )
	pitch = (2*math.pi/delta_t) * math.sqrt(math.pow(d,2) - math.pow(r0*delta_t,2))
	
	a = math.degrees(math.atan(2*math.pi*r0 / pitch))
	w0 = math.degrees((d * math.sin(math.radians(a))) / r0)
	
	w0 = numpy.around(w0, decimals=3)
	a = numpy.around(a, decimals=3)
	w1 = numpy.around(w1, decimals=3)
	pitch = numpy.around(pitch, decimals=3)

	if not justprint:
		CCCPgenerate(name, num_chains, chain_length, r0, r1, w0, w1, a, ph1, c_dir, dph0, z_off, z_type)
	
	return ' '.join([str(i[0])+"_"+i[1] for i in zip([name, P, d, r0, w0, a, r1, w1, ph1, dph0, z_off, c_dir, pitch, num_chains], ['name', 'P', 'd', 'r0', 'w0', 'a', 'r1', 'w1', 'ph1', 'dph0', 'z_off', 'c_dir', 'pitch', 'chains'])])


class FitParser():


	short_key = {'starting heptad position':'heptad',
				 'r1 (a)':'r1',
				 'r0 (a)':'r0',
				 'rise per residue (a)': 'd',
				 'alpha (rad)':'a',
				 'w0 (rad/res)':'w0',
				 'w1 (rad/res)':'w1',
				 'error':'fit_rmsd',
				 'pitch (a)':'PITCH',
				 'ph1 (rad)':'ph1',
				 'ph0 (rad)':'ph0',
				 'absolute ap z_off (a)':'ab_ap_z_off',
				 'dph0_p_ap (rad)':'dph0_p_ap',
				 'dph0_2 (rad)':'dph0_2',
				 'dph0_3 (rad)':'dph0_3',
				 'dph0_4 (rad)':'dph0_4',
				 'absolute z_off_4 (a)':'abs_z_off_4',
				 'absolute z_off_2 (a)':'abs_z_off_2',
				 'absolute z_off_3 (a)':'abs_z_off_3',
				 'message':'message',
				 'Z_aa for chain 2 (a)':'zaa2',
				 'Z_aa for chain 3 (a)':'zaa3',
				 'Z_aa for chain 4 (a)':'zaa4',
				 }

	def __init__(self, mode=1):
		pass

	def parse(self, filename, todegs = False, prefix = ""):


		# GENERAL (mode=0)

		#5 r0 (a) = 7.359923
		#6 r1 (a) = 2.259965
		#7 w0 (rad/res) = -0.042699
		#8 w1 (rad/res) = 1.795185
		#9 alpha (rad) = -0.209592
		#10 ph1 (rad) = -0.645646
		#11 dph0_2 (rad) = 1.570807
		#12 dph0_3 (rad) = -3.141592
		#13 dph0_4 (rad) = -1.570786
		#14 starting heptad position = b
		#15 pitch (a) = 217.397019
		#16 rise per residue (a) = 1.510448
		#17 register z_off for chain 2 (a) = 1.785075
		#18 Z_aa for chain 2 (a) = 1.694279
		#19 register z_off for chain 3 (a) = -0.000002
		#20 Z_aa for chain 3 (a) = -0.000002
		#21 register z_off for chain 4 (a) = 1.785069
		#22 Z_aa for chain 4 (a) = 1.694273
		#23 absolute z_off_2 (a) = 0.502193
		#24 absolute z_off_3 (a) = 0.000005
		#25 absolute z_off_4 (a) = 0.502185
		#26 error = 0.000563


		# SYMMETRIC (mode=1)

		#5 r0 (a) = 6.531083
		#6 r1 (a) = 2.265139
		#7 w0 (rad/res) = -0.044128
		#8 w1 (rad/res) = 1.794917
		#9 alpha (rad) = -0.191279
		#10 ph1 (rad) = -1.340822
		#11 dph0_p_ap (rad) = -1.851448
		#12 starting heptad position = b
		#13 pitch (a) = 211.912218
		#14 rise per residue (a) = 1.515943
		#15 register z_off for chain 2 (a) = 2.369028
		#16 Z_aa for chain 2 (a) = 2.365918
		#17 register z_off for chain 3 (a) = 2.369028
		#18 Z_aa for chain 3 (a) = 2.365918
		#19 register z_off for chain 4 (a) = 0.000000
		#20 Z_aa for chain 4 (a) = 0.000000
		#21 absolute ap z_off (a) = -3.507711
		#22 error = 0.700876


		f = open(filename)
		l = f.readline()
		while l.find('structure summary:') != -1:
			l=f.readline()
		
		x = dict([x[:-1].split(' = ') for x in f.readlines() + [l]])
		
		for old_key in x.keys():
			data = x.pop(old_key)
			new_key = prefix + self.short_key[old_key]
			try:
				data = float(data)
				if todegs and new_key in ['w0', 'w1', 'a', 'ph1']:
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
			print('%s' % val),

	def renameweb_cccp(self, filename, phisteps, zsteps):
		assert phisteps == zsteps
		number = int(filename.split('.')[0])-1
	
		axial_rotation, axial_shift = number / zsteps, number % phisteps
	
		newname = "cccpweb_%s_%s.pdb" % (axial_shift, axial_rotation)
	
		return newname
	


def tetramer(ap=False):

	num_chains = 4
	delta_phi0 = "[90,-180,-90]" 
	pre_z_off = "[{0},0,{0}]"

	R0_range = numpy.linspace(7.36, 7.36, 1) # CCCP default 7.36 

	if not ap:
		orient = "[1,1,1]" # parallel
	else:
		orient = "[0,1,0]" # antiparallel
		
	return num_chains, delta_phi0, pre_z_off, orient, R0_range
	
if __name__ == "__main__":	


	#Zoff_type = "registerzoff"
	z_type = "zoffaa"
	
	P = 3.5
	pmax = 16
	#Zoff_raw = 0
	length = 7*3
	
	f = open('bb_gen.log', 'w')
	
	for num_chains, delta_phi0, pre_z_off, orient, R0_range in [tetramer(ap=True)]:
		r0 = R0_range[0]
		assert len(R0_range) == 1
	
		for z_off_raw in numpy.linspace(1, 4, 7):
	
			z_off = pre_z_off.format(z_off_raw)

			pos = 0
			for phi1 in numpy.linspace(-pmax, pmax, pmax/2 + 1):
				phi1 = phi1 - 6.5		

				name = "-".join([str(num_chains), str(z_off_raw), str(pos)])

				res = cccp_generate_periodicity(P, name, num_chains, length, r0, "[%s]" % ','.join(itertools.repeat(str(phi1), num_chains)),\
													   orient, delta_phi0, z_off, z_off_type, 3.62705, ahelix_d, ahelix_R1, justprint=False)
				print >>f, res
				print(res)
				pos += 1



	#print CCCPgenerate_periodicity(3.631, 'prosta', 2, 36, 5, -9, '[1]', '[180]', '[0]', 'zoffaa', 3.63, 1.51, 2.26)
	#print CCCPgenerate_periodicity(3.5, 'krzywa', 2, 36, 5, -9, '[1]', '[180]', '[0]', 'zoffaa', 3.63, 1.51, 2.26)


	#print calcPER(3.55, tolerance=0.5)
	#sys.exit(-1)
	#print w0_P_to_p(-2.45, 3.5)
	#sys.exit(-1)

	#print in_the_core(44.331749962, 360/11.0, 11)
	
	#pass


	#print round(w0_to_P((np.degrees(float(sys.argv[1]))), 3.62125352),3)
	#print calcPER(4.45, tolerance=0.025)