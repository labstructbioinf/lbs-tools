import re, subprocess

def TMalign(pdb1, pdb2, bin="", args=""):

	assert bin!="", "Please provide a full path to TMalign binary"

	cmd = f'{bin} {pdb1} {pdb2} {args}'	

	res_rmsd, res_tmscore = None, None
	
	for l in str(subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE).stdout.read()).split('\\n'):
		if l.startswith('Aligned length'):
			temp = re.search('RMSD=(.*),', l)
			if temp: 
				res_rmsd = float(temp.group(1))
			else:
				print('TMalign warning. Could not parse RMSD!')
		elif l.find('if normalized by average length of chains')>-1:
			temp = re.search('TM-score=(.*)\(', l)	
			if temp:
				res_tmscore = float(temp.group(1))
			else:
				print('TMalign warning. Could not parse TM-Score!')	
			break
			
	return res_rmsd, res_tmscore
	