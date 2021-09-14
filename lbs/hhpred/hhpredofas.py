

from csb.bio.io.hhpred import HHOutputParser

def generate_msa(hhr_file, queryseq, hitslist, maxevalue=1e-3, ident_cut=0.5, qcov_cut=0.5, eval_cut=1e-3):
	
	assert len(hitslist) > 0, 'provide at least one hit id in `hitlist`'

	fasta = [queryseq]

	for hit in HHOutputParser(alignments=True).parse_file(hhr_file):
	
		hit_id = f'{hit.id}_{hit.qstart}_{hit.qend}'
		if not hit_id in hitslist: continue
	
		query_cov = 1.*len(hit.alignment.subject.replace('-', ''))
		if hit.identity < ident_cut: continue
		if query_cov / len(queryseq) < qcov_cut: continue
		if hit.evalue > eval_cut: continue
		
		temp = ''
		mpos = 0
		sbjct = "-"*(hit.qstart-1) + hit.alignment.subject
		
		# for each aa in sbjct
		for i in range(len(sbjct)):
		
			# no insertion at this position
			if i - hit.qstart + 1 < 0 or hit.alignment.query[i - hit.qstart+ 1 ] != '-': 
				while fasta[0][mpos] == '-':
					mpos = mpos + 1
					temp = temp + '-'      
				temp = temp + str(sbjct[i])
				mpos = mpos + 1    
				
			# insertion present
			else: 
				if fasta[0][mpos] != '-':
					for f in range(len(fasta)): # we need to add a gap
						fasta[f] = fasta[f][:mpos] + "-" + fasta[f][mpos:]		
				temp = temp + str(sbjct[i])	
				mpos = mpos + 1
			
		fasta.append(temp)

	# fill gaps at the N terminus
	for f in range(len(fasta)):
		if len(fasta[f])<len(fasta[0]):
			fasta[f]=fasta[f]+ "-" * (len(fasta[0])-len(fasta[f]))

	return fasta    

	
