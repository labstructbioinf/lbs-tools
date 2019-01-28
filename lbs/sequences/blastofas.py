from Bio import SeqIO
from Bio.Blast import NCBIXML
from Bio.Alphabet import IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import numpy as np
import argparse
import re

def __correct_alignment(master_query, hsp):
    if '-' in hsp.query: # insertions
        indices = [match.start() for match in re.finditer('-', hsp.query)]
        tmp_sbjct = ''.join([hsp.sbjct[i] for i in range(0, len(hsp.sbjct)) if i not in indices])
        new_sbjct = "-"*(hsp.query_start-1) + tmp_sbjct + "-"*(len(master_query)-hsp.query_end)
    else: # no insertions
        starts = [(match.start(), match.end()) for match in re.finditer(str(hsp.query), master_query)]
        new_sbjct = "-"*starts[0][0] + hsp.sbjct+ "-"*(len(master_query)-starts[0][1])
    assert len(new_sbjct) == len(master_query)
    return new_sbjct

def blastofas(xmlfile, queryFasta, evalue=1e-3, iternr = -1, maxgaps=0.25):
	"""
	generates MSA based on XML PSI-BLAST results
	Important note: gaps in the query sequence are collapsed, thus
	hit sequences may not be contiguous
	"""
	query_seq = next(SeqIO.parse(queryFasta, "fasta"))
	psiblast_records = list(NCBIXML.parse(open(xmlfile)))[iternr]
	fasta = []
	fasta.append(query_seq)
	query_seq = str(query_seq.seq)
	for alignment in psiblast_records.alignments:
			for hsp in alignment.hsps:
				if hsp.expect <= float(evalue):
					corr_seq = __correct_alignment(query_seq, hsp)
					seq = Seq(corr_seq, IUPAC.protein)
					fasta.append(SeqRecord(seq, id=alignment.title, description=""))
	mat = np.array([list(seq.seq) for seq in fasta])
	mat = mat[:, ~(mat[0] == '-')]
	new_fasta = []
	for seq, fas in zip(mat, fasta):
		seq = ''.join(seq)
		if not (float(seq.count('-')) / mat.shape[1]) >= maxgaps:
			seq = SeqRecord(Seq(seq), description="", id=fas.id)
			new_fasta.append(seq)
	return new_fasta
    
def blastofasComplete(xmlfile, queryFasta, evalue=1e-3, iternr = -1, maxgaps=0.25):
	"""
	:param str xmlfile: input xml file from blast or psiblast
	:param str queryFasta: query sequence (fasta file)
	:param float evalue: max e-value of hsp
	:param int iternr: number of psiblast iteration to parse (-1 = last)
	:param float maxgaps: used only if mask=True. maximal fraction of gaps in a sequence
	"""

	queryseq = next(SeqIO.parse(open(queryFasta, 'r'), "fasta"))
	b_record = list(NCBIXML.parse(open(xmlfile)))[iternr]

	fasta = []

	masterquery = str(queryseq)
	fasta.append(queryseq)


	for alignment in b_record.alignments:
		for hsp in alignment.hsps:

			if hsp.expect <= float(evalue):
			
				temp = ''
				hsp.sbjct = "-"*(hsp.query_start-1) + hsp.sbjct
				mpos=0
				# dla kazdego znaku w sbjct
				for i in range(len(hsp.sbjct)): ## -hsp.query_start
					if i-hsp.query_start+1<0 or hsp.query[i-hsp.query_start+1]!='-': # nie ma insercji
					
						while fasta[0].seq[mpos]=='-':
							mpos = mpos + 1
							temp = temp + '-'      
						temp = temp + str(hsp.sbjct[i])

						mpos = mpos + 1    
					else: # jest insercja
					
						if fasta[0].seq[mpos]!='-':
							for f in range(len(fasta)): # trzeba dodac gap
					
								fasta[f].seq=fasta[f].seq[:mpos]+"-"+fasta[f].seq[mpos:]
						
						temp = temp + str(hsp.sbjct[i])
					
						mpos = mpos + 1
					
				simple_seq = Seq(temp, IUPAC.protein)

				fasta_rec = SeqRecord(simple_seq, id=alignment.title, description="")
				
				fasta.append(fasta_rec)


	# uzupelnienie sekwencji na N koncu

	for f in range(len(fasta)):
		if len(fasta[f])<len(fasta[0]):
			fasta[f].seq=fasta[f].seq+"-"*(len(fasta[0])-len(fasta[f]))

	return fasta    


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Blastofas')
    parser.add_argument('-xml',
		                help='Psiblast XML outfile',
		                required=True,
		                metavar='XML')
    parser.add_argument('-query',
		                help='Query sequence in fasta format',
		                metavar='QUERY',
		                required=True)
    parser.add_argument('-maxgaps',
		                help='Max percentage of gaps',
		                metavar='MAXGAPS',
		                default=0.25,
                        type=float)
    parser.add_argument('-iternr',
		                help='Number of iteration',
		                metavar='ITERNR',
		                default=-1,
                        type=int)
    parser.add_argument('-evalue',
		                help='Hit evalue threshold',
		                metavar='EVALUE',
		                default=1e-3,
                        type=float)
    parser.add_argument('-fastmode',
		                help='Use fast mode',
		                default=False,
		                action='store_true')                        
                                            
    args = parser.parse_args()
    
    if args.fastmode:
    	func = blastofas
    else:
    	func = blastofasComplete
    
    fasta = func(args.xml, args.query, evalue=args.evalue, iternr=args.iternr, maxgaps=args.maxgaps)
    for f in fasta:
        print(f.format('fasta'))

