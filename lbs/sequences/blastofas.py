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
    args = parser.parse_args()
    fasta = blastofas(args.xml, args.query, evalue=args.evalue, iternr=args.iternr, maxgaps=args.maxgaps)
    for f in fasta:
        print(f.seq)

