import subprocess

def run_segmasker(sequence, window=12, locut=2.2, hicut=2.5, debug=False,
				  segmasker_path = "/opt/apps/ncbi-blast+/bin/segmasker"):
				  
    """
    Run SEG algorithm using segmasker and return a list of ranges.
    
    It scans a sequence using a window of a given length (W
	parameter) and checks whether the first threshold (locut) is met. If it is then the
	fragment is extended to the left and to the right until it reaches the second threshold
	(hicut). 

    Args:
    - sequence (str): The input sequence to be analyzed.
    - window (int): SEG window size.
    - locut (float): SEG locut threshold.
    - hicut (float): SEG hicut threshold.

    Returns:
    - list: A list of tuples representing the ranges of low complexity regions.
    """
    
    # Run SEG algorithm with specified parameters
    cmd = [
        segmasker_path,
        "-window", str(window),
        "-locut", str(locut),
        "-hicut", str(hicut),
       # "-outfmt", "fasta"
    ]
        
    if debug:
    	print(' '.join(cmd))
        
    # Run the command with the sequence input
    process = subprocess.Popen(cmd, stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    stdout, stderr = process.communicate(f'>Sequence\n{sequence}')
        
    if debug:
    	print(stdout)
    
    masked_seq = sequence
    lcr_residues = 0
    
    # Parse the output to extract low complexity regions
    ranges = []
    for line in stdout.split('\n'):
        if not line.startswith('>') and line!='':
            start, end = line.split(' - ')
            start, end = int(start), int(end)
            ranges.append((start, end))
            masked_seq = masked_seq[:start] + masked_seq[start:end+1].lower() + masked_seq[end+1:]
            lcr_residues += end-start+1
                        
    assert len(sequence) == len(masked_seq)
                        
    return ranges, masked_seq, lcr_residues / len(sequence)
    
if __name__ == '__main__':

	example = 'AGVEEVAASGSHLNGDLDPDDREEGAASTAEEAAKKKRRKKKKSKGPSAAGEQEPDKESGASVDEVARQLERSALEDKE'

	ranges, masked_ex, lcr_res = run_segmasker(example, debug=False)

	for r in ranges:
		print(r)
		
	print(example)
	print(masked_ex)
		
	