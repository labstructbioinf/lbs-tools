def parse_pdb_ranges(pdb_ranges):
	"""
		parses ecod domain pdb_ranges such as
			A:4,A:6-188
			A:474-512,A:603-675
			etc.
	"""
    
    for pdb_range in pdb_ranges.split(','):
        tmp = pdb_range.split(":")[1]

        if tmp[0] == '-':
            minus_start = True
            tmp = tmp[1:]
        else:
            minus_start = False

        if tmp.find('--')>0:
            tmp = tmp.replace('--', '-')
            minus_end = True
        else:
            minus_end = False
            
        if tmp.find('-') != -1:
        	ecod_start, ecod_end = tmp.split("-")
		else:
			ecod_start = ecod_end = tmp

        if minus_start:
            ecod_start = '-' + ecod_start
            
        if minus_end:
            ecod_end = '-' + ecod_end
            
        yield ecod_start, ecod_end
        
        

